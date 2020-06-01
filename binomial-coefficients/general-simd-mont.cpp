#pragma GCC optimize ("O3")
#pragma GCC target ("avx2")
#include <bits/stdc++.h>
#include <immintrin.h>
using namespace std;

#define rep(i,a,b) for (int i = (a); i < (b); ++i)
#define all(x) x.begin(), x.end()
typedef long long ll;
typedef uint32_t u32;
typedef uint64_t u64;

ll euclid(ll a, ll b, ll &x, ll &y) {
	if (b) { ll d = euclid(b, a % b, y, x);
		return y -= a/b * x, d; }
	return x = 1, y = 0, a;
}
ll crt(ll a, ll m, ll b, ll n) {
	if (n > m) swap(a, b), swap(m, n);
	ll x, y, g = euclid(m, n, x, y);
	assert((a - b) % g == 0); // else no solution
	x = (b - a) % n * x % n / g * m + a;
	return x < 0 ? x + m*n/g : x;
}

template<class F>
void factor(ll n, F f) {
	for (ll p = 2; p*p <= n; p++) {
		int a = 0;
		while (n % p == 0) n /= p, a++;
		if (a > 0) f(p, a);
	}
	if (n > 1) f(n, 1);
}

u32 invp2(u32 x) { // x^-1 mod 2^32
	unsigned xinv = 1;
	rep(i,0,5) xinv = xinv * (2 - x * xinv);
	return xinv;
}

ll modpow(ll a, ll e, ll m) {
	if (e == 0) return 1;
	ll x = modpow(a * a % m, e >> 1, m);
	return e & 1 ? x * a % m : x;
}

struct Mont {
	u32 m, npr;

	Mont(u32 mod) : m(mod), npr(-invp2(mod)) {}

	u32 redc(u64 a) {
		u32 b = (u32)a * npr;
		u64 c = a + (u64)b * m;
		return (u32)(c >> 32);
	}
};

// N choose K modulo p^a
ll solve(ll N, ll K, int p, int a) {
	int mod = 1;
	rep(i,0,a) mod *= p;

	ll trailingZeroes = 0;
	map<int, int> ivs;
	int accMult = 0;
	ll res = 1, resdiv = 1;
	rep(i,0,3) {
		ll x = (i == 0 ? N : i == 1 ? K : N-K);
		int mult = (i == 0 ? 1 : -1);
		while (x > 0) {
			// Include the product (1 * 2 * ... * x) ^ mult in the answer,
			// with the numbers divisible by p excluded.
			ll lim = x + 1;
			ivs[(int)(lim % mod)] += mult;
			if (lim / mod % 2 == 1 && (mod == 4 || p > 2))
				res *= -1;
			accMult += mult;
			x /= p;
			trailingZeroes += x * mult;
		}
	}

	u32 pinv = invp2(p);
	u32 plim = 0xFFFFFFFF / p;

	Mont mont(mod);

	int cur = 0, nextp = 0;
	for (auto pa : ivs) {
		int lim = pa.first;
		// The numbers in the range [cur, lim) that aren't divisible by 'p' get
		// included 'accMult' times in the answer -- 'pa.second' times for the
		// interval ending at 'lim', and also for all larger intervals.
		ll prod = 1;
		if (p == 2) {
			int ilim = lim / 2;
			const int PAR = 8;
			u32 subprod[PAR] = {1,1,1,1,1,1,1,1};
			u32 cur2 = cur * 2 + 1;
			for (; cur + PAR < ilim; cur += PAR, cur2 += PAR*2) {
				rep(i,0,PAR) subprod[i] *= cur2 + 2*i;
			}
			u32 uprod = 1;
			rep(i,0,PAR) uprod *= subprod[i];
			for (; cur < ilim; cur++) {
				uprod *= cur * 2 + 1;
			}
			prod = uprod % mod;
		} else {
			const int PAR = 16;
			typedef __m256i mi;
			mi msubprod[PAR];
			rep(i,0,PAR) msubprod[i] = _mm256_set1_epi64x(1);
			mi mcur = _mm256_setr_epi64x(cur, cur + 1, cur + 2, cur + 3);
			mi mnpr = _mm256_set1_epi64x(mont.npr);
			mi mm = _mm256_set1_epi64x(mod);
			int iters = lim - cur;

			while (nextp < cur) nextp += p;
			mi mnextp = _mm256_set1_epi64x(nextp - 1);
			mi mp = _mm256_set1_epi64x(p);

			int maxStep = (p == 3 ? 6 : 5) * PAR;
			while (cur + maxStep <= lim) {
				rep(j,0,PAR) {
					mi mnext = _mm256_add_epi64(mcur, _mm256_set1_epi64x(4));
					cur += 4;

					if (nextp <= cur) {
						// (save half a cycle by using a 32-bit compare)
						mcur = _mm256_sub_epi32(mcur, _mm256_cmpgt_epi32(mcur, mnextp));
						mnext = _mm256_add_epi64(mnext, _mm256_set1_epi64x(1));
						mnextp = _mm256_add_epi64(mnextp, mp);
						nextp += p;
						cur++;
						iters--;
						if (p == 3) {
							mcur = _mm256_sub_epi32(mcur, _mm256_cmpgt_epi32(mcur, mnextp));
							mnext = _mm256_add_epi64(mnext, _mm256_set1_epi64x(1));
							mnextp = _mm256_add_epi64(mnextp, mp);
							nextp += p;
							cur++;
							iters--;
						}
					}

					mi ma = _mm256_mul_epu32(msubprod[j], mcur);
					mi mb = _mm256_mul_epu32(ma, mnpr);
					mi mc = _mm256_add_epi64(ma, _mm256_mul_epu32(mb, mm));
					msubprod[j] = _mm256_srli_epi64(mc, 32);
					mcur = mnext;
				}
			}

			rep(i,0,PAR) {
				union { u64 i[4]; mi m; } u;
				u.m = msubprod[i];
				rep(j,0,4)
					prod = mont.redc((u64) prod * u.i[j]);
			}

			for (; cur < lim; cur++) {
				u32 cur2 = (u32)cur * pinv > plim ? cur : 1;
				prod = mont.redc((u64) prod * cur2);
			}
			prod = prod * modpow(2, 32LL * (iters + 4 * PAR), mod) % mod;
		}
		if (accMult > 0) res = res * modpow(prod, accMult, mod) % mod;
		else resdiv = resdiv * modpow(prod, -accMult, mod) % mod;
		accMult -= pa.second;
	}

	res = res * modpow(p, trailingZeroes, mod) % mod;
	res = res * modpow(resdiv, mod - mod/p - 1, mod) % mod;
	return res;
}

// N choose K modulo M
ll solve(ll N, ll K, ll M) {
	ll res = 0, prod = 1;
	factor(M, [&](ll p, int a) {
		ll r = solve(N, K, (int)p, a), m = 1;
		rep(i,0,a) m *= p;
		res = crt(res, prod, r, m);
		prod *= m;
	});
	return res;
}

void test() {
	vector<vector<int>> binom(100, vector<int>(100, -1));
	rep(m,1,100) rep(n,0,100) rep(k,0,n+1) {
		if (k == 0 || k == n) binom[n][k] = 1 % m;
		else binom[n][k] = (binom[n-1][k-1] + binom[n-1][k]) % m;
		ll r = solve(n, k, m), r2 = binom[n][k];
		if (r != r2) {
			cerr << n << ' ' << k << ' ' << m << ' ' << r << ' ' << r2 << endl;
			abort();
		}
	}
	exit(0);
}

int main(int argc, char** argv) {
	if (argc > 1 && argv[1] == string("--test")) test();
	ll N, K, M;
	cin >> N >> K >> M;
	cout << solve(N, K, M) << endl;
}
