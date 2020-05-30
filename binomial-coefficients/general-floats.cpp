#pragma GCC optimize ("O3")
#include <bits/stdc++.h>
using namespace std;

#define rep(i,a,b) for (int i = (a); i < (b); ++i)
#define all(x) x.begin(), x.end()
typedef long long ll;
typedef uint32_t u32;

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

ll modpow(ll a, ll e, ll m) {
	if (e == 0) return 1;
	ll x = modpow(a * a % m, e >> 1, m);
	return e & 1 ? x * a : x;
}

int modmul(int a, int b, int m) {
	ll ab = (ll)a * b;
	return (int)(ab - (u32)(int)(1.0 / m * ab) * m);
}

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

	u32 pinv = (u32)modpow(p, (1LL << 31) - 1, 1LL << 32);
	u32 plim = 0xFFFFFFFF / p;

	int cur = 1;
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
			const int PAR = 8;
			int subprod[PAR];
			rep(i,0,PAR) subprod[i] = 1;
			while (cur + PAR <= lim) {
				rep(j,0,PAR) {
					int cur2 = (u32)cur * pinv > plim ? cur : 1;
					subprod[j] = modmul(subprod[j], cur2, mod);
					cur++;
				}
			}
			rep(i,0,PAR) prod = modmul((int) prod, subprod[i], mod);
			for (; cur < lim; cur++) {
				int cur2 = (u32)cur * pinv > plim ? cur : 1;
				prod = modmul((int)prod, cur2, mod);
			}
			prod %= mod;
			if (prod < 0) prod += mod;
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
