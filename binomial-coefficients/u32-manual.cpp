// $ time echo 4294967294 123456789 | ./a.out 
// 1985794154
// 
// real 0m0.429s

#pragma GCC optimize ("O3")
#pragma GCC target ("sse4.1")
#pragma GCC optimize ("unroll-loops")
#include <bits/stdc++.h>
#include <immintrin.h>
using namespace std;

typedef __m128i M;

#define rep(i,a,b) for (int i = (a); i < (b); ++i)
typedef long long ll;
typedef uint32_t u32;

u32 modpow(u32 a, ll e) {
	if (e == 0) return 1;
	u32 x = modpow(a * a, e >> 1);
	return e & 1 ? x * a : x;
}

u32 solve(ll N, ll K) {
	ll trailingZeroes = 0;
	map<u32, int> ivs;
	int accMult = 0;
	rep(i,0,3) {
		ll x = (i == 0 ? N : i == 1 ? K : N-K);
		int mult = (i == 0 ? 1 : -1);
		while (x > 0) {
			// Include the product (1 * 3 * ... * (x % 2^32)) ^ mult in the answer.
			ivs[(u32)x + 1] += mult;
			accMult += mult;
			x /= 2;
			trailingZeroes += x * mult;
		}
	}

	u32 cur = 0, res = 1, resdiv = 1;
	for (auto pa : ivs) {
		u32 lim = pa.first, ilim = lim / 2;
		// The odd numbers in the range [last lim, lim) get included 'accMult'
		// times in the answer -- 'pa.second' times for the interval ending at
		// 'lim', and also for all larger intervals. We divide by 2 to avoid
		// potential for overflow.
		const int PAR = 8;
		M madd = _mm_set1_epi32(PAR * 8);
		M msubprod[PAR], mcur2[PAR];
		rep(i,0,PAR) {
			msubprod[i] = _mm_set1_epi32(1);
			mcur2[i] = _mm_setr_epi32(
				2*cur + 8*i + 1,
				2*cur + 8*i + 3,
				2*cur + 8*i + 5,
				2*cur + 8*i + 7);
		}
		for (; cur + PAR * 4 < ilim; cur += PAR * 4) {
			rep(i,0,PAR) {
				msubprod[i] = _mm_mullo_epi32(msubprod[i], mcur2[i]);
				mcur2[i] = _mm_add_epi32(mcur2[i], madd);
			}
		}
		u32 prod = 1;
		union {
			M mprod;
			u32 parts[4];
		} u;
		u.mprod = _mm_set1_epi32(1);
		rep(i,0,PAR) u.mprod = _mm_mullo_epi32(u.mprod, msubprod[i]);
		rep(i,0,4) prod *= u.parts[i];
		for (; cur < ilim; cur++) {
			prod *= cur * 2 + 1;
		}
		if (accMult > 0) res *= modpow(prod, accMult);
		else resdiv *= modpow(prod, -accMult);
		accMult -= pa.second;
	}

	res *= modpow(2, trailingZeroes);
	res *= modpow(resdiv, (1LL << 31) - 1);
	return res;
}

int main() {
	ll N, K;
	cin >> N >> K;
	cout << solve(N, K) << endl;
}
