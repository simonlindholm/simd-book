#include <bits/stdc++.h>
#include <immintrin.h>
using namespace std;

#define rep(i, from, to) for (int i = from; i < (to); ++i)
#define trav(a, x) for (auto& a : x)
#define all(x) x.begin(), x.end()
#define sz(x) (int)(x).size()
typedef long long ll;
typedef pair<int, int> pii;
typedef vector<int> vi;

ll modpow(ll a, ll e, ll mod) {
	if (e == 0) return 1 % mod;
	ll x = modpow(a, e >> 1, mod);
	x = x * x % mod;
	if (e & 1)
		x = x * a % mod;
	return x;
}

unsigned modmul(unsigned a, unsigned b, unsigned M) {
	int ret = a * b - M * unsigned(1.0 / M * a * b);
	return ret + M * (ret < 0) - M * (ret >= (ll)M);
}

int modmul2(int a, int b, int M) {
	unsigned ret = (unsigned)a * (unsigned)b - (unsigned)M * (unsigned)(int)(1.0 / M * a * b);
	return (int)ret;
}

typedef uint64_t u64;
typedef __uint128_t u128;
struct Barrett {
	u64 b, m;
	Barrett(u64 b) : b(b), m(-1ULL / b) {}
	u64 reduce(u64 a) {
		u64 q = (u64)((u128(m) * a) >> 64), r = a - q * b;
		return r - (r >= b) * b;
	}
};

struct RelaxedBarrett {
	u64 b, m;
	RelaxedBarrett(u64 b) : b(b), m(-1ULL / b) {}
	u64 reduce(u64 a) {
		return a - (u64)((u128(m) * a) >> 64) * b;
	}
};

struct N {
	int x = 0;
};

struct Mont {
	int Mod, R1Mod, R2Mod, NPrime;

	Mont(int mod);

	N redc(int a, int b);
	N raw(int x) { N r; r.x = x; return r; }
	N from(int x) { assert (x < Mod); return redc(x, R2Mod); }
	N one() { return raw(R1Mod); }
	int get(N a) { return redc(a.x, 1).x; }

	N mul(N a, N b) { return redc(a.x, b.x); }
	N add(N a, N b) { int x = a.x + b.x; if (x >= Mod) x -= Mod; return raw(x); }
	N sub(N a, N b) { int x = a.x - b.x; if (a.x < b.x) x += Mod; return raw(x); }
};

Mont::Mont(int mod) : Mod(mod) {
	const ll B = (1LL << 32);
	assert((mod & 1) != 0);
	ll R = B % mod;
	ll xinv = 1, bit = 2;
	for (int i = 1; i < 32; i++, bit <<= 1) { // Hensel lifting!
		ll y = xinv * mod;
		if ((y & bit) != 0)
			xinv |= bit;
	}
	assert(((mod * xinv) & (B-1)) == 1);
	R1Mod = (int)R;
	R2Mod = (int)(R * R % mod);
	NPrime = (int)(B - xinv);
}

N Mont::redc(int a, int b) {
	ll T = (ll)a * b;
	ll m = (unsigned)T * (unsigned)NPrime;
	T += m * Mod;
	T >>= 32;
	if (T >= Mod)
		T -= Mod;
	return raw((int)T);
}

struct RelaxedMont {
	int Mod, R1Mod, R2Mod, NPrime;

	RelaxedMont(int mod);

	N redc(int a, int b);
	N raw(int x) { N r; r.x = x; return r; }
	N from(int x) { assert (x < Mod); return redc(x, R2Mod); }
	N one() { return raw(R1Mod); }
	int get(N a) { return redc(a.x, 1).x; }

	N mul(N a, N b) { return redc(a.x, b.x); }
};

RelaxedMont::RelaxedMont(int mod) : Mod(mod) {
	const ll B = (1LL << 32);
	assert((mod & 1) != 0);
	ll R = B % mod;
	ll xinv = 1, bit = 2;
	for (int i = 1; i < 32; i++, bit <<= 1) { // Hensel lifting!
		ll y = xinv * mod;
		if ((y & bit) != 0)
			xinv |= bit;
	}
	assert(((mod * xinv) & (B-1)) == 1);
	R1Mod = (int)R;
	R2Mod = (int)(R * R % mod);
	NPrime = (int)(B - xinv);
}

N RelaxedMont::redc(int a, int b) {
	ll T = (ll)a * b;
	ll m = (unsigned)T * (unsigned)NPrime;
	T += m * Mod;
	T >>= 32;
	return raw((int)T);
}

const int M = 1'000'000'007;
int M_dynamic = M;

int main(int argc, char** argv) {
	cin.sync_with_stdio(false);
	cin.exceptions(cin.failbit);
	int method = atoi(argv[1]);
	if (method == 0) {
		// Really naive, with dynamic modulo. 19.161s.
		ll prod = 1;
		for (int i = 1; i < M; i++) {
			prod = prod * i % M_dynamic;
		}
		cout << prod << endl;
	}
	if (method == 1) {
		// Naive. 5.569s.
		ll prod = 1;
		for (int i = 1; i < M; i++) {
			prod = prod * i % M;
		}
		cout << prod << endl;
	}
	else if (method == 2) {
		// Parallel really naive. 10.423s.
		const int PAR = 8;
		ll prods[PAR];
		rep(i,0,PAR) prods[i] = 1;
		int i = 1;
		for (; i + PAR <= M;) {
			rep(j,0,PAR)
				prods[j] = prods[j] * i % M_dynamic, i++;
		}
		ll prod = 1;
		rep(i,0,PAR) prod = prod * prods[i] % M_dynamic;
		while (i < M) {
			prod = prod * i % M_dynamic; i++;
		}
		cout << prod << endl;
	}
	else if (method == 3) {
		// Parallel, to avoid latency bottlenecks. 1.453s.
		const int PAR = 8;
		ll prods[PAR];
		rep(i,0,PAR) prods[i] = 1;
		int i = 1;
		for (; i + PAR <= M;) {
			rep(j,0,PAR)
				prods[j] = prods[j] * i % M, i++;
		}
		ll prod = 1;
		rep(i,0,PAR) prod = prod * prods[i] % M;
		while (i < M) {
			prod = prod * i % M; i++;
		}
		cout << prod << endl;
	}
	else if (method == 4) {
		// Floating-point modmul (like KACTL's version but with doubles). 3.088s.
		const int PAR = 8;
		unsigned prods[PAR];
		rep(i,0,PAR) prods[i] = 1;
		int i = 1;
		for (; i + PAR <= M;) {
			rep(j,0,PAR)
				prods[j] = modmul(prods[j], i, M), i++;
		}
		ll prod = 1;
		rep(i,0,PAR) prod = prod * prods[i] % M;
		while (i < M) {
			prod = prod * i % M; i++;
		}
		cout << prod << endl;
	}
	else if (method == 5) {
		// Relaxed floating-point modmul (without the final reduction in the
		// function, allowing for negative and out-of-range numbers). 2.024s.
		const int PAR = 8;
		int prods[PAR];
		rep(i,0,PAR) prods[i] = 1;
		int i = 1;
		for (; i + PAR <= M;) {
			rep(j,0,PAR)
				prods[j] = modmul2(prods[j], i, M), i++;
		}
		ll prod = 1;
		rep(i,0,PAR) prod = prod * (prods[i] % M) % M;
		if (prod < 0) prod += M;
		while (i < M) {
			prod = prod * i % M; i++;
		}
		cout << prod << endl;
	}
	else if (method == 6) {
		// SIMD using relaxed floating-point modmul. 0.691s.
		const int PAR = 8;
		typedef __m128i mi;
		typedef __m256d md;
		mi ones = _mm_set1_epi32(1);
		mi prods[PAR];
		mi iaccs[PAR];
		rep(i,0,PAR) {
			prods[i] = ones;
			iaccs[i] = _mm_setr_epi32(4*i+1, 4*i+2, 4*i+3, 4*i+4);
		}
		mi iaccadd = _mm_set1_epi32(4 * PAR);
		mi ms = _mm_set1_epi32(M);
		md minv = _mm256_set1_pd(1.0 / M);
		int i = 1;
		for (; i + 4 * PAR <= M; i += 4 * PAR) {
			rep(j,0,PAR) {
				mi a = prods[j];
				mi b = iaccs[j];
				iaccs[j] = _mm_add_epi32(b, iaccadd);
				mi ab = _mm_mullo_epi32(a, b);
				mi fltresult = _mm256_cvtpd_epi32(
					_mm256_mul_pd(
						_mm256_mul_pd(minv,
							_mm256_cvtepi32_pd(b)),
						_mm256_cvtepi32_pd(a)
					)
				);
				mi res = _mm_sub_epi32(ab, _mm_mullo_epi32(ms, fltresult));
				prods[j] = res;
			}
		}
		ll prod = 1;
		rep(i,0,PAR) {
			union {
				int i[4];
				mi m;
			} u;
			u.m = prods[i];
			rep(j,0,4) prod = prod * (u.i[j] % M) % M;
		}
		if (prod < 0) prod += M;
		while (i < M) {
			prod = prod * i % M; i++;
		}
		cout << prod << endl;
	}
	else if (method == 7) {
		// Barrett reduction. 1.637s.
		const int PAR = 8;
		Barrett ba(M);
		u64 prods[PAR];
		rep(i,0,PAR) prods[i] = 1;
		int i = 1;
		for (; i + PAR <= M;) {
			rep(j,0,PAR)
				prods[j] = ba.reduce(prods[j] * i), i++;
		}
		u64 prod = 1;
		rep(i,0,PAR) prod = ba.reduce(prod * prods[i]);
		while (i < M) {
			prod = ba.reduce(prod * i), i++;
		}
		cout << prod << endl;
	}
	else if (method == 8) {
		// Relaxed Barrett reduction. 1.195s.
		const int PAR = 8;
		RelaxedBarrett ba(M);
		u64 prods[PAR];
		rep(i,0,PAR) prods[i] = 1;
		int i = 1;
		for (; i + PAR <= M;) {
			rep(j,0,PAR)
				prods[j] = ba.reduce(prods[j] * i), i++;
		}
		u64 prod = 1;
		rep(i,0,PAR) prod = ba.reduce(prod * prods[i]);
		while (i < M) {
			prod = ba.reduce(prod * i), i++;
		}
		prod %= M;
		cout << prod << endl;
	}
	else if (method == 9) {
		// Montgomery multiplication. 1.668s.
		const int PAR = 8;
		Mont mont(M);
		N prods[PAR];
		rep(i,0,PAR) prods[i] = mont.one();
		int i = 1;
		for (; i + PAR <= M;) {
			rep(j,0,PAR)
				prods[j] = mont.mul(prods[j], mont.raw(i)), i++;
		}
		N prod = mont.one();
		rep(i,0,PAR) prod = mont.mul(prod, prods[i]);
		while (i < M) {
			prod = mont.mul(prod, mont.raw(i)), i++;
		}
		// We ought to multiply by R^(M-1) to account for the non-Montgomery
		// form numbers that got multiplied in. But that's 1, so no need.
		cout << mont.get(prod) << endl;
	}
	else if (method == 10) {
		// Relaxed Montgomery multiplication. 1.232s.
		const int PAR = 8;
		RelaxedMont mont(M);
		N prods[PAR];
		rep(i,0,PAR) prods[i] = mont.one();
		int i = 1;
		for (; i + PAR <= M;) {
			rep(j,0,PAR)
				prods[j] = mont.mul(prods[j], mont.raw(i)), i++;
		}
		N prod = mont.one();
		rep(i,0,PAR) prod = mont.mul(prod, prods[i]);
		while (i < M) {
			prod = mont.mul(prod, mont.raw(i)), i++;
		}
		// We ought to multiply by R^(M-1) to account for the non-Montgomery
		// form numbers that got multiplied in. But that's 1, so no need.
		cout << mont.get(prod) << endl;
	}
	else if (method == 11) {
		// SIMD Montgomery multiplication. 0.493s.
		const int PAR = 8;
		typedef __m256i mi;
		Mont mont(M);
		mi prods[PAR];
		mi accs[PAR];
		rep(i,0,PAR) {
			prods[i] = _mm256_set1_epi64x(mont.one().x);
			accs[i] = _mm256_setr_epi64x(4*i + 1, 4*i + 2, 4*i + 3, 4*i + 4);
		}
		mi iaccadd = _mm256_set1_epi64x(4 * PAR);
		mi mnprime = _mm256_set1_epi64x(mont.NPrime);
		mi mmod = _mm256_set1_epi64x(mont.Mod);
		int i = 1;
		for (; i + PAR * 4 <= M; i += PAR * 4) {
			rep(j,0,PAR) {
				mi a = prods[j];
				mi b = accs[j];
				accs[j] = _mm256_add_epi64(b, iaccadd);
				mi T = _mm256_mul_epu32(a, b);
				mi m = _mm256_mul_epu32(T, mnprime); // uses lo 32 bits of T
				T = _mm256_add_epi64(T, _mm256_mul_epu32(m, mmod)); // uses lo 32 bits of m
				T = _mm256_srli_epi64(T, 32);
				T = _mm256_sub_epi64(T,
					_mm256_andnot_si256(
						_mm256_cmpgt_epi32(mmod, T),
						mmod));
				prods[j] = T;
			}
		}
		N prod = mont.one();
		rep(i,0,PAR) {
			union {
				u64 i[4];
				mi m;
			} u;
			u.m = prods[i];
			rep(j,0,4)
				prod = mont.mul(prod, mont.raw((int)u.i[j]));
		}
		while (i < M) {
			prod = mont.mul(prod, mont.raw(i)), i++;
		}
		// We ought to multiply by R^(M-1) to account for the non-Montgomery
		// form numbers that got multiplied in. But that's 1, so no need.
		cout << mont.get(prod) << endl;
	}
	else if (method == 12) {
		// SIMD Relaxed Montgomery multiplication. 0.435s.
		const int PAR = 8;
		typedef __m256i mi;
		Mont mont(M);
		mi prods[PAR];
		mi accs[PAR];
		rep(i,0,PAR) {
			prods[i] = _mm256_set1_epi64x(mont.one().x);
			accs[i] = _mm256_setr_epi64x(4*i + 1, 4*i + 2, 4*i + 3, 4*i + 4);
		}
		mi iaccadd = _mm256_set1_epi64x(4 * PAR);
		mi mnprime = _mm256_set1_epi64x(mont.NPrime);
		mi mmod = _mm256_set1_epi64x(mont.Mod);
		int i = 1;
		for (; i + PAR * 4 <= M; i += PAR * 4) {
			rep(j,0,PAR) {
				mi a = prods[j];
				mi b = accs[j];
				accs[j] = _mm256_add_epi64(b, iaccadd);
				mi T = _mm256_mul_epu32(a, b);
				mi m = _mm256_mul_epu32(T, mnprime); // uses lo 32 bits of T
				T = _mm256_add_epi64(T, _mm256_mul_epu32(m, mmod)); // uses lo 32 bits of m
				T = _mm256_srli_epi64(T, 32);
				prods[j] = T;
			}
		}
		N prod = mont.one();
		rep(i,0,PAR) {
			union {
				u64 i[4];
				mi m;
			} u;
			u.m = prods[i];
			rep(j,0,4)
				prod = mont.mul(prod, mont.raw((int)u.i[j]));
		}
		while (i < M) {
			prod = mont.mul(prod, mont.raw(i)), i++;
		}
		// We ought to multiply by R^(M-1) to account for the non-Montgomery
		// form numbers that got multiplied in. But that's 1, so no need.
		cout << mont.get(prod) % M << endl;
	}

	exit(0);
}
