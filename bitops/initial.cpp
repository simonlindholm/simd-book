#pragma GCC target ("avx2")
#pragma GCC optimize ("O3")
#pragma GCC optimize ("unroll-loops")
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

struct GC {
    char buf[1 << 16];
    size_t bc = 0, be = 0;
    char peek() {
        if (bc >= be) {
            bc = 0, buf[0] = 0;
            be = fread(buf, 1, sizeof(buf), stdin);
        }
        return buf[bc]; // returns 0 on EOF
    }
    char operator()() {
        if (bc >= be) {
            bc = 0, buf[0] = 0;
            be = fread(buf, 1, sizeof(buf), stdin);
        }
        return buf[bc++]; // returns 0 on EOF
    }
} gc;
int read() {
    char c;
    while ((c = gc()) < 40);
    if (c == '-') return -read();
    int a = c - '0';
    while (isdigit(c = gc())) a = a * 10 + c -'0';
    return a;
}

typedef __m256i M;
union U {
    M m;
    uint64_t words[4];
    uint16_t shorts[16];
};

// m >> 1
U srl1(U u) {
    rep(i,0,3)
        u.words[i] = (u.words[i] >> 1) | (u.words[i + 1] << 63);
    u.words[3] >>= 1;
    return u;
}

// m << 1
U sll1(U u) {
    for (int i = 3; i >= 1; i--)
        u.words[i] = (u.words[i] << 1) | (u.words[i - 1] >> 63);
    u.words[0] <<= 1;
    return u;
}

void setbit(U& u, unsigned bit) {
    u.words[bit >> 6] |= 1ULL << (bit & 63);
}

bool getbit(U& u, unsigned bit) {
    return u.words[bit >> 6] & (1ULL << (bit & 63));
}

ostream& operator<<(ostream& os, U u) {
    rep(i,0,256) os << (getbit(u, i) ? '1' : '0');
    return os;
}

int main() {
    // Given a bit array A of size N and Q operations of the types:
    // 1. for every i in [a, b), set A[i] = 0
    // 2. for every i in [a, b), set A[i] = 1
    // 3. for every i in [a, b) set A[i] = A[i] | A[i+1]
    // 4. for every i in [a, b) set A[i] = A[i] | A[i-1]
    // 5. for every i in [a, b) set A[i] = A[i] & A[i+1]
    // 6. for every i in [a, b) set A[i] = A[i] & A[i-1]
    // 7. compute popcount of [a, b)
    // Computes the array A after all operations have been performed.
    cin.tie(0);
    // cin.tie(0)->sync_with_stdio(false);
    // cin.exceptions(cin.failbit);
    int N = read(), Q = read();
    string str(N, '0');
    while (!isdigit(gc.peek())) gc();
    rep(i,0,N) {
        str[i] = gc();
#ifndef DENSE_INPUT
        gc();
#endif
    }
    int pad = 1;
    int nblocks = max((N + 255) / 256, pad);

    // Order the bits so that block i -- a 256-bit AVX vector -- consists of
    // bits nblocks*j + i for j = 0..255. With that representation, shifting
    // left/right by one bit becomes shifting left/right by one block, except
    // at the edges. Picture (with W = nblocks):
    //
    // 0         | 1         | 2         | 3      4      ... | W-1
    // W         | W+1       | W+2       | W+3    W+4        | 2*W-1
    // ...
    // 256*W,    | ...       |           | N-1,   0,     ... | 0
    // blocks[0] | blocks[1] | blocks[2] |               ... | blocks[W-1]

    U *blocks = (U*)_mm_malloc((nblocks + 2 * pad) * sizeof(U), sizeof(U));
    rep(i,0,N) {
        int bit = i / nblocks, block = i % nblocks;
        assert(bit < 256);
        if (str[i] == '1') setbit(blocks[pad + block], bit);
    }

    const M one = _mm256_set1_epi32(-1);
    const M popcntLookup = _mm256_setr_epi8(
        /* 0 */ 0, /* 1 */ 1, /* 2 */ 1, /* 3 */ 2,
        /* 4 */ 1, /* 5 */ 2, /* 6 */ 2, /* 7 */ 3,
        /* 8 */ 1, /* 9 */ 2, /* a */ 2, /* b */ 3,
        /* c */ 2, /* d */ 3, /* e */ 3, /* f */ 4,

        /* 0 */ 0, /* 1 */ 1, /* 2 */ 1, /* 3 */ 2,
        /* 4 */ 1, /* 5 */ 2, /* 6 */ 2, /* 7 */ 3,
        /* 8 */ 1, /* 9 */ 2, /* a */ 2, /* b */ 3,
        /* c */ 2, /* d */ 3, /* e */ 3, /* f */ 4
    );

    int evs[4];
    rep(qi,0,Q) {
        int tp = read(), a = read(), b = read();
        --a;
        if (tp == 3 || tp == 5) b--;
        if (tp == 4 || tp == 6) a++;
        assert(0 <= a && a <= b && b <= N);
        if (tp == 3 || tp == 5) assert(b < N);
        if (tp == 4 || tp == 6) assert(a > 0);

        evs[0] = 0;
        evs[1] = a % nblocks;
        evs[2] = b % nblocks;
        evs[3] = nblocks;
        sort(evs, evs + 4);

        // Explicitly compute the edge wrap-around and store it as padding
        // around the actual blocks.
        rep(i,0,pad) {
            blocks[i] = sll1(blocks[nblocks + i]);
            blocks[nblocks + pad + i] = srl1(blocks[pad + i]);
        }

        int sum = 0;
        auto process = [&](int from, int to) {
            if (from == to) return;
            assert(from < to);
            // Compute mask for which of the entries corresponding to bits
            // 0..63 should be updated for the range [from, to). Because of
            // how "evs" was constructed, this will be constant over the range
            // [from, to). Thus, compute the constant for "from".
            U u;
            u.m = _mm256_setzero_si256();
            rep(wi,0,4) {
                int i0 = wi*64 + 0;
                int i1 = wi*64 + 63;
                int x0 = i0 * nblocks + from;
                int x1 = i1 * nblocks + from;
                if (x1 < a || b <= x0) {
                    u.words[wi] = 0;
                } else if (a <= x0 && x1 < b) {
                    u.words[wi] = -1ULL;
                } else {
                    rep(si,0,4) {
                        int i0 = wi*64 + si*16 + 0;
                        int i1 = wi*64 + si*16 + 15;
                        int x0 = i0 * nblocks + from;
                        int x1 = i1 * nblocks + from;
                        if (x1 < a || b <= x0) {
                            u.shorts[wi*4 + si] = 0;
                        } else if (a <= x0 && x1 < b) {
                            u.shorts[wi*4 + si] = (uint16_t)-1U;
                        } else {
                            rep(j,0,16) {
                                int i = wi*64 + si*16 + j;
                                int x = i * nblocks + from;
                                u.shorts[wi*4 + si] |= (uint16_t)((a <= x && x < b) << j);
                            }
                        }
                    }
                }
            }
            M affected = u.m;
            M naffected = _mm256_andnot_si256(affected, one);

            from += pad;
            to += pad;
            if (tp == 1) {
                for (int i = from; i < to; ++i)
                    blocks[i].m &= naffected;
            } else if (tp == 2) {
                for (int i = from; i < to; ++i)
                    blocks[i].m |= affected;
            } else if (tp == 3) {
                for (int i = from; i < to; ++i)
                    blocks[i].m |= blocks[i+1].m & affected;
            } else if (tp == 4) {
                for (int i = to-1; i >= from; --i)
                    blocks[i].m |= blocks[i-1].m & affected;
            } else if (tp == 5) {
                for (int i = from; i < to; ++i)
                    blocks[i].m &= blocks[i+1].m | naffected;
            } else if (tp == 6) {
                for (int i = to-1; i >= from; --i)
                    blocks[i].m &= blocks[i-1].m | naffected;
            } else if (tp == 7) {
                // pshufb-based popcnt, see http://0x80.pl/articles/sse-popcount.html
                const M lowMask = _mm256_set1_epi8(0xf);
                M zero = _mm256_setzero_si256(), acc = zero;
#define ITER { \
    const M vec = blocks[i].m & affected; \
        const M lo  = _mm256_and_si256(vec, lowMask); \
        const M hi  = _mm256_and_si256(_mm256_srli_epi16(vec, 4), lowMask); \
        const M popcnt1 = _mm256_shuffle_epi8(popcntLookup, lo); \
        const M popcnt2 = _mm256_shuffle_epi8(popcntLookup, hi); \
        local = _mm256_add_epi8(local, popcnt1); \
        local = _mm256_add_epi8(local, popcnt2); \
        i++; \
    }
                int i = from;
                while (i + 8 <= to) {
                    M local = zero;
                    ITER ITER ITER ITER
                    ITER ITER ITER ITER
                    acc = _mm256_add_epi64(acc, _mm256_sad_epu8(local, zero));
                }
                M local = zero;
                while (i < to) {
                    ITER;
                }
                acc = _mm256_add_epi64(acc, _mm256_sad_epu8(local, zero));
                sum += (int) _mm256_extract_epi64(acc, 0);
                sum += (int) _mm256_extract_epi64(acc, 1);
                sum += (int) _mm256_extract_epi64(acc, 2);
                sum += (int) _mm256_extract_epi64(acc, 3);
#undef ITER
            }
        };

        if (tp == 4 || tp == 6) {
            for (int i = 2; i >= 0; i--)
                process(evs[i], evs[i+1]);
        }
        else {
            for (int i = 0; i <= 2; i++)
                process(evs[i], evs[i+1]);
        }
        if (tp == 7) printf("%d\n", sum);
    }

#ifdef PRINT_FINAL
    rep(i,0,N) {
        int bit = i / nblocks, block = i % nblocks;
        assert(bit < 256);
        cout << (getbit(blocks[pad + block], bit) ? '1' : '0');
    }
    cout << endl;
#endif

    fflush(stdout);
    _Exit(0);
}
