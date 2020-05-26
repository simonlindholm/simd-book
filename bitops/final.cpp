#ifdef USE_AVX256
#pragma GCC target ("avx2")
#else
#define NDEBUG
#pragma GCC target ("avx512f,avx512bw")
#endif
#pragma GCC optimize ("unroll-loops")
#pragma GCC optimize ("O4")
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
    char buf[1 << 17];
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
    int a = c - '0';
    while ((c = gc()) >= 40) a = a * 10 + c - '0';
    return a;
}

static char printbuf[1'000'000 * 8];
static int printbufi = 0;
static void printchar(char c) {
    printbuf[printbufi++] = c;
}
static void print(int x) {
    char buf[20];
    int qr = 0;
    if (!x) printchar('0');
    else {
        if (x < 0) printchar('-'), x = -x;
        while (x) buf[qr++] = (char)(x % 10 + '0'), x /= 10;
        while (qr) printchar(buf[--qr]);
    }
    printchar('\n');
}

#ifdef USE_AVX256
typedef __m256i M;
#define MM(op) _mm256_ ## op
#define MM_ANDNOT(a,b) _mm256_andnot_si256(a,b)
#define ZERO MM(setzero_si256)()
#define MM_SHUFFLE_MASK(a,b,c,d) \
    _mm256_set_epi32(a,b,c,d, a,b,c,d)
#define M_WORDS 4
#else
typedef __m512i M;
#define MM(op) _mm512_ ## op
#define MM_ANDNOT(a,b) _mm512_andnot_si512(a,b)
#define ZERO MM(setzero_si512)()
#define MM_SHUFFLE_MASK(a,b,c,d) \
    _mm512_set_epi32(a,b,c,d, a,b,c,d, a,b,c,d, a,b,c,d)
#define M_WORDS 8
#endif

#define ONE MM(set1_epi32)(-1)
#define M_BITS (M_WORDS * 64)

union U {
    M m;
    uint64_t words[M_WORDS];
    uint16_t shorts[4*M_WORDS];
};

// m >> 1
U srl1(U u) {
    rep(i,0,M_WORDS-1)
        u.words[i] = (u.words[i] >> 1) | (u.words[i + 1] << 63);
    u.words[M_WORDS-1] >>= 1;
    return u;
}

// m << 1
U sll1(U u) {
    for (int i = M_WORDS-1; i >= 1; i--)
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
    rep(i,0,M_BITS) os << (getbit(u, i) ? '1' : '0');
    return os;
}

template <class F>
void generic_for_each(F&& f) {}

template <class F, class Arg, class... Args>
void generic_for_each(F&& f, Arg& arg, Args&... args) {
    f(arg);
    generic_for_each(f, args...);
}

struct State {
    int a, b;
    int evs[8];
    int evi;
    int pos;
    int tp;
    M prev, acc;
    M affected;
    M ones, twos, fours, eights;
};

struct Temp {
    // Persistent state from State
    M prev, acc, affected;
    M ones, twos, fours, eights;
    // Temporary state
    M cur, next, local;
    Temp() {}
    Temp(State& s) :
        prev(s.prev), acc(s.acc),
        affected(s.affected),
        ones(s.ones), twos(s.twos), fours(s.fours), eights(s.eights)
    {}
    void writeBack(State& s) {
        s.prev = prev;
        s.acc = acc;
        s.affected = affected;
        s.ones = ones;
        s.twos = twos;
        s.fours = fours;
        s.eights = eights;
    }
};

struct Op {
    enum { IS_POPCNT = 0, NEED_NEXT = 0 };
};

struct Op1 : Op {
    void handle(Temp& t) {
        t.cur = MM_ANDNOT(t.affected, t.cur);
    }
};

struct Op2 : Op {
    void handle(Temp& t) {
        t.cur |= t.affected;
    }
};

struct Op3 : Op {
    enum { NEED_NEXT = 1 };
    void handle(Temp& t) {
#ifdef USE_AVX256
        t.cur |= t.next & t.affected;
#else
        t.cur = _mm512_ternarylogic_epi32(t.cur, t.next, t.affected, 0xF8);
#endif
    }
};

struct Op4 : Op {
    void handle(Temp& t) {
        M old = t.cur;
#ifdef USE_AVX256
        t.cur |= t.prev & t.affected;
#else
        t.cur = _mm512_ternarylogic_epi32(t.cur, t.prev, t.affected, 0xF8);
#endif
        t.prev = old;
    }
};

struct Op5 : Op {
    enum { NEED_NEXT = 1 };
    void handle(Temp& t) {
#ifdef USE_AVX256
        t.cur = MM_ANDNOT(MM_ANDNOT(t.next, t.affected), t.cur);
#else
        t.cur = _mm512_ternarylogic_epi32(t.cur, t.next, t.affected, 0xD0);
#endif
    }
};

struct Op6 : Op {
    void handle(Temp& t) {
        M old = t.cur;
#ifdef USE_AVX256
        t.cur = MM_ANDNOT(MM_ANDNOT(t.prev, t.affected), t.cur);
#else
        t.cur = _mm512_ternarylogic_epi32(t.cur, t.prev, t.affected, 0xD0);
#endif
        t.prev = old;
    }
};

// Harley-Seal popcount from http://0x80.pl/articles/sse-popcount.html
void CSA(M& h, M& a, M b, M c) {
#ifdef USE_AVX256
    const M u = a ^ b;
    h = (a & b) | (u & c);
    a = u ^ c;
#else
    M a2 = a;
    a = _mm512_ternarylogic_epi32(c, b, a2, 0x96);
    h = _mm512_ternarylogic_epi32(c, b, a2, 0xE8);
#endif
}

M popcount(const M v) {
    const M m1 = MM(set1_epi8)(0x55);
    const M m2 = MM(set1_epi8)(0x33);
    const M m4 = MM(set1_epi8)(0x0F);
    const M t1 = MM(sub_epi8)(v,       (MM(srli_epi16)(v,  1) & m1));
    const M t2 = MM(add_epi8)(t1 & m2, (MM(srli_epi16)(t1, 2) & m2));
    const M t3 = MM(add_epi8)(t2, MM(srli_epi16)(t2, 4)) & m4;
    return MM(sad_epu8)(t3, ZERO);
}

struct Op7 : Op {
    enum { IS_POPCNT = 1 };

    void handle16(Temp& t, U* data) {
        M twosA, twosB, foursA, foursB, eightsA, sixteensA;
        M ones = data[0].m, twos, fours, eights;
        CSA(twos, ones, data[1].m, data[2].m);
        CSA(twosA, ones, data[3].m, data[4].m);
        CSA(twosB, ones, data[5].m, data[6].m);
        CSA(fours, twos, twosA, twosB);
        CSA(twosA, ones, data[7].m, data[8].m);
        CSA(twosB, ones, data[9].m, data[10].m);
        CSA(foursA, twos, twosA, twosB);
        CSA(twosA, ones, data[11].m, data[12].m);
        CSA(twosB, ones, data[13].m, data[14].m);
        CSA(foursB, twos, twosA, twosB);
        CSA(eights, fours, foursA, foursB);

        CSA(twosA, t.ones, ones & t.affected, data[15].m & t.affected);
        CSA(foursA, t.twos, twos & t.affected, twosA);
        CSA(eightsA, t.fours, fours & t.affected, foursA);
        CSA(sixteensA, t.eights, eights & t.affected, eightsA);
        t.acc = MM(add_epi64)(t.acc, MM(slli_epi64)(popcount(sixteensA), 4));
    }

    void handle16_2(Temp& t1, Temp& t2, U* data) {
        M twosA, twosB, foursA, foursB, eightsA, sixteensA;
        M ones = data[0].m, twos, fours, eights;
        CSA(twos, ones, data[1].m, data[2].m);
        CSA(twosA, ones, data[3].m, data[4].m);
        CSA(twosB, ones, data[5].m, data[6].m);
        CSA(fours, twos, twosA, twosB);
        CSA(twosA, ones, data[7].m, data[8].m);
        CSA(twosB, ones, data[9].m, data[10].m);
        CSA(foursA, twos, twosA, twosB);
        CSA(twosA, ones, data[11].m, data[12].m);
        CSA(twosB, ones, data[13].m, data[14].m);
        CSA(foursB, twos, twosA, twosB);
        CSA(eights, fours, foursA, foursB);

        CSA(twosA, t1.ones, ones & t1.affected, data[15].m & t1.affected);
        CSA(foursA, t1.twos, twos & t1.affected, twosA);
        CSA(eightsA, t1.fours, fours & t1.affected, foursA);
        CSA(sixteensA, t1.eights, eights & t1.affected, eightsA);
        t1.acc = MM(add_epi64)(t1.acc, MM(slli_epi64)(popcount(sixteensA), 4));

        CSA(twosA, t2.ones, ones & t2.affected, data[15].m & t2.affected);
        CSA(foursA, t2.twos, twos & t2.affected, twosA);
        CSA(eightsA, t2.fours, fours & t2.affected, foursA);
        CSA(sixteensA, t2.eights, eights & t2.affected, eightsA);
        t2.acc = MM(add_epi64)(t2.acc, MM(slli_epi64)(popcount(sixteensA), 4));
    }

    // pshufb-based popcnt, see http://0x80.pl/articles/sse-popcount.html
    void handle(Temp& t) {
        const M lowMask = MM(set1_epi8)(0xF);
        const M popcntLookup = MM_SHUFFLE_MASK(
            0x04030302,
            0x03020201,
            0x03020201,
            0x02010100
        );
        const M vec = t.cur & t.affected;
        const M lo  = vec & lowMask;
        const M hi  = MM(srli_epi16)(vec, 4) & lowMask;
        const M popcnt1 = MM(shuffle_epi8)(popcntLookup, lo);
        const M popcnt2 = MM(shuffle_epi8)(popcntLookup, hi);
        t.local = MM(add_epi8)(t.local, popcnt1);
        t.local = MM(add_epi8)(t.local, popcnt2);
    }

    void handleNth(Temp& t) {
        t.acc = MM(add_epi64)(t.acc, MM(sad_epu8)(t.local, ZERO));
    }
};

template<class F>
static void withOp(int tp, F&& f) {
    switch (tp) {
        case 1: f(Op1()); break;
        case 2: f(Op2()); break;
        case 3: f(Op3()); break;
        case 4: f(Op4()); break;
        case 5: f(Op5()); break;
        case 6: f(Op6()); break;
        case 7: f(Op7()); break;
        default: assert(0);
    }
}

constexpr int pad = 2;
constexpr int PREFIX_PAD = 8;
constexpr int N_MAX = 1'000'000;
constexpr int nblocks = max((N_MAX + M_BITS-1) / M_BITS, pad * 2 + 5);
static State states[pad];
static U prefixes[M_BITS+1 + PREFIX_PAD*2];
static U* blocks;
// static int nblocks;
static int nops;

void computeAffected(State& s) {
    int at = s.pos;
    if (at < s.evs[s.evi])
        return;

    while (at >= s.evs[s.evi])
        s.evi++;

    // Compute mask for which of the entries corresponding to bits
    // 0..63 should be updated when iterating over the range
    // [s.evs[s.evi-1], s.evs[s.evi]). Because of how "evs" was
    // constructed, this will be constant over that range (if we
    // interpret infinities correctly) and we can sample a random
    // point "at" in the interval for computing it.
    //
    // The mask should have bit i set if:
    // a <= i * nblocks + at < b
    // (a - at) / nblocks <= i < (b - at) / nblocks
    // ceil((a - at) / nblocks) <= i < ceil((b - at) / nblocks)
    //
    // We assume this range is contained in [-PREFIX_PAD, M_BITS + PREFIX_PAD].
    auto ceildiv = [](int a) {
        // a can be negative, causing a badly rounded division, but not by very
        // much. Compensate by adding PREFIX_PAD*b.
        return (a + PREFIX_PAD*nblocks + nblocks-1) / nblocks - PREFIX_PAD;
    };
    int low = ceildiv(s.a - at);
    int high = ceildiv(s.b - at);
    assert(low <= high);
    assert(-PREFIX_PAD <= low);
    assert(high <= M_BITS + PREFIX_PAD);
    s.affected = MM_ANDNOT(prefixes[PREFIX_PAD + low].m, prefixes[PREFIX_PAD + high].m);
}

void advanceOneSmall(State& s, int goal) {
    int at = s.pos;
    assert(at < goal);
    computeAffected(s);

    int from = at, to = min(goal, s.evs[s.evi]);
    assert(from < to);
    Temp t(s);
    withOp(s.tp, [&](auto op) {
        rep(i,from,to) {
            t.local = ZERO;
            t.cur = blocks[i].m;
            if constexpr (op.NEED_NEXT) t.next = blocks[i+1].m;
            op.handle(t);
            if constexpr (op.IS_POPCNT) op.handleNth(t);
            else blocks[i].m = t.cur;
        }
    });
    t.writeBack(s);
    s.pos = to;
}

void advanceSmall(State& s, int goal) {
    while (s.pos != goal)
        advanceOneSmall(s, goal);
}

void advanceOneLargePopcounts(State& s1, State& s2, int goal) {
    int at = s1.pos;
    assert(at < goal);
    assert(at == s2.pos);
    computeAffected(s1);
    computeAffected(s2);

    int from = at, to = min({goal, s1.evs[s1.evi], s2.evs[s2.evi]});
    assert(from < to);
    Temp t1(s1);
    Temp t2(s2);
    Op7 op;

    int i = from;
    while (i + 16 <= to) {
        op.handle16_2(t1, t2, &blocks[i]);
        i += 16;
    }

    // Remainder (TODO: use Harley-Seal here as well)
    t1.local = ZERO;
    t2.local = ZERO;
    while (i < to) {
        t1.cur = t2.cur = blocks[i].m;
        op.handle(t1);
        op.handle(t2);
        i++;
    }
    op.handleNth(t1);
    op.handleNth(t2);

    t1.writeBack(s1);
    t2.writeBack(s2);
    s1.pos = to;
    s2.pos = to;
}

void advanceLargePopcounts(State& s1, State& s2, int goal) {
    while (s1.pos != goal)
        advanceOneLargePopcounts(s1, s2, goal);
}

constexpr bool onlyPopcount() { return true; }
template<class Op, class... Ops>
constexpr bool onlyPopcount(Op op, Ops... ops) {
    return op.IS_POPCNT && onlyPopcount(ops...);
}

template<class... Ops>
void runLockstepLoop(int start, int goal, Ops... ops) {
    constexpr int nops = (int) sizeof...(ops);
    assert(start < goal);
    // Ideally we'd use an array here:
    // Temp temps[nops];
    // but that confuses GCC's alias analysis. Work around this by declaring
    // two separate variables.
    Temp temp1(states[0]), temp2;
    if constexpr (nops == 2) temp2 = Temp(states[1]);
#define TEMP(ind) (ind ? temp2 : temp1)
    rep(i,0,nops) assert(states[i].pos == start + (nops-1 - i));

    constexpr int CHUNK = 16;
    int i = start;
    for (; i + CHUNK <= goal; i += CHUNK) {
        int ind = 0;
        generic_for_each([&](auto op) {
            Temp& t = TEMP(ind);
            if constexpr (op.IS_POPCNT) {
                op.handle16(t, &blocks[i + (nops-1 - ind)]);
            } else {
                rep(it,0,CHUNK) {
                    int j = i + it + (nops-1 - ind);
                    t.cur = blocks[j].m;
                    if constexpr (op.NEED_NEXT) t.next = blocks[j+1].m;
                    op.handle(t);
                    blocks[j].m = t.cur;
                }
            }
            ind++;
        }, ops...);
    }

    temp1.local = ZERO;
    temp2.local = ZERO;
    for (; i < goal; ++i) {
        int ind = 0;
        generic_for_each([&](auto op) {
            int j = i + (nops-1 - ind);
            Temp& t = TEMP(ind);
            t.cur = blocks[j].m;
            if constexpr (op.NEED_NEXT) t.next = blocks[j+1].m;
            op.handle(t);
            if constexpr (!op.IS_POPCNT) blocks[j].m = t.cur;
            ind++;
        }, ops...);
    }
    {
        int ind = 0;
        generic_for_each([&](auto op) {
            if constexpr (op.IS_POPCNT) op.handleNth(TEMP(ind));
            ind++;
        }, ops...);
    }

    temp1.writeBack(states[0]);
    if constexpr(nops == 2) temp2.writeBack(states[1]);
    rep(i,0,nops) states[i].pos = goal + (nops-1 - i);
}

template<class... Ops>
void loopLockstep(int goal, Ops&... ops) {
    constexpr int nops = sizeof...(ops);
    if constexpr (nops == 0) return;
    int pos = states[nops-1].pos;
    while (pos < goal) {
        int to = goal;
        rep(i,0,nops) {
            computeAffected(states[i]);
            to = min(to, states[i].evs[states[i].evi] - (nops-1 - i));
        }
        runLockstepLoop(pos, to, ops...);
        pos = to;
    }
}

template<class... Ops>
static void recLockstep(int goal, Ops&&... ops) {
    constexpr int i = sizeof...(ops);
    if constexpr (i > pad) {
        assert(0);
    } else if (i == nops) {
        loopLockstep(goal, ops...);
    } else if constexpr (i < pad) {
        withOp(states[i].tp, [&](auto op) {
            recLockstep(goal, ops..., op);
        });
    } else {
        assert(0);
    }
}

int main() {
    // Given a bit array A of size N, performs Q operations of the types:
    // 1. for every i in [a, b), set A[i] = 0
    // 2. for every i in [a, b), set A[i] = 1
    // 3. for every i in [a, b) set A[i] = A[i] | A[i+1]
    // 4. for every i in [a, b) set A[i] = A[i] | A[i-1]
    // 5. for every i in [a, b) set A[i] = A[i] & A[i+1]
    // 6. for every i in [a, b) set A[i] = A[i] & A[i-1]
    // 7. compute popcount of [a, b)
    cin.tie(0)->sync_with_stdio(false);
    int N = read(), Q = read();
    string str(N, '0');
    while (!isdigit(gc.peek())) gc();
    rep(i,0,N) {
        str[i] = gc();
#ifndef DENSE_INPUT
        gc();
#endif
    }
    // nblocks = max((N + M_BITS-1) / M_BITS, pad * 2 + 5);

    prefixes[0+PREFIX_PAD].m = ZERO;
    rep(i,0,M_BITS) {
        prefixes[i+1+PREFIX_PAD] = prefixes[i+PREFIX_PAD];
        setbit(prefixes[i+1+PREFIX_PAD], i);
    }
    rep(i,0,PREFIX_PAD) prefixes[i].m = ZERO;
    rep(i,0,PREFIX_PAD) prefixes[PREFIX_PAD+M_BITS+1 + i].m = ONE;

    // Order the bits so that block i -- a 256-bit AVX vector -- consists of
    // bits nblocks*j + i for j = 0..255. With that representation, shifting
    // left/right by one bit becomes shifting left/right by one block, except
    // at the edges. Picture: (with W = nblocks)
    //
    // 0         | 1         | 2         | 3      | 4     | ... | W-1
    // W         | W+1       | W+2       | W+3    | W+4   |     | 2*W-1
    // ...
    // 255*W     | ...       |           | N-1    | 0     | ... | 0
    // blocks[0] | blocks[1] | blocks[2] |        |       | ... | blocks[W-1]

    blocks = (U*)_mm_malloc((nblocks + 2 * pad) * sizeof(U), sizeof(U)) + pad;
    rep(i,0,N) {
        int bit = i / nblocks, block = i % nblocks;
        assert(bit < M_BITS);
        if (str[i] == '1') setbit(blocks[block], bit);
    }

    int qi = 0;
    while (qi < Q) {
        nops = min(Q - qi, pad);
        qi += nops;

        // Explicitly compute the edge wrap-around and store it as padding
        // around the actual blocks.
        rep(i,0,pad) {
            blocks[i - pad] = sll1(blocks[nblocks + i - pad]);
            blocks[nblocks + i] = srl1(blocks[i]);
        }

        rep(i,0,nops) {
            int tp = read(), a = read(), b = read();
            --a;
            if (tp == 3 || tp == 5) b--;
            if (tp == 4 || tp == 6) a++;
            assert(0 <= a && a <= b && b <= N);
            if (tp == 3 || tp == 5) assert(b < N);
            if (tp == 4 || tp == 6) assert(a > 0);

            // tp = 7;

            State& s = states[i];
            s.tp = tp;
            s.a = a;
            s.b = b;
            // s.affected = ZERO;
            if (tp == 7) {
                s.acc = ZERO;
                s.ones = ZERO;
                s.twos = ZERO;
                s.fours = ZERO;
                s.eights = ZERO;
                s.pos = 0;
            } else {
                // One too much, just to compute "prev". It only results in unused
                // padding being incorrect.
                s.prev = ZERO;
                s.pos = -nops + i;
            }
            s.evi = 0;
            int* evs = s.evs;
            evs[0] = INT_MIN;
            evs[1] = 0;
            evs[2] = a % nblocks;
            evs[3] = b % nblocks;
            evs[4] = nblocks;
            int nevs = 5;
            if (tp != 7) for (int x : {evs[2], evs[3]}) for (int delta : {-nblocks, nblocks}) {
                int y = x + delta;
                if (-pad <= y && y <= nblocks + pad) evs[nevs++] = y;
            }
            assert(nevs <= 7);
            evs[nevs++] = INT_MAX;
            sort(evs + 1, evs + nevs - 1);
        }

        if (nops == 2 && states[0].tp == 7 && states[1].tp == 7) {
            advanceLargePopcounts(states[0], states[1], nblocks);
        } else {
            // Initial part
            rep(i,0,nops) {
                advanceSmall(states[i], nops-1 - i);
            }

            // Middle part, in lock-step, until where states with tp == 7 need to stop
            int goal = nblocks - nops + 1;
            if (nops == 2 && states[0].tp != 7)
                goal++;
            recLockstep(goal);

            // Final part
            rep(i,0,nops) {
                int to = nblocks + (states[i].tp == 7 ? 0 : nops-1 - i);
                advanceSmall(states[i], to);
            }
        }

        rep(i,0,nops) if (states[i].tp == 7) {
            State& s = states[i];
            M acc = s.acc;
            acc = MM(add_epi64)(acc, MM(slli_epi64)(popcount(s.eights), 3));
            acc = MM(add_epi64)(acc, MM(slli_epi64)(popcount(s.fours), 2));
            acc = MM(add_epi64)(acc, MM(slli_epi64)(popcount(s.twos), 1));
            acc = MM(add_epi64)(acc, popcount(s.ones));
            ll sum = 0;
            U u;
            u.m = acc;
            rep(j,0,M_WORDS) sum += u.words[j];
            print((int) sum);
        }
    }

    fwrite(printbuf, printbufi, 1, stdout);
    fflush(stdout);
    _Exit(0);
}
