#include <bits/stdc++.h>
using namespace std;

#define rep(i, from, to) for (int i = from; i < (to); ++i)
#define trav(a, x) for (auto& a : x)
#define all(x) x.begin(), x.end()
#define sz(x) (int)(x).size()
typedef long long ll;
typedef pair<int, int> pii;
typedef vector<int> vi;

int read() {
    int x;
    cin >> x;
    return x;
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
    cin.exceptions(cin.failbit);
    int N = read(), Q = read();
    vi bits(N);
    rep(i,0,N) cin >> bits[i];
    rep(qi,0,Q) {
        int tp = read(), a = read(), b = read();
        --a;
        if (tp == 3 || tp == 5) b--;
        if (tp == 4 || tp == 6) a++;
        assert(0 <= a && a <= b && b <= N);
        if (tp == 3 || tp == 5) assert(b < N);
        if (tp == 4 || tp == 6) assert(a > 0);

        if (tp == 1) {
            for (int i = a; i < b; ++i)
                bits[i] = 0;
        } else if (tp == 2) {
            for (int i = a; i < b; ++i)
                bits[i] = 1;
        } else if (tp == 3) {
            for (int i = a; i < b; ++i)
                bits[i] |= bits[i+1];
        } else if (tp == 4) {
            for (int i = b-1; i >= a; --i)
                bits[i] |= bits[i-1];
        } else if (tp == 5) {
            for (int i = a; i < b; ++i)
                bits[i] &= bits[i+1];
        } else if (tp == 6) {
            for (int i = b-1; i >= a; --i)
                bits[i] &= bits[i-1];
        } else if (tp == 7) {
            int sum = 0;
            for (int i = a; i < b; ++i)
                sum += bits[i];
            cout << sum << '\n';
        }
    }

    // rep(i,0,N) cout << (bits[i] ? '1' : '0');
    // cout << endl;

    exit(0);
}
