#include <bits/stdc++.h>
using namespace std;

#define rep(i, from, to) for (int i = from; i < int(to); ++i)
#define trav(it, x) for (auto it = x.begin(); it != x.end(); ++it)
#define all(x) x.begin(), x.end()
typedef long long ll;
typedef pair<int, int> pii;
typedef vector<int> vi;
typedef unsigned int u32;

u32 modpow(u32 a, ll e) {
    if (e == 0) return 1;
    u32 x = modpow(a*a, e >> 1);
    if (e & 1) x *= a;
    return x;
}

u32 solve(ll N, ll K) {
    ll P = 0;
    vector<pair<int, int> > ms;
    for (int i = 0; i < 3; ++i) {
        ll n = (i == 0 ? N : i == 1 ? K : N-K);
        int mult = (i == 0 ? 1 : -1);
        int p2 = 0;
        P -= n * mult;
        while (n > 0) {
            ll m = (n+1) / 2;
            if (m > 1)
                ms.push_back(make_pair((int)(m % (1LL << 31)), mult));
            P += n * mult;
            n /= 2;
            p2++;
        }
    }
    if (P >= 32) {
        return 0;
    }

    sort(all(ms));
    int tot = 0;
    for (int i = 0; i < ms.size(); ++i) {
        tot += ms[i].second;
    }
    int cur = 0;
    u32 curv = 1;
    u32 res = 1, resdiv = 1;
    for (int i = 0; i < ms.size(); ++i) {
        int m = ms[i].second;
        int lim = ms[i].first;
        u32 prod = 1;
        while (cur < lim - 3) {
            prod *= (curv * (curv + 2)) * ((curv + 4) * (curv + 6));
            curv += 8;
            cur += 4;
        }
        while (cur < lim) {
            prod *= curv;
            curv += 2;
            cur++;
        }
        int abtot = tot; if (abtot < 0) abtot = -abtot;
        u32 rprod = 1;
        for (int i = 0; i < abtot; ++i)
            rprod *= prod;
        if (tot >= 0)
            res *= rprod;
        else
            resdiv *= rprod;
        tot -= m;
    }

    res <<= P;
    res *= modpow(resdiv, (1LL << 31) - 1);
    return res;
}

int main() {
    ll N, K;
    cin >> N >> K;
    cout << solve(N, K) << endl;
}
