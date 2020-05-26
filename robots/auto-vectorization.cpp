// initial version is useless, runs in >10s
// && -> & triggers auto-vectorization, 0m1.748s
// looking at asm with -S -fno-asynchronous-unwind-tables
// -> using lots of pcmpgtw, pcmpeqw, punpckhwd, pcmpeqd
// changing <= to > and negating everything gets rid of pcmpeqw and makes it run in 0m1.587s
// changing x + y to (short)(x + y) makes it start to use pcmpgtw and runs in 0m0.955s
// still punpckhwd, apparently to add to the 'res' accumulator. Try making 'res' a short -> 0m4.938s
// unsigned short -> 0m0.840s
#pragma GCC target ("mmx,sse,sse2")
#pragma GCC optimize ("O3")
#include <bits/stdc++.h>
using namespace std;

#define rep(i, from, to) for (int i = from; i < int(to); ++i)
#define trav(a, x) for (auto& a : x)
#define all(x) x.begin(), x.end()
#define sz(x) (int)(x).size()
typedef long long ll;
typedef pair<int, int> pii;
typedef vector<int> vi;

struct Counter {
	vector<short> xs, ys, ups;
	void insert(short x, short y, short len) {
		xs.push_back(x);
		ys.push_back(y);
		ups.push_back((short)(x + y + len));
	}
	int count(short x, short y) {
		// int res = sz(xs);
		unsigned short res = (unsigned short)sz(xs);
		rep(i,0,sz(xs)) {
			// if (ys[i] <= y && xs[i] <= x && x + y <= ups[i]) res++;
			// res += (ys[i] <= y && xs[i] <= x && x + y <= ups[i]);
			// res += ((ys[i] <= y) & (xs[i] <= x) & (x + y <= ups[i]));
			// res -= ((ys[i] > y) | (xs[i] > x) | (x + y > ups[i]));
			// res -= ((ys[i] > y) | (xs[i] > x) | ((short)(x + y) > ups[i]));
			res = (unsigned short)(res - ((ys[i] > y) | (xs[i] > x) | ((short)(x + y) > ups[i])));
		}
		return res;
	}
};

int main() {
	cin.sync_with_stdio(false);
	int N, Q;
	cin >> N >> Q;
	Counter d1, d2, d3, d4;
	rep(qi,0,Q) {
		int type;
		cin >> type;
		if (type == 1) {
			short dir, x, y, len;
			cin >> dir >> x >> y >> len;
			if (dir == 1) d1.insert((short)x, (short)y, len);
			if (dir == 2) d2.insert((short)x, (short)-y, len);
			if (dir == 3) d3.insert((short)-x, (short)y, len);
			if (dir == 4) d4.insert((short)-x, (short)-y, len);
		}
		else {
			assert(type == 2);
			short x, y;
			cin >> x >> y;
			int ans = 0;
			ans += d1.count((short)x, (short)y);
			ans += d2.count((short)x, (short)-y);
			ans += d3.count((short)-x, (short)y);
			ans += d4.count((short)-x, (short)-y);
			cout << ans << '\n';
		}
	}
}
