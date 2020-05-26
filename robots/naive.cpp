// 0m10.622s; 0m4.245s with && -> &
// -O3 and target pragma doesn't help; no auto-vectorization happens here
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
	struct Tri {
		int x, y, ups;
	};
	vector<Tri> tris;
	void insert(int x, int y, int len) {
		tris.push_back({x, y, x + y + len});
	}
	int count(int x, int y) {
		int res = 0;
		rep(i,0,sz(tris)) {
			if (tris[i].y <= y && tris[i].x <= x && x + y <= tris[i].ups) res++;
			// res += (tris[i].y <= y && tris[i].x <= x && x + y <= tris[i].ups);
			// res += ((tris[i].y <= y) & (tris[i].x <= x) & (x + y <= tris[i].ups));
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
			int dir, x, y, len;
			cin >> dir >> x >> y >> len;
			if (dir == 1) d1.insert(x, y, len);
			if (dir == 2) d2.insert(x, -y, len);
			if (dir == 3) d3.insert(-x, y, len);
			if (dir == 4) d4.insert(-x, -y, len);
		}
		else {
			assert(type == 2);
			short x, y;
			cin >> x >> y;
			int ans = 0;
			ans += d1.count(x, y);
			ans += d2.count(x, -y);
			ans += d3.count(-x, y);
			ans += d4.count(-x, -y);
			cout << ans << '\n';
		}
	}
}
