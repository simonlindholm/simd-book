import random
import sys

random.seed(sys.argv[1])

n = 10**6 # random.randint(1, 600)
q = 10**6
print(n, q)
print(' '.join(random.choice("01") for _ in range(n)))
for qi in range(q):
    tp = random.randint(1,7)
    lo = 2 if tp in [4,6] else 1
    hi = n-1 if tp in [3,5] else n
    a = random.randint(lo, hi)
    b = random.randint(lo, hi)
    if a > b:
        a,b = b,a
    if tp in [4,6]:
        a -= 1
    if tp in [3,5]:
        b += 1
    print(tp,a,b)
