import random
n = 5000
q = 10**5
print(n, q)
for i in range(q):
    if i < q//2:
        dir = random.randint(1, 4)
        le = random.randint(1, n)
        x = random.randint(1, n - le + 1)
        y = random.randint(1, n - le + 1)
        if dir >= 3:
            x += le
        if dir % 2 == 0:
            y += le
        print(1,dir,x,y,le)
    else:
        x = random.randint(1, n)
        y = random.randint(1, n)
        print(2,x,y)
