import numpy as np
# A = np.array(([0.0, 1.0, 1.0], [1.0, 2.0, 3.0], [1.0, 1.0, 1.0]))
# A = np.array(([1.0, 2, 3, 5], [1, 2, 1, 1], [2, 5, 9, 7], [4, 3, 6, 1]))
# A = np.array(([1.0, 2, 3], [7, 8, 3], [2, 4, 1]))
A = np.array(([5.0, 2, 3, 1, 2], [4, 2, 8, 4, 5], [9, 2, 3, 4, 8], [1, 4, 1, 4, 5], [3, 2, 3, 4, 5]))
print("A\n", A)
m, n = np.shape(A)
s = np.zeros(m)
c = np.zeros(m)
ll = np.zeros(m)
k = np.zeros(m)
m1 = int(m*(m-1)/2)
j = np.zeros((m1, m, m))
i = 0
x = 0
Q = np.zeros((m-1, m, m))
for p in range(0, n-1):
    for q in range(p+1, m):
        j[i, :, :] = np.eye(m)
        if A[q, p] >= A[p, p]:
            t = A[p, p] / A[q, p]
            s[p] = 1/np.sqrt(1 + t ** 2)
            c[p] = s[p] * t
        else:
            t = A[q, p] / A[p, p]
            c[p] = 1/np.sqrt(1 + t ** 2)
            s[p] = c[p] * t
        ll[p] = q
        k[p] = p
        j[i, p, p] = c[p]
        j[i, p, q] = s[p]
        j[i, q, p] = -s[p]
        j[i, q, q] = c[p]
        print("\nj[", p, ",", q, "] In Step ", p+1, ":\n", j[i, :, :])
        A = np.dot(j[i, :, :], A)
        print("i", i)
        i += 1
    print("\nA After Step ", p+1, "Is:\n", A)
    Q[p, :, :] = np.eye(m)
    if p == 0:
        for t in range(i - 1, p-1, -1):
            print("t", t)
            Q[p, :, :] = Q[p, :, :].dot(j[t, :, :])
    elif x == i-1:
        Q[p, :, :] = Q[p, :, :].dot(j[i-1, :, :])
    else:
        for t in range(i-1, x, -1):
            print("x", x)
            print("t", t)
            Q[p, :, :] = Q[p, :, :].dot(j[t, :, :])
    x += n-2
    print("\nThe Orthogonal Matrix Q In Step", p+1,  "Is:\n", Q[p, :, :])
QQ = np.eye(m)
for i in range(m-2, -1, -1):
    QQ = np.dot(QQ, Q[i, :, :])
QQ = np.transpose(QQ)
print("\nThe Upper Triangular Matrix R Is:\n", A)
print("\nThe Final Orthogonal Matrix Q Is:\n ", QQ)
print("A = QR\n", np.dot(QQ, A))
