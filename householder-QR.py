import numpy as np
from numpy import linalg as la
# A = np.array(([-1.0, -1.0, 1.0], [1.0, 3.0, 3.0], [-1.0, -1.0, 5.0], [1.0, 3.0, 7.0]))   # m > n
# A = np.array(([1.0, -4.0], [2.0, 3.0], [2.0, 2.0]))                                      # m > n
# A = np.array(([1.0, 1.0, 2.0], [1.0, 2.0, 1.0]))                                           # m < n
A = np.array(([4.0, 3.0, 0.0], [2.0, 1.0, 2.0], [4.0, 7.0, 0.0]))                        # m = n
print("A\n", A)
m, n = np.shape(A)
s = min(m-1, n)
u = np.zeros(m)
beta = 0
alpha = 0
R = np.zeros((m, n))
Q = np.eye(m)
for k in range(0, s):
    v = np.zeros(m)
    v[k:m] = A[k:m, k]
    alpha = - np.sign(A[k, k]) * la.norm(v)
    v[k] = v[k] - alpha                         # w = v - alpha
    u[k] = v[k]
    print("\nu:", u)
    A[k, k] = alpha
    print("\nv:", v)
    for i in range(k+1, m):
        A[i, k] = v[i]
    for j in range(k+1, n):
        s = 0
        beta = 2/(la.norm(v)**2)
        for t in range(k, m):
            s += v[t] * A[t, j]
        beta = beta * s
        print("beta", beta)
        for i in range(k, m):
            ff = A[i, j] - beta * v[i]
            A[i, j] = ff
    if k == m-2:
        print("\nMatrix A in The Last Step Is:\n", A)
        R = np.triu(A)
        print("\nTHe Upper Triangular Matrix R Is:\n", R)
    else:
        print("\nMatrix A After Step ", k+1, "Is:\n", A)
    Q = Q.dot(np.eye(m) - (2/(la.norm(v)**2))*np.outer(v, v))

if m < n:
    S = A[:, n-m+1:]
    print("Rows Are Less Than Columns:\n Matrix S(", m, "by", n-m, ") Is:\n", S)

print("Q\n", Q)
print("A=QR\n", np.dot(Q, R))
