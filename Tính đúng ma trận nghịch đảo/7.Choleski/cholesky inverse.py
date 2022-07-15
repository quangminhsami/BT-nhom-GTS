import numpy as np 
import math

# Gói kiểm tra ma trận A có đối xứng hay không?
def isSymetric(A):
    n = A.shape[0]
    for i in range(n):
        for j in range(n):
            if A[i,j] != A[j,i]:
                return False
    return True

# Gói nhân AT với A 
def ATA(A):
    A1 = A.T @ A
    return A1

# Gói nhân AT với b
def ATb(A,b):
    B1 = A.T @ b
    return B1

# Gói tìm ma trận S theo phân tích Cholesky ST*S=A1
def cholesky(A1):
    n = A1.shape[0]
    sum = 0
    S = np.arange(n**2)
    S = S.reshape((n,n))
    S = np.zeros_like(S, dtype = float)

    for i in range(n):
        for j in range(i+1):
            sum = 0
            for k in range(j):
                sum += S[k,i] * S[k,j]
            if sum == A1[i,i]:
                break
            if i == j:
                S[j,i] = math.sqrt(A1[i,i] - sum)
            else:
                S[j,i] = (A1[j,i] - sum) / (S[j,j])
        if S[i,i] == 0:
            break 
    return S

def cholesky_inverse(S):
    n = S.shape[0]
    S_1 = np.arange(n**2)
    S_1 = S_1.reshape((n,n))
    S_1 = np.zeros_like(S_1, dtype = float)

    for i in range(n-1, -1 , -1):
        S_1[i,i] = 1 / S[i,i]
        for j in range(i, n):
            sum = 0
            for k in range(i+1, n):
                sum += S[i,k] * S_1[k,j]
                if j > i:
                    S_1[i,j] = - sum / S[i,i]

    return S_1

# =========================================================================
# chương trình chính
A = np.loadtxt('B:/BT nhóm GTS/Tính đúng ma trận nghịch đảo/7.Choleski/cholesky inverse det equal 0.txt', dtype = float)
print("A = \n", A)

if np.linalg.det(A) != 0:

    if isSymetric(A):
        A1 = A 
        print("A1 = \n", A1) 
    else:
        A1 = ATA(A)
        print("A1 = \n", A1)

    S = cholesky(A1)
    print("S = \n", S)

    S_1 = cholesky_inverse(S)
    print("S* = \n", S_1)

    A1_1 = S_1 @ S_1.T

    print("Nghịch đảo của ma trận A là: \n", A1_1 @ A.T)
else:
    print("Ma trận ko khả nghịch")




