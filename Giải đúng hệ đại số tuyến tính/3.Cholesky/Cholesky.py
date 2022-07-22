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

# gói kiểm tra ma trận xác định dương
def positive_definite_matrix(A):
    if np.all(np.linalg.eigvals(A) > 0) :
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

# Gói giải phương trình ma trận tam giác dưới Y = S*X   ST*Y = B1
def solve1(S,B1):
    n = S.shape[0]
    sum = 0
    Y = np.zeros(n, dtype = float)
    for i in range(n):
        if S[i,i] == 0:
            print("ko thỏa mãn điều kiện cholesky")
        sum = 0
        for k in range(i):
            sum += S[k,i] * Y[k]
        Y[i] = (B1[i] - sum) / S[i,i]
    return Y

# Gói giải phương trình ma trận tam giác trên S*X=Y
def solve2(S,Y):
    n = S.shape[0]
    X = np.zeros(n)
    sum = 0
    for i in range(n):
        sum = 0
        for k in range(i):
            sum += S[n-i-1, n-k-1] * X[n-k-1]
        
        X[n-i-1] = (Y[n-i-1] - sum) / S[n-i-1, n-i-1]
    return X


# =========================================================================
# chương trình chính

# Ma trận đối xứng:
# A = np.loadtxt('B:/BT nhóm GTS/Giải đúng hệ đại số tuyến tính/3.Cholesky/cholesky_doi_xung.txt', dtype = float)


# Ma trận ko đối xứng:
A = np.loadtxt('B:/BT nhóm GTS/Giải đúng hệ đại số tuyến tính/3.Cholesky/cholesky_ko_doi_xung.txt', dtype = float)
print("A = \n", A)

b = np.array([9.45, -12.20, 7.78, -8.1, 10.0])

if isSymetric(A) and positive_definite_matrix(A):
    print("Ma trận A đối xứng \n")
    A1 = A 
    print("A1 = \n", A1)
    B1 = b
    print("B1 = ", B1)   
else:
    print("Ma trận A ko đối xứng \n")
    A1 = ATA(A)
    print("A1 = \n", A1)
    B1 = ATb(A,b)
    print("B1 = ", B1)   

S = cholesky(A1)
print("S = \n", S)

Y = solve1(S, B1)
print("Y = ", Y)

X = solve2(S,Y)
print("Nghiệm của hệ PT là: ")
print("X = ", X)

# nghiệm theo numpy
# print("Nghiệm theo numpy:")
# print(np.linalg.solve(A,b))


