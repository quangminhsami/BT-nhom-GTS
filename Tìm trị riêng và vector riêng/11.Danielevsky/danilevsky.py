from mydata import *
import numpy as np
def specialCase1(A, k, j):
    copA = A.copy()
    copA[[k-1, j], :] = copA[[j, k-1], :]
    copA[:, [k-1, j]] = copA[:, [j, k-1]]
    return copA

def specialCase2(A, k):
    return A[k:, k:], A[0:k, 0:k]

def findSimpleA(A, k):
    n, _ = A.shape
    M = np.eye(n)
    M[k-1, :] = A[k, :]

    inverseM = np.eye(n)
    inverseM[k-1, :] = -A[k, :] / A[k, k-1]
    inverseM[k-1, k-1] = 1 / A[k, k-1]

    similarA = M.dot(A).dot(inverseM)
    return similarA, M, inverseM

def getCharPolynomial(A):

    n = A.shape[0]
    p = (-1)**n * np.ones(n+1)
    p[1:] = p[1: ] * (-1) *A[0, :]
    return p

def mulPoly(p):

    res = p[0]
    for i in range(1, len(p)):
        res = np.polymul(res, p[i])
    return res

def Danilevski(A):

    all_polynomial = []
    n = A.shape[0]
    back = np.eye(n)
    similar = A.copy()
    list_eigenvalues = []
    list_eigenvectors = []
    charFunc = []
    k = n-1
    m = n    
    while k > 0:
        if similar[k, k-1] != 0:
 
            similar, _, inverseM = findSimpleA(similar, k)
            back = back.dot(inverseM)
            
        else:
            non = False
            for j in range(0, k-1):
    
                if similar[k, j] != 0:
                    similar = specialCase1(similar, k , j)
                    
                    back[:, [j, k-1]] = back[:, [k-1, j]]
                    
                    non = True
                    k = k+1
                    break
                    
                    
            if not non:

                for j in range(k, m-1):
                    M = np.eye(m)
                    M[:k, j+1] = -similar[:k, j]
            
                    inverseM = np.eye(m)
                    inverseM[:k, j+1] = similar[:k, j]
                    similar = M.dot(similar).dot(inverseM)
                    back = back.dot(inverseM)
                    
                tt = False
         
                for j in range(k-1, -1, -1):
                    if similar[j, m-1] != 0:
                        M = np.eye(m)
                        x = M[k:m, :]
                        y = M[k-1, :]
                        M = np.vstack([M[0:k-1, :], x, y])
                        M1 = M.T
                        similar = M.dot(similar).dot(M1)
                        back = back.dot(M1)
                        k = m;
                        tt = True
                        break
                
                if not tt:
                    X = similar[k:, k:]
            
                    t = X.shape[0]
                    eigenValues = findValue(X)
                    for j in range(len(eigenValues)):
                        print("Vector rieng cua: ", eigenValues[j])
                        list_eigenvalues.append(eigenValues[j])
                        y = np.power(eigenValues[j], np.arange(t))[::-1].reshape((t, 1))
                        v = np.zeros((n, 1))
                        p = np.zeros((m, 1))
                        p[k:m] = y
                        p = back.dot(p)
                        v[:m] = p
                        list_eigenvectors.append(v)
                        print(v)
                    m = k
                    similar = similar[:k, :k]
                    back = np.eye(m)
#        
        k = k - 1
    X = similar
    t = X.shape[0]
    eigenValues = findValue(X)
    for j in range(len(eigenValues)):
        print("Vector rieng cua: ", eigenValues[j])
        list_eigenvalues.append(eigenValues[j])
        y = np.power(eigenValues[j], np.arange(t))[::-1].reshape((t, 1))
        p = np.zeros((m, 1))
        v = np.zeros((n, 1))
        p[k:m] = y
        p = back.dot(p)
        v[:m] = p
        list_eigenvectors.append(v)
        print(v)
    return list_eigenvalues, list_eigenvectors
def findValue(A):

    p = getCharPolynomial(A)
    p = list(p)[::-1]
    eigenValues = findRoots(p)
    return sorted(eigenValues)
	
# A = np.array([[1, 2, 3, 4], [2, 1, 2, 3], [3, 2, 1, 2], [4, 3, 2, 1]])
# A = np.loadtxt('B:/BT nhóm GTS/Tìm trị riêng và vector riêng/11.Danielevsky/TH1.txt', dtype = float)
# A = np.loadtxt('B:/BT nhóm GTS/Tìm trị riêng và vector riêng/11.Danielevsky/TH2.txt', dtype = float)
A = np.loadtxt('B:/BT nhóm GTS/Tìm trị riêng và vector riêng/11.Danielevsky/TH3.txt', dtype = float)
print("Ma trận A : \n", A)
eigenvalue, eigenvector = Danilevski(A)