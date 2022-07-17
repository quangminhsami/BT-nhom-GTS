# Eigenvalues using Power iteration

from cv2 import eigen
import numpy as np 
import pandas as pd

def relativeError(xnew, xold):
    return abs((xnew - xold) / xnew) * 100

result = pd.DataFrame(columns =['eigenvalue', 'error'])

A = np.array([[1,9,8],
              [2,7,6],
              [7,20,9]])

x = np.array([4,2,-8])
max_iterations = 10
e_tolerance = 0.5
oldEigenValue = 0

for i in range(max_iterations):
    x = np.dot(A,x)
    eigenValue = np.linalg.norm(x)
    error = relativeError(eigenValue, oldEigenValue)
    x = x / eigenValue
    result.loc[i] = [eigenValue, error]
    if error < e_tolerance:
        oldEigenValue = eigenValue
print(result)