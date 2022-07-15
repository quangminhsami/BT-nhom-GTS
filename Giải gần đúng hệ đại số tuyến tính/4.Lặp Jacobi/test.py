import numpy as np

a = np.array([[9,-1,2,3],
              [1,8,-1,1],
              [-2,4,10,-1],
              [1,2,3,7]], dtype = float)

print(np.linalg.norm(a, 1))