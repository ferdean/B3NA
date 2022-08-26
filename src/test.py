import sympy as sm
import numpy as np

xi = sm.Symbol('xi')
#x  = (x1 - x0) * xi + x0

phi1 = 1 - 3 * xi**2 + 2 * xi**3
phi2 = xi * (xi - 1)**2
phi3 = 3 * xi**2 - 2 * xi**3
phi4 = xi**2 * (xi - 1)
phi  = [phi1, phi2, phi3, phi4]

M = np.zeros((4,4))
S = np.zeros((4,4))

for i in range(4):
    for j in range(4):
    
        S[i,j] = float(sm.integrate(sm.diff(phi[i], xi, 2)*sm.diff(phi[j], xi ,2), (xi, 0, 1)))
                    
        M[i,j] = float(sm.integrate(phi[i]*phi[j],(xi, 0, 1))) 
                        
print(M)
print(S)