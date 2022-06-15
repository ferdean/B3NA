import sympy as sm
import numpy as np
def get_local_matrix():
    # <><><><><><>
    x = sm.Symbol('x')
    print("done")
    S = np.zeros((4,4))
    phi1 = 1 - 3* x**2 + 2 * x**3
    phi2 = x* (x-1)**2
    phi3 = 3*x**2 - 2*x**3
    phi4 = x**2 * (x-1)
    phi  = [phi1,phi2,phi3,phi4]
    for i in range(4):
        for j in range(4):
            S[i,j] = float(sm.integrate(sm.diff(phi[i],x,2)*sm.diff(phi[j],x,2),(x,0,1)))
    #a = 
    #a = sm.diff(phi1,x,4)
    print(S)
    return 0
get_local_matrix()
