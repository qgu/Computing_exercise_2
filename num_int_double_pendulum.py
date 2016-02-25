import numpy as np
import scipy.integrate.odeint as odeint

# Defining constants
m1 = 1
m2 = 1
g = 10
l1 = 1
l2 = 1
# Matrix G
G = [[(m1 + m2)*l1, m2*l2], [l1, l2]]

def evolution(y0):
    y1 = np.zeros(4)
    K = G
    y1[0] = y0[2]
    y1[1] = y0[3]
    # Matrix K
    cos_y1_y2 = np.cos(y0[0] - y0[1])
    sin_y1_y2 = np.sin(y0[0] - y0[1])
    K[0][1] = G[0][1]*cos_y1_y2
    K[1][0] = G[1][0]*cos_y1_y2
    W = np.matrix(K)
    [y1[2], y1[3]] = W.I * [-m2 * l2 * sin_y1_y2 * y0[3] **2 - ( m1 + m2 ) * g * np.sin(y0[0]),
           l1 * y0[2] **2 * np.sin(y0[0] - y0[1]) - g * np.sin(y0[1])] 
    
    return y1

    



