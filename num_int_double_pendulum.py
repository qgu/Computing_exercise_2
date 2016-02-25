import numpy as np
from scipy.integrate import odeint 
import matplotlib.pyplot as plt
from numpy import sin, cos
# Defining constants
m1 = 1
m2 = 1
g = 10
l1 = 1
l2 = 1

def evolution(y, t):
    dydt = np.zeros(4)
    dydt[0] = y[2]
    dydt[1] = y[3]

    cos_y1_y2 = cos(y[0] - y[1])
    sin_y1_y2 = sin(y[0] - y[1])
    sin_y1    = sin(y[0])
    sin_y2    = sin(y[1])
    
    # Determinant
    D = 1#(m1 + m2) * l1 * l2 - m2 * l1 * l2 * cos_y1_y2 ** 2
    
    dydt[2] = (- l2**2 * m2 * sin_y1_y2 * y[3]**2 - ( m1 + m2 ) * l2 * g * sin_y1 - m2 * l1 * l2 * y[2] **2 * sin_y1_y2 * cos_y1_y2 + m2 * l2 * g * sin_y2 * cos_y1_y2)/D
    dydt[3] = (l1 * l2 * m2 * sin_y1_y2 * cos_y1_y2 + ( m1 + m2 ) * l1 * g * sin_y1 * cos_y1_y2 + ( m1 + m2 ) * l1**2 * y[2] **2 * sin_y1_y2 - g * ( m1 + m2 ) * l1 * sin_y2)/D
    
    return dydt

def test_evolution(evolution):
    y0_test = [0,0,1,1]
    t = 0;
    y_out = evolution(y0_test, t)
    y_out_theory = [1,1,0,0]
    for ii in range(0,len(y_out_theory)):
        if np.abs(y_out[ii] - y_out_theory[ii]) > 1e-16 :
            print 'The' + str(ii) + ' component of the evolution function is wrong.'

test_evolution(evolution)

t = np.linspace(0,8,100)
y0 = [0.01, 0, .01, 0]
sol = odeint(evolution, y0, t)

plt.plot(t,sol[:,0])
plt.show()



