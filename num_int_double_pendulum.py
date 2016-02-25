import numpy as np
from scipy.integrate import odeint 
import matplotlib.pyplot as plt
from numpy import sin, cos
# Defining constants
m1 = 0.000001
m2 = 1.0
g = 10.0
l1 = 1.0
l2 = 1.0

def evolution(y, t):
    dydt = np.zeros(4)
    dydt[0] = y[2]
    dydt[1] = y[3]

    cos_y1_y2 = cos(y[0] - y[1])
    sin_y1_y2 = sin(y[0] - y[1])
    sin_y1    = sin(y[0])
    sin_y2    = sin(y[1])
    
    # Determinant
    D = (m1 + m2) * l1 * l2 - m2 * l1 * l2 * cos_y1_y2 ** 2
    
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

def energy(y):
    x2dot = l1 * cos( y[:,0] ) * y[:,2] + l2 * cos( y[:,1] ) * y[:,3]
    y2dot = l1 * sin( y[:,0] ) * y[:,2] + l2 * sin( y[:,1] ) * y[:,3]
    
    T1 = m1 * ( l1 * y[:,2] ) **2 /2.0
    T2 = m2 * ( x2dot**2 + y2dot**2) / 2.0
    V1 = m1 * g * ( l1 * ( 1 - cos(y[:,0]) ) ) 
    V2 = m2 * g * ( l1 * ( 1 - cos(y[:,0]) ) + l2 * ( 1 - cos(y[:,1]) ) )

    return (T1 + T2 + V1 + V2)

t = np.linspace(0,4,1000)
y0 = [0.0, 0, 0.0, 1]

E0 = energy(np.array([y0]))

y = odeint(evolution, y0, t)

E = [E0 for i in t]

plt.plot(t, y[:,3])

plt.show()



