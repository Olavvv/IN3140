import sympy as sp
import numpy as np
from sympy import cos, sin, Symbol, symbols, diff
from sympy.matrices import Matrix
from sympy.physics.vector import dynamicsymbols, ReferenceFrame, time_derivative

sp.init_printing(use_unicode=True)

# Define weight of masses
m1 = 0.3833
m2 = 0.2724
m3 = 0.1406

dot = '\u0307'

# Define symbols for the link lengths and theta
L1 = Symbol("L1")
L2 = Symbol("L2")
L3 = Symbol("L3")
th1 = dynamicsymbols("theta1", real=True)
th2 = dynamicsymbols("theta2", real=True)
th3 = dynamicsymbols("theta3", real=True)

q = Matrix([th1,
            th2,
            th3])
qdot = Matrix([diff(th1,Symbol('t')),
               diff(th2,Symbol('t')),
               diff(th3,Symbol('t'))])

sp.pprint(qdot)
#Column vector containing gravitational constant
g = Matrix([0, 
            0, 
            9.81])

# Define the inertia tensors for each mass.
I1x, I1y, I1z, I2x, I2y, I2z, I3x, I3y, I3z = symbols('I1x I1y I1z I2x I2y I2z I3x I3y I3z', positive=True)

I1 = Matrix([[I1x, 0, 0],           
            [0, I1y, 0], 
            [0, 0, I1z]])

I2 = Matrix([[I2x, 0, 0], 
            [0, I2y, 0], 
            [0, 0, I2z]])

I3 = Matrix([[I3x, 0, 0], 
            [0, I3y, 0], 
            [0, 0, I3z]])

# Define vectors to height of mass.
rc1 = Matrix([0,
              0, 
              L1/2])

rc2 = Matrix([0,
              0,
              (L2*sin(th2) + 2*L1)/2])

rc3 = Matrix([0,
              0,
              (L2*sin(th2) + 2*L1)/2])

J = Matrix([[-sin(th1)*(L2*cos(th2) + L3*cos(th2+th3)), -cos(th1)*(L2*sin(th2) + L3*sin(th2+th3)), -cos(th1)*(L3*sin(th2+th3))],
            [cos(th1)*(L2*cos(th2) + L3*cos(th2+th3)), -sin(th1)*(L2*sin(th2) + L3*sin(th2+th3)), -sin(th1)*(L3*sin(th2+th3))],
            [0, (L2*cos(th2) + L3*cos(th2+th3)), L3*cos(th2+th3)],
            [0, sin(th1), sin(th1)],
            [0, -cos(th1), -cos(th1)],
            [1, 0, 0]])

#Define the Jvi and Jw matrices
Jv1 = Matrix([[0, 0, 0],
              [0, 0, 0],
              [0, 0, 0]])
Jv2 = Matrix([[-sin(th1)*(L2*cos(th2)), -cos(th1)*(L2*sin(th2)), 0],
              [cos(th1)*(L2*cos(th2)), -sin(th1)*(L2*sin(th2)), 0],
              [0, (L2*cos(th2)), 0]])
Jv3 = Matrix([[-sin(th1)*(L2*cos(th2) + L3*cos(th2+th3)), -cos(th1)*(L2*sin(th2) + L3*sin(th2+th3)), -cos(th1)*(L3*sin(th2+th3))],
            [cos(th1)*(L2*cos(th2) + L3*cos(th2+th3)), -sin(th1)*(L2*sin(th2) + L3*sin(th2+th3)), -sin(th1)*(L3*sin(th2+th3))],
            [0, (L2*cos(th2) + L3*cos(th2+th3)), L3*cos(th2+th3)]])

Jw = Matrix([[0, sin(th1), sin(th1)],
              [0, -cos(th1), -cos(th1)],
              [1, 0, 0]])

# Calculate the rotational matrix from base frame to all mass centers.
RM1 = Matrix([[1, 0, 0],
              [0, 1, 0],
              [0, 0, 1]])

RM2 = RM1*Matrix([[cos(th2), -sin(th2), 0],
                  [sin(th2), cos(th2), 0],
                  [0, 0, 1]])

RM3 = RM2*Matrix([[cos(th3), -sin(th3), 0],
                  [sin(th3), cos(th3), 0],
                  [0, 0, 1]])

# Calculate the potential energy for the masses.
P1 = m1*g.T*rc1
P2 = m2*g.T*rc2
P3 = m2*g.T*rc3
P = P1+P2+P3

print("The potential energy for the dust crawler: \n")
sp.pprint(sp.simplify(P))

# Kinetic energy terms
K1 = (m1*Jv1.T*Jv1) + (Jw.T*RM1*I1*RM1.T*Jw)
K2 = (m1*Jv2.T*Jv2) + (Jw.T*RM2*I2*RM2.T*Jw)
K3 = (m1*Jv3.T*Jv3) + (Jw.T*RM3*I3*RM3.T*Jw)

Dq = K1+K2+K3

# The finished kinetic energy for the robot.
K = 0.5*qdot.T*Dq*qdot 

# Calculate the coriolis/sentripedal matrix
def SumCkj(k,j):
    Ckj = 0
    for i in range(3): # Implementation of  equation (7.60) in book
        Ckj += 0.5*(diff(Dq[k,j],q[j]) + diff(Dq[k,i], q[j]) - diff(Dq[i,j],q[k]))
    return Ckj

C = Matrix([[SumCkj(0,0), SumCkj(0,1), SumCkj(0,2)],
            [SumCkj(1,0), SumCkj(1,1), SumCkj(1,2)],
            [SumCkj(2,0), SumCkj(2,1), SumCkj(2,2)]])



# Calculate the gravitational matrix.
g = Matrix([diff(P,th1),
            diff(P,th2),
            diff(P,th3)])

# The complete dynamic model
tau = Dq + C + g
