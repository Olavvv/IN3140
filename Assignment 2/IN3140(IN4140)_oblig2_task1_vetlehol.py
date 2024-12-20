import numpy as np
from numpy import cos, sin, arctan
# Declaring variables 

L1 = 100.9 #[mm] Length of link 1.
L2 = 222.1 #[mm] Length of link 2.
L3 = 136.2 #[mm] Length of link 3.




def forward(joint_angles):
    """
    The forward kinematics function takes x sets of joint angles as input,
    and gives the corresponding cartesian coordinates for the tip of the arm
    as output, takes in joint anglesas: [theta1, theta2, theta3], returns as array: [x,y,z]
    """
    th1 = joint_angles[0]
    th2 = joint_angles[1]
    th3 = joint_angles[2]
    x = cos(th1)*( (L2*cos(th2)) + (L3*cos(th2 + th3)) )
    y = sin(th1)*( (L2*cos(th2)) + (L3*cos(th2 + th3)) )
    z = (L1*L2*sin(th2)) + (L3*sin(th2 + th3))
    return np.array([x,y,z]) 

def inverse(cart_cord):
    """
    The inverse kinematics function takes the cartesian position of the tip
    of the pen as input, and gives the corresponding joint configurations as
    output
    """
    x = cart_cord[0]
    y = cart_cord[1]
    z = cart_cord[2]
    
    th1 = arctan(x/y)
    th2 = y
    th3 = z
    return [th1, th2, th3]



point = forward([90, -30, 45])
angles = inverse(point)

print(f"Angles fed into forward function: {[90, -30, 45]}, point given by forward function: {point}. Which is fed into inverse function and got these angles: {angles}")