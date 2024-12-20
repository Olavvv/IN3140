import numpy as np
from numpy import cos, sin
# Declaring variables 

L1 = 100.9 #[mm] Length of link 1.
L2 = 222.1 #[mm] Length of link 2.
L3 = 136.2 #[mm] Length of link 3.


# For this function, im assuming, we only want the linear velocities.

def jacobian(joint_angles, joint_velocities):
    """
    It takes the instant joint angles and joint velocities as input, and gives a
    3-dimensional vector of cartesian velocities of the tip of the pen as output
    """
    th1 = joint_angles[0]
    th2 = joint_angles[1]
    th3 = joint_angles[2]
    q1 = joint_velocities[0]
    q2 = joint_velocities[1]
    q3 = joint_velocities[2]

    j_v = [[sin(th1)*(L2*cos(th2) + L3*cos(th2+th3)), -cos(th1)*((L1*L2*sin(th2)) + (L3*sin(th2+th3)) - L1), -cos(th1)*((L1*L2*sin(th2)) + (L3*sin(th2+th3)) - (L2*sin(th2)) - L1)],
           [cos(th1)*(L2*cos(th2) + L3*cos(th2+th3)), -sin(th1)*((L1*L2*sin(th2)) + (L3*sin(th2+th3)) - L1), -sin(th1)*((L1*L2*sin(th2)) + (L3*sin(th2+th3)) - (L2*sin(th2)) - L1)],
           [0, -((cos(th1))**2)*((L2*cos(th2)) + (L3*cos(th2+th3))), -((L2*cos(th2)) + (L3*cos(th2+th3))) + (L2*cos(th2))]]

    vel_x = (j_v[0][0]*q1) + (j_v[0][1]*q2) + (j_v[0][2]*q3)
    vel_y = (j_v[1][0]*q1) + (j_v[1][1]*q2) + (j_v[1][2]*q3)
    vel_z = (j_v[2][0]*q1) + (j_v[2][1]*q2) + (j_v[2][2]*q3)

    return [vel_x, vel_y, vel_z]


jnt_angls = [90, -30, 45]
jnt_vels = [.1, .05, .05] # rad/s

print(jacobian(jnt_angls, jnt_vels))