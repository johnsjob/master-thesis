from numpy import array as mat, sin, arccos, dot, linspace
from numpy.linalg import norm

def rot_to_quat(R):

    r11, r12, r13,\
    r21, r22, r23,\
    r31, r32, r33 = R.reshape(9)

    q0 = ( r11 + r22 + r33 + 1.0) / 4.0
    q1 = ( r11 - r22 - r33 + 1.0) / 4.0
    q2 = (-r11 + r22 - r33 + 1.0) / 4.0
    q3 = (-r11 - r22 + r33 + 1.0) / 4.0
    if(q0 < 0.0):
        q0 = 0.0
    if(q1 < 0.0):
        q1 = 0.0
    if(q2 < 0.0):
        q2 = 0.0
    if(q3 < 0.0):
        q3 = 0.0
    q0 = sqrt(q0)
    q1 = sqrt(q1)
    q2 = sqrt(q2)
    q3 = sqrt(q3)
    if((q0 >= q1) and (q0 >= q2) and (q0 >= q3)):
        q0 *= +1.0
        q1 *= sign(r32 - r23)
        q2 *= sign(r13 - r31)
        q3 *= sign(r21 - r12)
    elif((q1 >= q0) and (q1 >= q2) and (q1 >= q3)):
        q0 *= sign(r32 - r23)
        q1 *= +1.0
        q2 *= sign(r21 + r12)
        q3 *= sign(r13 + r31)
    elif((q2 >= q0) and (q2 >= q1) and (q2 >= q3)):
        q0 *= sign(r13 - r31);
        q1 *= sign(r21 + r12);
        q2 *= +1.0
        q3 *= sign(r32 + r23);
    elif(q3 >= q0 and q3 >= q1 and q3 >= q2):
        q0 *= sign(r21 - r12);
        q1 *= sign(r31 + r13);
        q2 *= sign(r32 + r23);
        q3 *= +1.0
    else:
        print("coding error\n")

    q = mat( (q0, q1, q2, q3) )
    q = q / norm(q)

    return q

def quat_to_rot(q):

    q0, q1, q2, q3 = q

    r1 = [
        q0**2 + q1**2 - q2**2 - q3**2,
        2*(q1*q2 + q0*q3),
        2*(q1*q3 - q0*q2)
        ]

    r2 = [
        2*(q1*q2 - q0*q3),
        q0**2 - q1**2 + q2**2 - q3**2,
        2*(q2*q3 + q0*q1)
        ]

    r3 = [
        2*(q1*q3 + q0*q2),
        2*(q2*q3 - q0*q1),
        q0**2 - q1**2 - q2**2 + q3**2
        ]

    R = mat(zip(r1,r2,r3))
    return R

def quat_conj(q):
    return mat((q[0], -q[1], -q[2], -q[3]))

def quat_unit(q):
    return (mat(q) / norm(q))

def quat_inv(q):
    return(quat_conj(q) / (norm(q)**2))

def quat_mul(q1, q2):
    """
    q1*q2 = (w1,x1,y1,z1)*(w2,x2,y2,z2)
          = ( w1*w2-x1*x2-y1*y2-z1*z2 ,
              w1*x2+x1*w2+y1*z2-z1*y2 ,
              w1*y2+y1*w2+z1*x2-x1*z2 ,
              w1*z2+z1*w2+x1*y2-y1*x2  )
    """
    result = ( q1[0] * q2[0] - q1[1] * q2[1] - q1[2] * q2[2] - q1[3] * q2[3],
               q1[0] * q2[1] + q1[1] * q2[0] + q1[2] * q2[3] - q1[3] * q2[2],
               q1[0] * q2[2] + q1[2] * q2[0] + q1[3] * q2[1] - q1[1] * q2[3],
               q1[0] * q2[3] + q1[3] * q2[0] + q1[1] * q2[2] - q1[2] * q2[1] )
    return mat(result)

def quat_rot_v(q,v):
    v = mat([0] + list(v))
    products = (quat_inv(q), v , q)
    result = reduce(quat_mul, products)
    return mat(result)

def slerp(q1, q2, u):
    """
    u : [0, ..., 1]
    """
    theta = arccos(dot(q1, q2))
    V = sin((1-u)*theta) / sin(theta)
    W = sin(u*theta) / sin(theta)
    return mat(V*mat(q1) + W*mat(q2))

def lerp(p1, p2, u):
    """
    u : [0, ..., 1]
    """
    return mat(p1)*(1-u) + u*mat(p2)

def point_lerp(p1, p2, num_points=50):
    interp = linspace(0, 1, num_points)
    result = mat(map(lambda u: lerp(p1, p2, u), interp))
    return result

def quat_slerp(q1, q2, num_points=50):
    interp = linspace(0, 1, num_points)
    result = mat(map(lambda u: slerp(q1, q2, u), interp))
    return result
