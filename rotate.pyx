from vector import *
import math

def quaternion_multiply(a, b):
    wa, xa, ya, za = a
    wb, xb, yb, zb = b
    w = wa * wb - xa * xb - ya * yb - za * zb
    x = wa * xb + xa * wb + ya * zb - za * yb
    y = wa * yb + ya * wb + za * xb - xa * zb
    z = wa * zb + za * wb + xa * yb - ya * xb
    return [w, x, y, z]

def quaternion_conjugate(a):
    w, x, y, z = a
    return [w, -x, -y, -z]

def quaternion_from(axis, angle):
    x, y, z = unit(axis)
    theta = angle / 2.0
    w = math.cos(theta)
    sinth = math.sin(theta)
    x *= sinth
    y *= sinth
    z *= sinth
    return [w, x, y, z]

def rotate_v(vector, axis, angle):
    qv = [0.0, vector[0], vector[1], vector[2]]
    qr = quaternion_from(axis, angle)
    qrc = quaternion_conjugate(qr)
    return quaternion_multiply(quaternion_multiply(qr,qv),qrc)[1:]
