"""Quick vector operations"""

import math

def cross(u, v):
    ux, uy, uz = u
    vx, vy, vz = v
    return [uy*vz-vy*uz, uz*vx-vz*ux, ux*vy-vx*uy]

def dot(u, v):
    return sum([ui*vi for ui, vi in zip(u, v)])

def mag(u):
    return math.sqrt(sum([xi*xi for xi in u]))

def vangle(u, v):
    return math.acos(dot(u,v)/(mag(u)*mag(v)))

def sub(u, v):
    return [ui - vi for ui, vi in zip(u, v)]

def add(u, v):
    return [ui + vi for ui, vi in zip(u, v)]

def scale(u, a):
    return [ui * a for ui in u]

def unit(u):
    return scale(u, 1.0/mag(u))


def transform(x_t, x_r, radius, normal):
    """Calculate 3D location"""
    ur = unit(radius)
    ut = unit(cross(normal, radius))

    return add(scale(ur, x_r), scale(ut, x_t))
