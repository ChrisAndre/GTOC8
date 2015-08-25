from vector import *
import math

def circ_vis_viva(r, mu):
    return math.sqrt(mu/mag(r))

def circ_v(normal, r, mu):
    v = math.sqrt(mu /mag(r))
    vdir = unit(cross(normal, r))
    return scale(vdir, v)

def circ_tof(r, mu, psi):
    r = mag(r)
    v = math.sqrt(mu / r)
    return psi * r / v

def circ_psi(r, mu, tof):
    r = mag(r)
    v = math.sqrt(mu / r)
    return v * tof / r

def period(a, mu):
    return math.pi * 2 * math.sqrt(a**3 / mu)
