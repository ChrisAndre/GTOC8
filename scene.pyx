'''
Contains the specifics of the GTOC8 problem description
'''

import PyKEP as pkp
import math
import parse_radio
import nearest_angular

RE = 6378.14e3
RM = 1737.5e3
G = 9.80665
YR = 365.25
MUE = 398600.4329 * 1000**3
MUM = 4902.8006 * 1000**3

T0 = pkp.epoch(58849.0, "mjd")
X0 = RE + 400e3
A = X0
E = 0.0
I = 0.0
IspC = 450.0
MaxImpulseC = 3000.0
IspLT = 5000.0
MaxThrustLT = 0.1
T0Max = pkp.epoch(58880.0, "mjd")
TOFMax = 3.0 * YR
MI = 4000.0
MF = 1890.0
MinR = 6578.14e6
MaxR = 1e9
FlyByMinR = RM + 50.0e3
FlyByMinVInf = 0.25e3
ObsDTMin = 15.0
AltMin = 10000.0e3
ObsTol = 0.1 * pkp.DEG2RAD

s0 = pkp.planet.keplerian(T0, [X0, 0.0, 0.0, 0.0, 1.0e-5, 0.0], MUE, 0, 1, 2)
moon = pkp.planet.keplerian(T0, [383500.0e3,0.04986,5.2586,98.0954*pkp.DEG2RAD,69.3903*pkp.DEG2RAD,164.35025*pkp.DEG2RAD],MUE,MUM,RM,FlyByMinR,'moon')

points, radecl, ids = parse_radio.parse_radio()
constellation = nearest_angular.NearestAngular(points)

def perf_index(P, h, decl):
    if h < 10000.0e3:
        return 0
    return P * h * (0.2 + math.cos(decl)**2)

def get_P(n, h1, h2, h3):
    if n == 1:
        return 1
    if n == 2:
        hmax = max(h1, h2)
        hmin = min(h1, h2)
        if hmax / hmin >= 3:
            return 3
        else:
            return 1
    if n == 3:
        hmax = max(h1, h2, h3)
        hmin = min(h1, h2, h3)
        hmid = h1 + h2 + h3 - hmax - hmin
        if hmax / hmid >= 3 and hmid / hmin > 3:
            return 6
        else:
            if hmax / hmid >= 3 and max(h1, h2) / min(h1, h2) < 3:
                return 3
            else:
                return 1
    return 0

def get_mf(mi, dv):
    return mi * math.exp(-dv/(9.80665*450.0))