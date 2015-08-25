from libc.math cimport sin, cos
from math import pi

def parse_radio():
    f = open('gtoc8_radiosources.txt')
    lines = f.read().split('\n')
    list = []
    for ln in lines:
        if len(ln) is not 0 and ln[0] is not '%':
            split = ln.split()
            list.append((split[0], split[1], split[2]))
    xyz = []
    radecl = []
    id = []
    cdef double rai, decli, x, y, z
    for source in list:
        rai = float(source[1]) * pi / 180.0
        decli = float(source[2]) * pi / 180.0
        radecl.append((rai, decli))
        x = cos(decli) * cos(rai)
        y = cos(decli) * sin(rai)
        z = sin(decli)
        xyz.append([x, y, z])
        id.append(int(source[0]))
        #print '{' + source[1] + ',' + source[2] + '},'

    return xyz, radecl, id