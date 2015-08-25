cimport octtree
from octtree import OctTree
from libc.math cimport sin
from vector import unit, vangle

class NearestAngular(object):

    def __init__(self, pts):
        self._pts = pts
        self._oct = OctTree([0.0, 0.0, 0.0], [3.0, 3.0, 3.0])
        i = 0
        for pt in pts:
            self._oct.place(unit(pt), i)
            i+=1

    def find(self, pt, angle):
        cdef double chord = sin(angle/2.0)*2.0
        ret = self._oct.find(unit(pt), chord)
        return ret