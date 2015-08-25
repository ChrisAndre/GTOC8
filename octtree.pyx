from octtree cimport OctTree
from vec3 cimport mag3

cdef class OctTree:
    def __cinit__(self, center, lwh):
        for i in range(3):
            self._center[i] = center[i]
            self._lwh[i] = lwh[i]
        self._leaf_node = 1
        self._has_pt = 0

    cdef inline int _bin(self, double pt[3]):
        cdef int flag = 0b000
        if pt[0] >= self._center[0]:
            flag |= 0b001
        if pt[1] >= self._center[1]:
            flag |= 0b010
        if pt[2] >= self._center[2]:
            flag |= 0b100
        return flag

    cdef inline OctTree _unbin(self, int flag):
        if flag == 0b000:
            return self._sub000
        elif flag == 0b001:
            return self._sub001
        elif flag == 0b010:
            return self._sub010
        elif flag == 0b011:
            return self._sub011
        elif flag == 0b100:
            return self._sub100
        elif flag == 0b101:
            return self._sub101
        elif flag == 0b110:
            return self._sub110
        elif flag == 0b111:
            return self._sub111
        else:
            return None

    cdef inline void _ensure_sub(self, int bin):
        lwh2 = [self._lwh[i]/2.0 for i in range(3)]
        cdef double cx, cy, cz
        cx = self._center[0] - 0.5 * lwh2[0] * (-1.0)**((bin & 0b001) > 0)
        cy = self._center[1] - 0.5 * lwh2[1] * (-1.0)**((bin & 0b010) > 0)
        cz = self._center[2] - 0.5 * lwh2[2] * (-1.0)**((bin & 0b100) > 0)
        pt = [cx, cy, cz]
        if bin == 0b000 and self._sub000 == None:
            self._sub000 = OctTree(pt, lwh2)
        elif bin == 0b001 and self._sub001 == None:
            self._sub001 = OctTree(pt, lwh2)
        elif bin == 0b010 and self._sub010 == None:
            self._sub010 = OctTree(pt, lwh2)
        elif bin == 0b011 and self._sub011 == None:
            self._sub011 = OctTree(pt, lwh2)
        elif bin == 0b100 and self._sub100 == None:
            self._sub100 = OctTree(pt, lwh2)
        elif bin == 0b101 and self._sub101 == None:
            self._sub101 = OctTree(pt, lwh2)
        elif bin == 0b110 and self._sub110 == None:
            self._sub110 = OctTree(pt, lwh2)
        elif bin == 0b111 and self._sub111 == None:
            self._sub111 = OctTree(pt, lwh2)

    cdef inline void _place(self, double pt[3], object metadata):
        cdef int ptflag
        cdef OctTree sub
        if self._leaf_node == 1:
            if self._has_pt == 1:
                self._leaf_node = 0
                self._has_pt = 0

                ptflag = self._bin(self._pt)
                self._ensure_sub(ptflag)
                sub = self._unbin(ptflag)
                sub._place(self._pt, self._metadata)

                ptflag = self._bin(pt)
                self._ensure_sub(ptflag)
                sub = self._unbin(ptflag)
                sub._place(pt, metadata)
            else:
                self._has_pt = 1
                for i in range(3):
                    self._pt[i] = pt[i]
                self._metadata = metadata
        else:
            ptflag = self._bin(pt)
            self._ensure_sub(ptflag)
            sub = self._unbin(ptflag)
            sub._place(pt, metadata)

    cdef inline void _find(self, double pt[3], double r, object list):
        cdef double d[3]
        cdef double dist
        cdef double x,y,z
        cdef OctTree unbinned
        for i in range(3):
            d[i] = pt[i] - self._center[i]
        dist = mag3(d)
        if dist > 1.7320509 * 0.5 * max(self._lwh[0], self._lwh[1], self._lwh[2]) + r:
            return
        else:
            if self._has_pt:
                for i in range(3):
                    d[i] = self._pt[i] - pt[i]
                if mag3(d) <= r:
                    x = self._pt[0]
                    y = self._pt[1]
                    z = self._pt[2]
                    list.append([x,y,z,self._metadata])
            for i in range(8):
                unbinned = self._unbin(i)
                if unbinned != None:
                    unbinned._find(pt, r, list)
            return

    def place(self, pt, meta):
        cdef double _pt[3]
        for i in range(3):
            _pt[i] = pt[i]
        self._place(_pt, meta)

    def find(self, pt, double r):
        cdef double _pt[3]
        for i in range(3):
            _pt[i] = pt[i]
        ret = []
        self._find(_pt, r, ret)
        return ret