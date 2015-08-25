cdef class OctTree:
    cdef double _center[3]
    cdef double _lwh[3]
    cdef double _pt[3]
    cdef object _metadata
    cdef int _leaf_node
    cdef int _has_pt

    # Cython won't allow an array:
    cdef OctTree _sub000
    cdef OctTree _sub001
    cdef OctTree _sub010
    cdef OctTree _sub011
    cdef OctTree _sub100
    cdef OctTree _sub101
    cdef OctTree _sub110
    cdef OctTree _sub111

    cdef inline int _bin(self, double pt[3])
    cdef inline OctTree _unbin(self, int flag)
    cdef inline void _ensure_sub(self, int bin)
    cdef inline void _place(self, double pt[3], object metadata)
    cdef inline void _find(self, double pt[3], double r, object list)