'''
Utility for converting between Cython double arrays and Python array types
'''

cdef inline void pack_double_array(double *u, object double_list):
    for i in range(len(double_list)):
        u[i] = double_list[i]

cdef inline object unpack_double_array(double *u, int n):
    return [u[i] for i in range(n)]