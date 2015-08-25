"""Quick vector operations"""

from libc.math cimport sqrt, acos

cdef inline void cross3(double *out, double u[3], double v[3]):
    out[0] = u[1]*v[2]-v[1]*u[2]
    out[1] = u[2]*v[0]-v[2]*u[0]
    out[2] = u[0]*v[1]-v[0]*u[1]

cdef inline double dot3(double u[3], double v[3]):
    return u[0]*v[0] + u[1]*v[1] + u[2]*v[2]

cdef inline double mag3(double u[3]):
    return sqrt(u[0]**2+u[1]**2+u[2]**2)

cdef inline double vangle3(double u[3], double v[3]):
   return acos(dot3(u,v)/(mag3(u)*mag3(v)))

cdef inline void sub3(double *out, double u[3], double v[3]):
    out[0] = u[0]-v[0]
    out[1] = u[1]-v[1]
    out[2] = u[2]-v[2]

cdef inline void add3(double *out, double u[3], double v[3]):
    out[0] = u[0]+v[0]
    out[1] = u[1]+v[1]
    out[2] = u[2]+v[2]

cdef inline void scale3(double *out, double u[3], double a):
    out[0] = u[0]*a
    out[1] = u[1]*a
    out[2] = u[2]*a

cdef inline void unit3(double *out, double u[3]):
   scale3(out, u, 1.0/mag3(u))

# cdef inline object transform3(x_t, x_r, radius, normal):
#    """Calculate 3D location"""
#    ur = unit3(radius)
#    ut = unit3(cross3(normal, radius))
#
#    return add3(scale3(ur, x_r), scale3(ut, x_t))
