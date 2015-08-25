from vec3 cimport sub3, vangle3, mag3, cross3, scale3, add3, unit3
from libc.math cimport sin
cimport marshalling

cdef inline double _altitude(object _a, object _b, object _c):
    cdef double a[3], b[3], c[3]
    cdef double l1[3], l2[3], l3[3]
    cdef double aa, ab, ac
    cdef double l1m, l2m, l3m
    cdef double ha, hb, hc
    marshalling.pack_double_array(a, _a)
    marshalling.pack_double_array(b, _b)
    marshalling.pack_double_array(c, _c)
    sub3(l1, a, b)
    sub3(l2, b, c)
    sub3(l3, c, a)
    aa = vangle3(l1, l2)
    ab = vangle3(l2, l3)
    ac = vangle3(l3, l1)
    l1m = mag3(l1)
    l2m = mag3(l2)
    l3m = mag3(l3)
    ha = sin(aa) * l1m
    hb = sin(ab) * l2m
    hc = sin(ac) * l3m
    return min(ha, hb, hc)

cdef inline void _dcross3(double *out, double u[3], double du[3], double k[3], double dk[3]):
    cdef double a[3], b[3]
    cross3(a, du, k)
    cross3(b, u, dk)
    add3(out, a, b)

cdef inline void _tri_norm_omega(double *out, double a[3], double da[3], double b[3], double db[3], double c[3], double dc[3]):
    cdef double u[3], du[3], v[3], dv[3]
    cdef double norm[3], dnorm[3]
    sub3(u, b, a)
    sub3(du, db, da)
    sub3(v, c, a)
    sub3(dv, dc, da)
    cross3(norm, u, v)
    _dcross3(dnorm, u, du, v, dv)
    _omega3(out, norm, dnorm)

cdef inline void _omega3(double *out, double r[3], double dr[3]):
    cdef double m2 = mag3(r)**2.0
    cdef double cr[3]
    cross3(cr, r, dr)
    scale3(out, cr, 1.0/m2)

cdef inline object _tri_norm(object r1, object r2, object r3):
    cdef double _r1[3], _r2[3], _r3[3]
    cdef double a[3], b[3]
    cdef double norm[3]
    marshalling.pack_double_array(_r1, r1)
    marshalling.pack_double_array(_r2, r2)
    marshalling.pack_double_array(_r3, r3)
    sub3(a, _r2, _r1)
    sub3(b, _r3, _r1)
    unit3(a, a)
    unit3(b, b)
    cross3(norm, a, b)
    unit3(norm, norm)
    return marshalling.unpack_double_array(norm, 3)

cdef inline double _get_dt(object rv1, object rv2, object rv3, double max_norm_psi, double max_body_psi):
    cdef double _r1[3], _v1[3]
    cdef double _r2[3], _v2[3]
    cdef double _r3[3], _v3[3]
    cdef double omt[3], om1[3], om2[3], om3[3]
    cdef double mt, m1, m2, m3
    cdef double max_body_dt, max_norm_dt
    cdef double ret
    r1, v1 = rv1
    r2, v2 = rv2
    r3, v3 = rv3
    marshalling.pack_double_array(_r1, r1)
    marshalling.pack_double_array(_r2, r2)
    marshalling.pack_double_array(_r3, r3)
    marshalling.pack_double_array(_v1, v1)
    marshalling.pack_double_array(_v2, v2)
    marshalling.pack_double_array(_v3, v3)
    _omega3(om1, _r1, _v1)
    _omega3(om2, _r2, _v2)
    _omega3(om3, _r3, _v3)
    _tri_norm_omega(omt, _r1, _v1, _r2, _v2, _r3, _v3)
    m1 = mag3(om1)
    m2 = mag3(om2)
    m3 = mag3(om3)
    mt = mag3(omt)
    max_body_dt = max_body_psi / max(m1, m2, m3)
    max_norm_dt = max_norm_psi / mt
    ret = min(max_body_dt, max_norm_dt)
    return ret