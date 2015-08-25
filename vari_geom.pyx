cimport vari_geom

cpdef double altitude(object a, object b, object c):
    return vari_geom._altitude(a,b,c)

cpdef object tri_norm(object r1, object r2, object r3):
    return vari_geom._tri_norm(r1,r2,r3)

cpdef double get_dt(object rv1, object rv2, object rv3, double max_norm_psi, double max_body_psi):
    return vari_geom._get_dt(rv1,rv2,rv3,max_norm_psi,max_body_psi)
