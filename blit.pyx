'''
Utility for 'blitting' spacecraft states to GTOC8's required format
'''

class BlitSPC(object):

    def __init__(self, fname, i):
        self.fname = fname
        self.f = open(fname, 'a+')
        self.f.write('# Spacecraft %i\n' % i)
        self.f.write('# MJD, x, y, z, vx, vy, vz, mass, tx, ty, tz\n')

    def finish(self):
        self.f.close()

    def blit_state(self, mjd, r, v, t, m):
        self.f.write(' %.11f ' % mjd)  # t
        self.f.write(' %.15e ' % (r[0] / 1000))  # r
        self.f.write(' %.15e ' % (r[1] / 1000))
        self.f.write(' %.15e ' % (r[2] / 1000))
        self.f.write(' %.15e ' % (v[0] / 1000))  # v
        self.f.write(' %.15e ' % (v[1] / 1000))
        self.f.write(' %.15e ' % (v[2] / 1000))
        self.f.write(' %.15e ' % m) # m
        self.f.write(' %.15e ' % 0.0)  # T
        self.f.write(' %.15e ' % 0.0)
        self.f.write(' %.15e ' % 0.0)
        self.f.write('\n')

    def blit_impulse(self, mjd, r, vi, vf, t, mi, mf):
        self.blit_state(mjd, r, vi, t, mi)
        self.blit_state(mjd, r, vf, t, mf)

class BlitOBS(object):

    def __init__(self, fname):
        self.fname = fname
        self.f = open(fname, 'a+')
        self.f.write('# Observations\n')
        self.f.write('# MJD, ID, P, h, s1x, s1y, s1z, s2x, s2y, s2z, s3x, s3y, s3z, nx, ny, nz\n')

    def finish(self):
        self.f.close()

    def blit_observation(self, mjd, id, P, h, sr1, sr2, sr3, norm):
        self.f.write(' %.11f ' % mjd)  # t
        self.f.write(' %3i ' % id)  # ID
        self.f.write(' %i ' % P)  # P
        self.f.write(' %.15e ' % (h / 1000))  # h
        self.f.write(' %.15e ' % (sr1[0] / 1000))  # s1
        self.f.write(' %.15e ' % (sr1[1] / 1000))
        self.f.write(' %.15e ' % (sr1[2] / 1000))
        self.f.write(' %.15e ' % (sr2[0] / 1000))  # s1
        self.f.write(' %.15e ' % (sr2[1] / 1000))
        self.f.write(' %.15e ' % (sr2[2] / 1000))
        self.f.write(' %.15e ' % (sr3[0] / 1000))  # s1
        self.f.write(' %.15e ' % (sr3[1] / 1000))
        self.f.write(' %.15e ' % (sr3[2] / 1000))
        self.f.write(' %.15e ' % (norm[0]))  # n
        self.f.write(' %.15e ' % (norm[1]))
        self.f.write(' %.15e ' % (norm[2]))
        self.f.write('\n')