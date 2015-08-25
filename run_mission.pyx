'''
Main script to generate feasible trajectories.
'''
########################################################################################################################
# STATIC INITIALIZATION
########################################################################################################################
import random
import PyKEP
import math
import orbit
import scene
import vector
import rotate
import vari_geom
import shelve
import blit
import mpi4py.MPI
# MPI Process ID
RANK = mpi4py.MPI.COMM_WORLD.rank
# Observation function; returns details (id, h, norm) for a radio source if found, defaults (0,0,0)
def observe(r1, v1, r2, v2, r3, v3):
    # Both normal vectors to formation
    norm_a = vari_geom.tri_norm(r1, r2, r3)
    norm_b = vector.scale(norm_a, -1.0)
    # Get sources within tolerance
    aa = scene.constellation.find(norm_a, scene.ObsTol)
    ab = scene.constellation.find(norm_b, scene.ObsTol)
    if len(aa) >= 1 or len(ab) >= 1:  # found any source (n > 1 is impossible)
        # Load found source into 'first' and corresponding normal
        if len(aa) >= 1:
            first = aa[0]
            fnorm = norm_a
        else:
            first = ab[0]
            fnorm = norm_b
        # Get source zero-based ID
        ordinal = first[3]
        # Get source ID as per radio source file
        id = ordinal + 1
        # Calculate minimum altitude for formation
        h = vari_geom.altitude(r1, r2, r3)
        return id, h, fnorm
    else:
        return 0, 0, 0
########################################################################################################################
# STEP
########################################################################################################################
def step():
    print "Node %i iteration!" % RANK
    # Initial spacecraft states
    r0, v0 = scene.s0.eph(scene.T0)
    # Randomize initial impulse delta-v
    dv = [random.uniform(2500.0, 3000.0) for i in range(3)]
    # Determine angle to apply impulse that keeps vf strictly greater than the pre-impulse vi
    maxang = [math.pi - math.acos(0.5 * dvi / vector.mag(v0)) for dvi in dv]
    ang = [random.uniform(-mxai, mxai) for mxai in maxang]
    # Randomize spacecraft impulse timing to be under a single orbit period for the initial state
    dt = [random.uniform(0.0, orbit.period(scene.X0, scene.MUE) * PyKEP.SEC2DAY) for i in range(3)]  # days
    #dt[0] = 0.0 # We get faster and slightly better results
    # Calculate final masses for impulses
    mf = [scene.get_mf(4000.0, dvi) for dvi in dv]
    # Generate planet objects that represent the spacecraft in their final states
    sf = []
    for i in range(3):
        rotang = 2 * math.pi * 86400 * dt[i] / orbit.period(scene.X0, scene.MUE)
        vp = rotate.rotate_v([0, dv[i] * math.cos(ang[i]), dv[i] * math.sin(ang[i])], [0, 0, 1], rotang)
        ti = PyKEP.epoch(scene.T0.mjd + dt[i], "mjd")
        r, v = scene.s0.eph(ti)
        vn = vector.add(v, vp)
        sf.append(PyKEP.planet.keplerian(ti, r, vn, scene.MUE, 0, 1, 2))
    # Minimum time for an observation
    startt = max(dt) + scene.T0
    # Maximum time for an observation (mission end)
    endt = min(dt) + scene.T0 + 3 * 365.25
    # Source dartboard to count individual observations
    sourcehitcount = [0] * 421
    sourceh = [[0] * 421, [0] * 421, [0] * 421]
    obs = []
    J = 0.0
    t = startt
    last_t = t
    while t < endt:
        # Current state of spacecraft formation
        ep = PyKEP.epoch(t, "mjd")
        r1, v1 = sf[0].eph(ep)
        r2, v2 = sf[1].eph(ep)
        r3, v3 = sf[2].eph(ep)
        id, h, fnorm = observe(r1, v1, r2, v2, r3, v3)
        if id != 0 and h >= 10000.0e3:
            # Count of times the formation has struck a source
            n = sourcehitcount[id]
            if n < 3:
                # Add current altitude to list of previous altitudes
                sourceh[n][id] = h
                P = scene.get_P(n + 1, sourceh[0][id], sourceh[1][id], sourceh[2][id])
                if P > 0:
                    dj = scene.perf_index(P, h, scene.radecl[id - 1][1])
                    J += dj
                    print 'Node %i observation! id: %3i, t: %8.3f, del-t: %8.3f, P: %i, hit: %i, J: %5.1e km, del-J: %5.1e km, status: %4.1f percent' % (
                    RANK, id, t - startt, t - last_t, P, n + 1, J / 1000, dj / 1000,
                    100 * (t - startt) / (endt - startt))
                    obs.append([t, id, P, h, r1, r2, r3, fnorm])
                    last_t = t
                    t += 15.0
                    sourcehitcount[id] += 1
        # Calculate highest timestep that preserves fidelity
        ddt = vari_geom.get_dt((r1, v1), (r2, v2), (r3, v3), 0.2 * scene.ObsTol, 10.0 * PyKEP.DEG2RAD) / 86400.0
        t += ddt
        if t - last_t > 60.0:  # quit
            print 'Node %i quitting early, mediocre initial conditions...' % RANK
            break
    ####################################################################################################################
    # STORAGE AND STATISTICS
    ####################################################################################################################
    print '\nNode %i finished; J := %.3e km, observations := %i\n' % (mpi4py.MPI.COMM_WORLD.rank, J / 1000.0, len(obs))

    folder = './solutions/'

    fname = 'J := %0.3e km' % (J / 1000.0)
    x = shelve.open(folder + fname + '-shelf', writeback=True)
    x["obs"] = obs
    x["dt"] = dt
    x["mf"] = mf
    x["sf"] = sf
    x["J"] = J
    x.close()

    f = open(folder + fname + '-blit', 'a+')
    f.write('# J := %.16e km\n' % (J / 1000.0))
    f.close()

    out_obs = blit.BlitOBS(folder + fname + '-blit')
    for ob in obs:
        out_obs.blit_observation(*ob)
    out_obs.finish()

    for i, sp in enumerate(sf):
        out_spc = blit.BlitSPC(folder + fname + '-blit', i)
        if dt[i] != min(dt):
            ep = PyKEP.epoch(scene.T0.mjd + min(dt), "mjd")
            out_spc.blit_state(ep.mjd, scene.s0.eph(ep)[0], scene.s0.eph(ep)[1], [0, 0, 0], 4000.0)
        ep = PyKEP.epoch(scene.T0.mjd + dt[i], "mjd")
        ri, vi = scene.s0.eph(ep)
        rf, vf = sp.eph(ep)
        out_spc.blit_impulse(ep.mjd, ri, vi, vf, [0, 0, 0], 4000.0, mf[i])
        out_spc.finish()
########################################################################################################################
# RUN
########################################################################################################################
while True:
    step()
