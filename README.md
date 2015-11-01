GTOC 8 Solution Generator
=========================
Team 31 - Chris Andre, chrisandre01@gmail.com
---------------------------------------------
Feasible solution generator for GTOC 8.

The solution method relies entirely on the single chemical impulse allotted to each spacecraft. The trajectories are propagated forward until the end of the mission, taking observations as the system state moves into alignment with a radio source.

This code is unlikely to run out-of-the-box, since I made an attempt to clean it up earlier, but never got around to testing it! The logic should be preserved fairly well. The main script is run_mission.pyx.

The driving algorithm is fairly simple, and should be readily visible from the main script. A lot of the code is dedicated to last-minute runtime optimizations (hence the Cython - it started out in Python!). Even now, a single mission propagation will take a minute or two (I am going of hazy memory). You will see mention of octtrees/quadtrees; I used these to speed up the radio source angular search, borrowing from a similar tactic used in n-body simulations. I also put some effort into getting an optimal timestep, which you will see in the vari_geom package, but that is simply tuning the timestep so that the formation normal path has enough granularity so as to not skip over potential observations.

Requires `mpi4py`, `cython`, and `pykep`.
