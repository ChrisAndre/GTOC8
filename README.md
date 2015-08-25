GTOC 8 Solution Generator
=========================
Team 31 - Chris Andre, chrisandre01@gmail.com
---------------------------------------------
Feasible solution generator for the GTOC 8 competition.

The solution method relies entirely on the single chemical impulse allotted to each spacecraft. The trajectories are propagated forward until the end of the mission, taking observations as the system state moves into alignment with a radio source.

Requires `mpi4py`, `cython`, and `pykep`.