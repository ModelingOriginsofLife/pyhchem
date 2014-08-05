# Copyright (C) 2014 nineties
# $Id: replicator.py 2014-08-05 15:16:33 nineties $

from hchem_v3 import *

sim = HChemSimulator(
    n = 500,    # number of particles
    types = ["a", "b", "c"],
    init = [
        ("a0", 0.6),
        ("b0", 0.4),
        ],
    rules = [
        "a0 a0 -> a1-a1",
        "a1 a0 -> a2-a1",
        ("a1 a1 -> a2-a2", 0.99),
        ("a2-a2 -> a1 a1", 0.001),
        ]
    )

# sim.pos
# sim.vel
# sim.typ
# sim.state
# sim.mass
# sim.bonds

#for k in range(10):
#    sim.pos[k,:] = [k, k]

# 0, 1
#sim.bonds[0] = [1]
#sim.bonds[1] = [0, 2]
#sim.bonds[2] = [1]

HChemViewer(sim).loop()
