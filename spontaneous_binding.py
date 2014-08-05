# -*- coding: utf-8 -*-
# Copyright (C) 2014 nineties
# $Id: replicator.py 2014-08-05 15:50:31 nineties $

from hchem_v3 import *

sim = HChemSimulator(
    n = 500,    # number of particles
    types = ["a", "b"],
    init = [
        ("a0", 0.5),
        ("b0", 0.5),
        ],
    rules = [
        ("a0 a1 -> a1 a1", 0.99),
        ("a1 a1 -> a1 a0", 0.01),
        ("a1 a2 -> a2 a2", 0.99),
        ("a2 a2 -> a2 a1", 0.01),
        ("b0 b1 -> b1 b1", 0.99),
        ("b1 b1 -> b1 b0", 0.01),
        ("b1 b2 -> b2 b2", 0.99),
        ("b2 b2 -> b2 b1", 0.01)
        ]
    )

sim.typ[0] = 'a'
sim.state[0] = 1

sim.typ[1] = "b"
sim.state[1] = 1

sim.typ[2] = "a"
sim.state[2] = 2

sim.typ[3] = "b"
sim.state[3] = 2




# Can change positions in this fashion
# The system is 1300 x 650 (default values)
# sim.pos[3, :] = [k, k]
# Analogously: sim.vel,  sim.typ, sim.state (future development: sim.mass)
# sim.bonds[i] is a list with the particles to which each particle is bonded: sim.bonds[0] = [1], sim.bonds[1] = [0]

HChemViewer(sim).loop()
