# Copyright (C) 2014 nineties
# $Id: replicator.py 2014-08-05 12:24:50 nineties $

from hchem_v3 import *

sim = HChemSimulator(
    n = 500,    # number of particles
    types = ["a", "b", "c"],
    state_max = 1, 
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

HChemViewer(sim).loop()
