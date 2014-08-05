# Copyright (C) 2014 nineties
# $Id: replicator.py 2014-08-05 11:19:10 nineties $

from hchem_v3 import *

sim = HChemSimulator(
    n = 100,    # number of particles
    types = ["a", "b", "c", "d", "e", "f"],
    wildcards = ["x", "y"],
    state_max = 37, 
    init = [
        ("a0", 3.0/8),
        ("b0", 1.0/8),
        ("c0", 1.0/8),
        ("d0", 1.0/8),
        ("e0", 1.0/8),
        ("f0", 1.0/8)
        ],
    rules = [
        "e1-a37  -> e5-a10",
        "a10 e6  -> a37-e3",
        "e6-e3   -> e2-e3",
        "x2-y1   -> x7-y4",
        "x4 y3   -> x5-y7",
        "x5 x0   -> x6-x6",
        "x6 y7   -> x3-y4",
        "x6-y4   -> x1 x2",
        "x7-y1   -> x2-y2",
        "f2-a37  -> f9-a11",
        "a11 f3  -> a11-f9",
        "x2-y8   -> x9-y1",
        "x9-y9   -> x8 y8",
        "a11-a36 -> a11-a12",
        "f1 a12  -> f13-a37",
        "x13-y1  -> x14-y15",
        "a11 x15 -> a11-x16",
        "x14-y16 -> x27-y16",
        "x27-a11 -> x17 a11",
        "x17-y16 -> x17-y13",
        "x13-e8  -> x14-e15",
        "e13-a37 -> e18-a19",
        "e13-a19 -> e18-a20",
        "a20 a11 -> a21-a22",
        "e18-a22 -> e32 a23",
        "e18-a21 -> e32 a24",
        "a24 a37 -> a26 a27",
        "a27-a23 -> a37 a28",
        "a26-a36 -> a29-a30",
        "a29-a36 -> a31-a30",
        "a30 a28 -> a25-a33",
        "a31-a25 -> a32 a36",
        "a32-a30 -> a34-a36",
        "a34-a33 -> a37 a37",
        ]
    )

HChemViewer(sim).loop()
