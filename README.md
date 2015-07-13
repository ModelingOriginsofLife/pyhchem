### Python version of Tim Hutton's artificial chemistry ###

Requires: numpy, pygame.

Usage:

    python hchem.py <rules_filename> [optional: particles_filename]
  
e.g.

    python hchem.py rules\programmable_enzyme.txt particles\programmable_enzyme.dat

This loads a programmable enzyme. Hit 'P' to run and watch it bond A0 and A1 atoms into A2-A3 pairs, if you're lucky.

The terminal window is used for input, so keep it handy.

Commands:
  * 'L' to load rules (.txt files in rules folder) or particles (.dat files in particles folder).
  * 'S' to save rules or particles.

![pyhchem_programmable_enzyme](https://cloud.githubusercontent.com/assets/647092/6994944/36f8596a-db23-11e4-849f-96dd1a940080.png)

Vimeo: https://vimeo.com/133396645
