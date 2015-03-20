### Python version of Tim Hutton's artificial chemistry ###

Requires: numpy, pygame.

Usage:

    python hchem.py <rules_filename> [optional: particles_filename]
  
e.g.

    python hchem.py rules\binding_enzyme4.txt particles\binding_enzyme4b.dat

(This loads a binding enzyme. Hit 'P' to run and watch it make A22-A23 pairs!)

The terminal window is used for input, so keep it handy.

Commands:
  * 'L' to load rules (.txt files in rules folder) or particles (.dat files in particles folder).
  * 'S' to save rules or particles.
