## Fluid dynamics project, reynold computation recursive fxn
## Harrison Leece, Buddha

import numpy as np

#takes in roughness divided by hydraulic diameter, and reynolds number
#and returns a new reynolds number (and roughness/dia for recusion)
def reyNum(epsilon_dia, reynolds, velocity, density, viscosity):
    oldRey = reynolds
    
    
    if ( abs(reynolds-oldRey) > .001 ):
        reyNum(epsilon_dia, reynolds)
    else:
        return reynolds
    
    
