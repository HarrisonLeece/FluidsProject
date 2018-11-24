## Fluid dynamics project, reynold computation recursive fxn
## Harrison Leece, Buddha

import numpy as np
from velocityFxn import velocityFxn
from fricFracFxn import coleBrook

#takes in roughness divided by hydraulic diameter, and reynolds number
#and returns a new reynolds number (and roughness/dia for recusion)
def reyNum(dia, rey, velocity, density, viscosity, elevation, pressureByGamma, pumpHead, frictionFac, length, minorLoss, g, ed, n):
    oldRey = rey
    reynolds = (density * velocity * dia )/ viscosity
    frictionFac = coleBrook(ed,reynolds)
    v = velocityFxn(elevation, pressureByGamma, pumpHead, frictionFac, length, dia, minorLoss,g)
    
    if ( abs(reynolds - oldRey) > 10 ):
        print('Iteration number: ' + str(n))
        n = n+1
        return reyNum(dia, reynolds, v, density, viscosity, elevation, pressureByGamma, pumpHead, frictionFac, length, minorLoss, g, ed, n)
    else:
        print('Else check: ReynoldsDiff = ' + str(reynolds - oldRey))
        print('Velocity ' + str(v))
        print('Reynolds Number ' + str(reynolds))
        print('Friction Factor ' + str(frictionFac))
        return reynolds, frictionFac, v
    
    
