## Fluid dynamics project, compute velocity
## Harrison Leece, Buddha

import numpy as np

def velocityFxn(elevation, pressureByGamma, pumpHead, frictionFac, length, dia, minorLoss,g):
    v = np.sqrt(((elevation - pressureByGamma - pumpHead)*2*g)/(1-(frictionFac*length/dia) - minorLoss))
    return v
