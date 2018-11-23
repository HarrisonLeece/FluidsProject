## Fluid dynamics project, colebrook forumla recursive fxn
## Harrison Leece, Buddha

import numpy as np

#takes in roughness divided by hydraulic diameter, and reynolds number
#and returns a new reynolds number (and roughness/dia for recusion)
def coleBrook(epsilonOverDia, reynoldsNum):
    #Gives an initial low friction factor.  Reinitializes input to be shorter
    f = .00001
    e = epsilonOverDia
    r = reynoldsNum
    
    while ((1-abs(-2 * (np.sqrt(f)) * np.log10(e/3.7 + 2.51/(r*np.sqrt(f))))) > .0001):
        f = f + .00001
    return f

#for my own curiosity...
def haaland(epsilonOverDia, reynoldsNum):
    #Gives an initial low friction factor kicker.  Reinitializes input to be shorter
    f = .0001
    e = epsilonOverDia
    r = reynoldsNum
    
    while ((1-abs(-1.8 * np.sqrt(f)* np.log10((e/3.7)**1.11 + 2.51/(r*f)))) > .001):
        f = f + .0001
    return f
    
