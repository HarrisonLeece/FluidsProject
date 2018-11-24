## Fluid dynamics project, Main
## Harrison Leece (December), Buddha (December)

from fricFracFxn import coleBrook
from velocityFxn import velocityFxn
from reyNumFxn import reyNum
import numpy as np
import matplotlib as plt

#Gravity m/s^2
g = 9.81
#density in kg/m^3
rho = 1000
#viscosity of water
mu = 8.4*10**(-4)

#Source volumetric flow rate in m^3/s
Qin = .02120
#Source static pressure in kPa
Pin = 15
#Roughness and hydraulic diameter
e = .00026
D = .1
#Roughness over diameter
ed1 = e/D
#initial guess at reynolds number (in this case 10^5 for water
reGuess = 1*10**5

################################################################################
#Pipe vars
#first length of main pipe vars
length1 = 44
elevation1 = 4
pHead = 30.58
d1 = 4
pressure1 = 15000
pressureByGamma1 = pressure/(g*rho)
#One elbow fitting between node 1 and 2
minorLoss1 = .3

#Second length of main pipe vars
length1 = 4
elevation1 = 4
pHead = 0
d1 = 4
pressure1 = 15000
pressureByGamma1 = pressure/(g*rho)
minorLoss1 = .3

#Second length of main pipe vars
length2 = 4
elevation1 = 4
pHead = 0
d1 = 4
pressure1 = 15000
pressureByGamma1 = pressure/(g*rho)
minorLoss1 = .3

#Third length of main pipe vars
length3 = 4
elevation1 = 4
pHead = 0
d1 = 4
pressure1 = 15000
pressureByGamma1 = pressure/(g*rho)
minorLoss1 = .3

#firstOutlet
length4 = 0
elevation1 = 0
pHead = 0
d1 = 2
pressure1 = 15000
pressureByGamma1 = pressure/(g*rho)
minorLoss1 = 30

#secondOutlet
length5 = 5
elevation1 = 0
pHead = 0
d1 = 2
pressure1 = 15000
pressureByGamma1 = pressure/(g*rho)
minorLoss1 = .3

#thirdOutlet
length5 = 4
elevation1 = 0
pHead = 0
d1 = 2
pressure1 = 15000
pressureByGamma1 = pressure/(g*rho)
minorLoss1 = .3

#############################################################################

vin = Qin/(3.1416 * (d1/2)**2)

#Kicker lines to obtain initial velocity
frFactor = coleBrook(ed1, reGuess)
print(frFactor)
v1 = velocityFxn(elevation, pressureByGamma, pHead, frFactor, length, D, minorLoss,g)

#recursivly calls reyNum until change in Reynolds number per fxn call is small
#returns reynolds, frictionFac, v
rey, ff, vf = reyNum(D, reGuess, v1, rho, mu, elevation, pressureByGamma, pHead, frFactor, length, minorLoss, g, ed1, 0)
vf = velocityFxn(elevation, pressureByGamma, pHead, ff, length, D, minorLoss,g)


print('\n' + 'Final velocity: ' +str(vf))
