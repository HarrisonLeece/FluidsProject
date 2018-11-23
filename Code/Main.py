## Fluid dynamics project, Main
## Harrison Leece (December), Buddha (December)

from fricFracFxn import coleBrook
from velocityFxn import velocityFxn
from reyNumFxn import reyNum
import numpy as np

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

#first length of pipe vars
length = 205
elevation = 25
pHead = 30.58
D = D
pressure = 150000
pressureByGamma = pressure/(g*rho)
minorLoss = 11.3

#Kicker lines
frFactor = coleBrook(ed1, reGuess)
print(frFactor)
v1 = velocityFxn(elevation, pressureByGamma, pHead, frFactor, length, D, minorLoss,g)

#recursivly calls reyNum until change in Reynolds number per fxn call is small
#returns reynolds, frictionFac, v
rey, ff, vf = reyNum(D, reGuess, v1, rho, mu, elevation, pressureByGamma, pHead, frFactor, length, minorLoss, g, ed1, 0)
vf = velocityFxn(elevation, pressureByGamma, pHead, ff, length, D, minorLoss,g)


print('\n' + 'Final velocity: ' +str(vf))
