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
diameter = D
pressure = 150000
pressureByGamma = pressure/(g*rho)
minorLoss = 11.3

frFactor = coleBrook(ed1, reGuess)
print(frFactor)


v1 = velocityFxn(elevation, pressureByGamma, pHead, frFactor, length, D, minorLoss,g)

print(v1)
