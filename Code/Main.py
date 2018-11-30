## Fluid dynamics project, Main
## Harrison Leece (December), Buddha (December)

from fricFracFxn import coleBrook
from velocityFxn import velocityFxn
from reyNumFxn import reyNum
from scipy.optimize import fsolve
import numpy as np
import matplotlib as plt

def equations(p, ffs):
	v1 = p[0]
	v2 = p[1]
	v3 = p[2]
	Hp = p[3]
	ff1 = ffs[0]
	ff2 = ffs[1]
	ff3 = ffs[2]

##	print('FF1: ' +str(ff1))
##	print('FF2: ' +str(ff2))
##	print('FF3: ' +str(ff3))
##	print('Pump head: ' +str(Hp))
##	print('v1: ' +str(v1))
##	print('v2: ' + str(v2))
##	print('v3: ' + str(v3))
	#Define constants
	#Gravity m/s^2
	g = 9.81
	#density in kg/m^3
	rho = 1000
	#viscosity of water
	mu = 8.4*10**(-4)
	#initial volumetric flow rate
	Qin = .02120
	#initial Pressure
	P0 = 15000
	p0g = P0/(g*rho)
	#length
	l1 = 36
	l2 =35
	l3 = 34
	#Diameter 
	D = 2/39.37
	A1 = 3.1416/4 * D**2
	A0 = 3.1416/4 * (4/39.37)**2
	#initial flow velocity
	v0 = Qin/A0

##	v1 = np.sqrt((2*g*(p0g+(v0**2/(2*g))+Hp-4))/(1+ff1*l1/D+.3))
##	v2 = np.sqrt((2*g*(p0g+(v0**2/(2*g))+Hp-8))/(1+ff2*l2/D+.3))
##	v3 = np.sqrt((2*g*(p0g+(v0**2/(2*g))+Hp-4))/(1+ff3*l3/D+.3))
##      v1+v2+v3 = v0*A0/A1
	
	return ( np.sqrt((2*g*(p0g+(v0**2/(2*g))+Hp-4))/(1+ff1*l1/D+.3)) -v1,np.sqrt((2*g*(p0g+(v0**2/(2*g))+Hp-8))/(1+ff2*l2/D+.3)) -v2,np.sqrt((2*g*(p0g+(v0**2/(2*g))+Hp-4))/(1+ff3*l3/D+.3)) -v3,v1 + v2 + v3 - (v0*A0/A1))


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
pressureByGamma1 = pressure1/(g*rho)
#One elbow fitting between node 1 and 2
minorLoss1 = .3

#firstOutlet
length4 = 0
elevation1 = 0
pHead = 0
d1 = 2
pressure1 = 15000
pressureByGamma1 = pressure1/(g*rho)
minorLoss1 = 30

#pipe rougness
e = 0.0000015
ed= e/d1

#############################################################################

#Kicker lines to obtain initial friction factors
re1 = reGuess
re2 = reGuess
re3 = reGuess
re4 = reGuess
ff1 = coleBrook(ed, reGuess)
ff2 = coleBrook(ed, reGuess)
ff3 = coleBrook(ed, reGuess)

nre1 = 0
nre2 = 0
nre3 = 0
nre4 = 0
Hp = 30
p = np.array([3, 3, 3, Hp])
ffs = np.array([ff1, ff2, ff3])
#Runs program with guessed friction factor until all changes in reynold numbers are small
while( abs(re1-nre1)>30 or abs(re2-nre2)>30 or abs(re3-nre3)>30):
	re1 = nre1
	re2 = nre2
	re3 = nre3
	
	v1, v2 ,v3, Hp = fsolve(equations,p,ffs)
	
	nre1 = reyNum(v1,d1,rho,mu)
	nre2 = reyNum(v2,d1,rho,mu)
	nre3 = reyNum(v3,d1,rho,mu)

	ff1 = coleBrook(ed, nre1)
	ff2 = coleBrook(ed, nre2)
	ff3 = coleBrook(ed, nre3)
	p = v1,v2,v3,Hp
	ffs = np.array([ff1,ff2,ff3])

print('velocity 1: ' + str(v1))
print('velocity 2: ' + str(v2))
print('velocity 3: ' + str(v3))
