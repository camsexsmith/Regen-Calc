import numpy as np
import math
from pint import UnitRegistry, get_application_registry
from bartz import Bartz

UR = UnitRegistry(system='mks')
UR = get_application_registry()
Q_ = UR.Quantity

g = Q_(32.2,"feet/second**2")

# Engine operating conditions
class Engine:
    pc = Q_(400,"pound/inch**2")
    cstar = Q_(2000,"meter/second")
    Dt = Q_(5,"centimeter") #throat diameter
    Dc = Q_(10,'centimeter')    #chamber diameter
    R = Q_(1,"centimeter")  #Radius of throat curvature
    At_Ac = 0.1
    wallThick = 5 * UR.millimeter
    wallK = Q_(16.3,"watt/meter/degK")

    chamberLength = 10 * UR.centimeter

# Combustion gas values taken from CEA
class Comb:
    visc = Q_(3,"Pa*second")    #CEA combustion viscosity

    cond = Q_(3.27,"milliwatts/centimeter/degK")    #CEA combustion conductivity

    cp = Q_(2.23e3,"joule/kilogram/degK")     #CEA combustion Cp
    Pr = 0.525  #Combustion Pr
    Mc = 0.1    #Arbitrary chamber Mach (fix later)
    gam = 1.25  #Combustion gas gamma

    Twg_Tc = 1  #Arbitrary hot wall to combustion temp ratio
    

# Coolant (RP-1) data taken from Aerojet - Properties and Performance of Liquid Rocket Propellants
class Coolant:

    cond = Q_(0.078,"Btu/feet/hour/degF")
    visc = Q_(100e-5,"pound/feet/second")
    cp = Q_(0.49,"Btu/pound/degF")

    Pr = (visc*cp/cond).to_base_units()

    hc = Q_(300,"watt/meter**2/degK")
    mdot = 1 * UR.kilogram / UR.second


dx = 0.01 * UR.meter

As = (math.pi*Engine.Dc*dx).to_base_units()


#Bartz Equation
hg = Bartz(Engine.Dt,Comb.visc,Comb.cp,Comb.Pr,Engine.pc,Comb.Mc,Comb.gam,Engine.R,Engine.cstar,Engine.At_Ac,Comb.Twg_Tc)

print(hg)


#Analysis
Tc0 = 300 * UR.degK

T0 = 2500 * UR.degK

Tc = Tc0

T = np.array(Tc)

#Range here represents cell number (0-19) with a delta length of dx meters
for n in range(20):

    R = 1/(hg*As) + (Engine.wallThick/(Engine.wallK*As)) + 1/(Coolant.hc*As)

    q = (T0-Tc)/R

    Tout = q/(Coolant.mdot*Coolant.cp) + Tc
    Tout = Tout.to_base_units()

    T = np.append(T,Tout / UR.degK)

    Tc = Tout

print(T)
