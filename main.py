import numpy as np
import math
from pint import UnitRegistry

UR = UnitRegistry(system='mks')
Q_ = UR.Quantity

siEnergy = "joule/kilogram/degK"
enEnergy = "Btu/pound/degR"

g = Q_(32.2,"feet/second**2")

# Engine operating conditions
class Engine:
    pc = Q_(400,"pound/inch**2")
    cstar = Q_(5000,"feet/second")
    Dt = Q_(5,"centimeter") #throat diameter
    dt = Dt.to("inch")
    Dc = Q_(10,'centimeter')    #chamber diameter
    R = Q_(1,"centimeter")  #Radius of throat curvature
    r = R.to("inch")
    At_Ac = 0.1
    wallThick = 5 * UR.millimeter
    wallK = Q_(16.3,"watt/meter/degK")

# Combustion gas values taken from CEA
class Comb:
    Visc = Q_(3,"Pa*second")    #CEA combustion viscosity
    visc = Visc.to("pound/inch/second")

    Cond = Q_(3.27,"milliwatts/centimeter/degK")    #CEA combustion conductivity
    cond = Cond.to("Btu/second/inch**2/degF*inch")

    Cp = Q_(2.23e3,"joule/kilogram/degK")     #CEA combustion Cp
    cp = Cp.to("Btu/pound/degF")
    Pr = 0.525  #Combustion Pr
    Mc = 0.1    #Arbitrary chamber Mach (fix later)
    gam = 1.25  #Combustion gas gamma

    Twg_Tc = 1  #Arbitrary hot wall to combustion temp ratio
    

# Coolant (RP-1)
class Coolant:

    cond = Q_(0.078,"Btu/feet/hour/degF")
    Cond = cond.to("watt/meter/degK")
    visc = Q_(100e5,"pound/feet/second")
    Visc = visc.to("pascal*second")
    cp = Q_(0.49,"Btu/pound/degF")
    Cp = cp.to("joule/kilogram/degK")

    Pr = Visc*Cp/Cond

    Hc = Q_(300,"watt/meter**2/degK")
    mdot = 1 * UR.kilogram / UR.second

dx = 0.01 * UR.meter


As = math.pi*Engine.Dc*dx
As = As.to_base_units()


# Analysis
sig = 1/((0.5*Comb.Twg_Tc*(1 + ((Comb.gam-1)/2)*Comb.Mc**2) + 0.5)**0.68*(1 + ((Comb.gam-1)/2)*Comb.Mc**2)**0.12)

hg = (0.026/Engine.dt**0.2)*(Comb.visc**0.2*Comb.cp/Comb.Pr**0.6)*(Engine.pc*g/Engine.cstar)**0.8*(Engine.dt/Engine.r)**0.1*(Engine.At_Ac)**0.9*sig

Hg = hg.to("watt/meter**2/degK")


Tc0 = 300 * UR.degK

T0 = 2500 * UR.degK

Tc = Tc0

T = np.array(Tc)
for n in range(20):

    R = 1/(Hg*As) + (Engine.wallThick/(Engine.wallK*As)) + 1/(Coolant.Hc*As)

    q = (T0-Tc)/R

    Tout = q/(Coolant.mdot*Coolant.Cp) + Tc
    Tout = Tout.to_base_units()

    T = np.append(T,Tout / UR.degK)

    Tc = Tout

print(T)


    