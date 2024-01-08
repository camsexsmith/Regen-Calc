import numpy as np
import math
from pint import UnitRegistry, get_application_registry
from bartz import Bartz
import matplotlib.pyplot as plt



UR = UnitRegistry(system='mks')
UR = get_application_registry()
UR.setup_matplotlib(True)
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
    rho = Q_(800,"kilogram/meter**3")
    mdot = Q_(1,"kilogram/second")

    Pr = (visc*cp/cond).to_base_units()

    passageID = (Engine.Dc + Engine.wallThick).to_base_units()
    passageThick = 1 * UR.centimeter

    passageOD = (passageID + passageThick*2).to_base_units()

    print(passageOD)

    hc = Q_(300,"watt/meter**2/degK")




cells = 100
dx = Engine.chamberLength/cells

lenVector = np.arange(0,Engine.chamberLength.to_base_units()/UR.meter,dx.to_base_units()/UR.meter) * UR.meter


#Discritized chamber surface area
As = (math.pi*Engine.Dc*dx).to_base_units()
print(As)

#Bartz Equation
hg = Bartz(Engine.Dt,Comb.visc,Comb.cp,Comb.Pr,Engine.pc,Comb.Mc,Comb.gam,Engine.R,Engine.cstar,Engine.At_Ac,Comb.Twg_Tc)

#Analysis
Tc0 = 300 * UR.degK

T0 = 2500 * UR.degK

Tcoolant = Tc0

TcoolantVec = np.array(Tcoolant)
TgwVec = np.array(Tcoolant)
TcwVec = np.array(Tcoolant)
#Range here represents cell number (0-19) with a delta length of dx meters
for n in lenVector:

    R = 1/(hg*As) + (Engine.wallThick/(Engine.wallK*As)) + 1/(Coolant.hc*As)

    q = (T0-Tcoolant)/R

    Tout = q/(Coolant.mdot*Coolant.cp) + Tcoolant
    Tout = Tout.to_base_units()

    TcoolantVec = np.append(TcoolantVec,Tout / UR.degK)

    Tcoolant = Tout

    Tgw = T0 - q/(As*hg)

    TgwVec = np.append(TgwVec,Tgw / UR.degK)

    Tcw = Tgw - Engine.wallThick*q/(Engine.wallK*As)

    TcwVec = np.append(TcwVec,Tcw / UR.degK)


TcoolantVec = TcoolantVec[1:] * UR.degK   #Removing initial cell and giving units back
TgwVec = TgwVec[1:] * UR.degK
TcwVec = TcwVec[1:] * UR.degK


fig, ax = plt.subplots()
ax.plot(lenVector,TcwVec)
ax.xaxis.set_units(UR.centimeter)
plt.grid(True)
plt.show()