import numpy as np
from pint import UnitRegistry

UR = UnitRegistry()
Q_ = UR.Quantity

siEnergy = "joule/kilogram/degK"
enEnergy = "Btu/pound/degR"

Pr = 0.525
Mc = 0.1
gam = 1.25
Twg_Tc = 1
At_Ac = 0.1

pc = Q_(400,"pound/inch**2")
cstar = Q_(5000,"feet/second")

Dt = Q_(5,"centimeter") #throat diameter
dt = Dt.to("inch")

R = Q_(1,"centimeter")
r = R.to("inch")

Visc = Q_(3,"Pa*second")
visc = Visc.to("pound/inch/second")

Cond = Q_(3.27,"milliwatts/centimeter/degK")
cond = Cond.to("Btu/second/inch**2/degF*inch")

Cp = Q_(2.23,"kilojoule/kilogram/degK")
cp = Cp.to("Btu/pound/degF")

g = Q_(32.2,"feet/second**2")

sig = 1/((0.5*Twg_Tc*(1 + ((gam-1)/2)*Mc**2) + 0.5)**0.68*(1 + ((gam-1)/2)*Mc**2)**0.12)

hg = (0.026/dt**0.2)*(visc**0.2*cp/Pr**0.6)*(pc*g/cstar)**0.8*(dt/r)**0.1*(At_Ac)**0.9*sig

Hg = hg.to("watt/meter**2/degK")

print(Hg)