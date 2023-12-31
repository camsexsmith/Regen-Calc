import numpy as np
from pint import UnitRegistry

UR = UnitRegistry()
Q_ = UR.Quantity

siEnergy = "joule/kilogram/degK"
enEnergy = "Btu/pound/degR"

Pr = 0.6

Dt = Q_(5,"centimeter") #throat diameter
dt = Dt.to("inch")

Visc = Q_(3,"Pa*second")
visc = Visc.to("pound/inch/second")

Cond = Q_(2,"watts/meter/degK")
cond = Cond.to("Btu/second/inch**2/degF*inch")

Cp = Q_(1,"kilojoule/kilogram/degK")
cp = Cp.to("Btu/pound/degF")