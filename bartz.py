from pint import UnitRegistry, get_application_registry

UR = UnitRegistry(system='mks')
UR = get_application_registry()

def Bartz(EngineDt,CombVisc,CombCp,CombPr,Pc,Mc,Gam,EngineR,cstar,At_Ac,Twg_Tc):

    g = 32.2 * UR.feet / UR.second**2

    #Converting inputs to correct units
    EngineDt = EngineDt.to("inch")
    CombVisc = CombVisc.to("pound/inch/second")
    CombCp = CombCp.to("Btu/pound/degF")
    Pc = Pc.to("pound/inch**2")
    EngineR = EngineR.to("inch")
    cstar = cstar.to("feet/second")

    #Bartz equation calc
    sig = 1/((0.5*Twg_Tc*(1 + ((Gam-1)/2)*Mc**2) + 0.5)**0.68*(1 + ((Gam-1)/2)*Mc**2)**0.12)

    hg = (0.026/EngineDt**0.2)*(CombVisc**0.2*CombCp/CombPr**0.6)*(Pc*g/cstar)**0.8*(EngineDt/EngineR)**0.1*(At_Ac)**0.9*sig

    #Convert to usefull units
    hg = hg.to("watt/meter**2/degK")

    return hg