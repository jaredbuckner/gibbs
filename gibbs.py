#!/bin/env python3
## Let's play with Gibbs Free Energy!

## DelH = DelG + T * DelS
##    where
## Del* = sum(products) - sum(reactants)

import formula

## Enthalpy of formation (kJ), Standard molar entropy(J)
data_table = {
    ##         
    'H2':     (0, 130.68),
    
    'O2':     (0, 205.15),
    'O3':     (142.67, 238.92),
    'H2O':    (285.83, 69.95),
    'H2O2':   (-187.78, 109.6),

    'C':      (0, 5.6),
    'CH2':    (386.39, 193.93),   ## Methylene
    'CH4':    (-74.87, 188.66),   ## Methane
    'C2H2':   (226.73, 200.93),   ## Acetylene
    'C2H4':   (52.47, 219.32),    ## Ethylene
    'C2H6':   (-84.0, 229.6),     ## Ethane
    'C3H8':   (-104.7, 269.91),   ## Propane
    'CO':     (-110.53, 197.66),
    'CO2':    (-393.51, 213.79),
    'C2O':    (286.60, 233.07),
    'CH2O':   (-115.9, 218.95),   ## Formaldehyde
    'CH2O2':  (-425.09, 131.84),  ## Formic acid
    'CH4O':   (-238.4, 127.19),   ## Methanol
    'C2H4O':  (-196.4, 117.3),    ## Acetaldehyde
    'C2H4O2': (-483.52, 158.0),   ## Acetic acid
    'C2H6O':  (-276.0, 159.86),   ## Ethanol
    'C3H6O':  (-249.4, 200.4),    ## Acetone
    'C3H6O2': (-510.8, 191.0),    ## Propanoic acid
    'C3H6O3': (-578.2, 180.0),    ## Glyceraldehyde, the basis of sugars!  Entropy estimated
    'C3H8O':  (-317.0, 180.58),   ## Isopropyl Alcohol
    'C3H8O3': (-669.6, 206.3),    ## Glycerol
    
    'N2':    (0, 191.6),
    'NH':    (376.56, 181.25),    ## Imidogen
    'NH3':   (-45.94, 192.77),    ## Ammonia
    'N2H4':  (50.63, 121.52),     ## Hydrazine
    'NO':    (90.29, 210.76),     ## Nitric Oxide
    'N2O':   (82.05, 219.96),     ## Nitrous Oxide
    'NO2':   (33.10, 240.04),     ## Nitrogen Dioxide
    'HNO2':  (-76.73, 249.41),    ## Nitrous Acid
    'HNO3':  (-134.31, 266.39),   ## Nitric Acid
    'C2N2':  (309.07, 241.57),    ## Cyanogen
    'HCN':   (135.14, 201.82),    ## Hydrogen Cyanide
    'CH3N':  (-47.3, 150.2),      ## Methylamine
    'HCNO':  (-101.67, 238.22)    ## Isocyanic acid
    'C2H5NO2': (-527.5, 103.5),   ## Glycine (Ammino acid)
    'CO(NH2)2': (-333.11, 104.26),  ## Urea

    'Ne':    (0, 146.33),

    'Mg':      (0, 32.67),
    'MgO':     (-601.6, 26.95),    ## Magnesia
    'Mg(OH)2': (-924.66, 63.18),   ## Brucite
    'MgCO3':   (-1111.69, 65.84),
    'Mg3N2':   (-461.08, 87.86),
    'MgS':     (-345.72, 50.30),   ## Niningerite
    'MgSO4':   (-1261.79, 91.46),  ## Epsomite (if hydrated)
    

    'Si':      (0, 18.82),
    'SiO2':    (-910.7, 41.46),    ## Quartz

    'Fe':      (0, 27.31),
    'FeO':     (-272.04, 60.75),   ## Wustite
    'Fe2O3':   (-1120.89, 145.2),  ## Magnetite
    'FeS':     (-101.67, 50.50),
    'FeS2':    (-167.36, 62.38),   ## Pyrite
    'FeSO4':   (-928.85, 120.93),  ## Melanterite (if hydrated)
    
    'S':       (0, 32.054),
    'H2S':     (-20.6, 205.81),
    'SO2':     (-296.81, 248.22),
    'SO3':     (-395.77, 256.77),
    'SO':      (5.01, 221.94),
    'S2O':     (-56.48, 266.89),
    'H2SO4':   (-814, 157),        ## Sulfates needed for life
    'COS':     (-138.41, 231.57),  ## Carbonyl sulfide

    'Ar':      (0, 154.85),

    'Al':      (0, 28.27),
    'AlO':     (66.94, 218.33),
    'AlO2':    (-86.18, 251.83),
    'Al2O':    (-145.19, 252.24),
    'Al2O2':   (-394.55, 280.90),
    'Al2O3':   (-1675.6, 50.92),   ## Corundum
    'AlN':     (-317.98, 20.14),
    
    'Ca':      (0, 41.59),
    'CaO':     (-634.92, 38.1),    ## Quicklime
    'Ca(OH)2': (-986.09, 83.36),   ## Slaked Lime
    'CaCO3':   (-1207, 93),        ## Calcite
    'CaSO4':   (-1433, 107),       ## Gypsum
    
    
    
    'P':       (-17.46, 41.09),    ## Black phosphorus.  NOTE:  The most stable state for Phosphorus is not the zero-reference state!
    'P2':      (144.0, 218.12),    ## Gas phase at high temperatures
    'P4':      (58.9, 280.01),     ## White phosphorus in gas phase
    'PH3':     (5.47, 210.24),     ## Phosphine
    'PO':      (-23.55, 222.78),
    'PO2':     (-314.52, 253.69),
    'P2O5':    (-2904.09/2, 403.96/2),  ## Note, actually a dimer (P4O10), so values are halved
    'H3PO4':   (-1271.66, 150.77)  ## Phosphoric Acid

    #### NOTE:  The connection of H3PO4 to P2O5 is something like this:
    ##
    ## P2O5 + 3 H2O <-> H2P2O6 + 2 H2O <-> H4P2O7 + H2O <-> 2 H3PO4
    ##
    #### But there isn't any data for the polymeric forms
    
    'Ti':    (0, 30.63),
    'TiO2':  (-944.7, 50.33),

    'Al':    (0, 28.33),
    'Al2O3': (-1675.7, 50.92),
    
    'Fe':    (0, 27.28),
    'Fe2O3': (-824.2, 87.40),
    'FeS':   (-100, 60.29),

    'Cu':    (0, 33.15),
    'CuO':   (-157.3, 42.63),
    'Cu2O':  (-168.6, 93.14),
    'CuSO4': (-771.36, 109),

    'Si':    (0, 18.83),
    'SiO2':  (-910, 41.84),
}

def str2reaction(s):
    reactStr, prodStr = s.split('->');

    reactFormulae, reactAtoms, dummy = formula.scanMixture(reactStr.strip())
    prodFormulae,  prodAtoms,  dummy = formula.scanMixture(prodStr.strip())

    participates = set(reactAtoms.keys())
    participates.union(prodAtoms.keys())

    for atom in participates:
        if atom not in reactAtoms:
            raise RuntimeError(f"No {atom} in the reactants")
        if atom not in prodAtoms:
            raise RuntimeError(f"No {atom} in the products")
        if reactAtoms[atom] != prodAtoms[atom]:
            raise RuntimeError(f"Different number of {atom} atoms in the reactants ({reactAtoms[atom]}) and the products ({prodAtoms[atom]})")
                
    return(reactFormulae, prodFormulae)


def react(eqStr, T=298.15):
    ## Returns (delH, delG, delQ)
    reactants, products = str2reaction(eqStr)
    
    delH = 0.0
    delS = 0.0
    for m, formula, dummy in reactants:
        formula_data = data_table[formula]
        delH -= m * formula_data[0]
        delS -= m * formula_data[1]
    
    for m, formula, dummy in products:
        formula_data = data_table[formula]
        delH += m * formula_data[0]
        delS += m * formula_data[1]

    delQ = T * delS * 0.001;
    delG = delH - delQ;

    return(delH, delG, delQ)


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description='Gibbs Free Energy')
    parser.add_argument('equation', help='A chemical equation with reactants and products')
    parser.add_argument('--temp', type=float, help='Temperature of reaction', default=298.15)

    args = parser.parse_args()

    delH, delG, delQ = react(args.equation, T=args.temp)

    print(f'{delH=:.4g}kJ, {delG=:.4g}kJ, {delQ=:.4g}kJ')

    exit(0)
