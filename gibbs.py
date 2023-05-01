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
    'C2H5NO2': (-527.5, 103.5),   ## Glycine (Ammino acid)

    'Ne':    (0, 146.33),

    'Mg':      (0, 32.67),
    'MgO':     (-601.6, 26.95),    ## Magnesia
    'Mg(OH)2': (-924.66, 63.18),   ## Brucite
    'MgCO3':   (-1111.69, 65.84),
    'Mg3N2':   (-461.08, 87.86),

    'Si':      (0, 18.82),
    'SiO2':    (-910.7, 41.46),    ## Quartz

    'Fe':      (0, 27.31),
    'FeO':     (-272.04, 60.75),   ## Wustite
    'Fe2O3':   (-1120.89, 145.2),  ## Magnetite
    
    'S':       (0, 32.054),
    'H2S':     (-20.6, 205.81),
    'SO2':     (-296.81, 248.22),
    'SO3':     (-395.77, 256.77),
    'SO':      (5.01, 221.94),
    'S2O':     (-56.48, 266.89),
    'H2SO4':   (-814, 157),
    
    

    'P':     (0, 41.09),
    'PH3':   (5.4, 210),
    'H3PO4': (-1279.0, 110.5),  ## Tri-phosphates needed for life
    
    'O2':    (0, 205.03),
    'H2O':   (-285.830, 69.91),
    'H2O2':  (-187.78, 109.6),

    'S':     (0, 31.8),
    'H2S':   (-20.63, 205.68),
    'SO2':   (-296.83, 248.11),
    'SO3':   (-395.72, 256.65),
    'H2SO4': (-813.99, 156.90),

    'C6H12O6':  (-1271, 209.2),   ## Glucose  (Standin for sugar)
    'C2H5NO2':  (-528.61, 103.51),  ## Glycine  (Amino acid, standin for protein)

    'Mg':    (0, 32.68),
    'MgO':   (-601.7, 26.94),

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
