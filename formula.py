import re
import unittest

def scanLiteral(s, lit, idx=0):
    end = idx+len(lit)
    if(s[idx:end] == lit):
        return(lit, end)
    else:
        return(None, idx)

spacePat = '(?:\s+)'
spaceRx  = re.compile(spacePat)
    
def scanSpace(s, idx=0):
    match = spaceRx.match(s, idx)
    if(match):
        return(match.group(), match.end())
    else:
        return(None, idx)

countPat = '(?:([0-9]+))'
countRx  = re.compile(countPat)

def scanCount(s, idx=0):
    match = countRx.match(s, idx)
    if(match):
        return(int(match.group()), match.end())
    else:
        return(1, idx)

    
formulaPartPat = '(?:(?:([A-Z][a-z]?)([0-9]*)))'
formulaPartRx  = re.compile(formulaPartPat)

def scanFormula(s, idx=0):
    formula = ""
    atoms   = {}

    match = formulaPartRx.match(s, idx)

    while(match):
        formula += match.group(0)
        atom = match.group(1)
        count = match.group(2)
        idx = match.end()
        
        if count == "":
            count = 1
        else:
            count = int(count)

        atoms.setdefault(atom, 0)
        atoms[atom] += count

        match = formulaPartRx.match(s, idx)
        
        
    if not formula:
        return (None, None, idx)
    else:
        return(formula, atoms, idx)


def scanMixture(s, idx=0):
    formulae = []
    atoms    = {}

    
    molMul, jdx = scanCount(s, idx)
    dummy, jdx = scanSpace(s, jdx)
    molForm, molAtoms, jdx = scanFormula(s, jdx)
    
    if(molForm is None):
        return(None, None, idx)
    
    formulae.append((molMul, molForm, molAtoms))
    for atom, count in molAtoms.items():
        atoms[atom] = molMul * count
    
    idx = jdx
    while(True):
        dummy, jdx = scanSpace(s, idx)
        sep,   jdx = scanLiteral(s, "+", jdx)
        if sep is None:
            return(formulae, atoms, idx)

        dummy, jdx = scanSpace(s, jdx)
        molMul, jdx = scanCount(s, jdx)
        dummy, jdx = scanSpace(s, jdx)
        molForm, molAtoms, jdx = scanFormula(s, jdx)

        if(molForm is None):
            return(formulae, atoms, idx)

        formulae.append((molMul, molForm, molAtoms))
        for atom, count in molAtoms.items():
            if atom in atoms:
                atoms[atom] += molMul * count
            else:
                atoms[atom] = molMul * count

        idx = jdx



class utFormula(unittest.TestCase):
    def testScanFormula(self):
        formula, atoms, idx = scanFormula("")
        self.assertIsNone(formula)
        self.assertIsNone(atoms)
        self.assertEqual(idx, 0)

        formula, atoms, idx = scanFormula("H")
        self.assertEqual(formula, "H")
        self.assertEqual(atoms, {"H": 1})
        self.assertEqual(idx, 1)

        formula, atoms, idx = scanFormula("H2")
        self.assertEqual(formula, "H2")
        self.assertEqual(atoms, {"H": 2})
        self.assertEqual(idx, 2)

        formula, atoms, idx = scanFormula("He")
        self.assertEqual(formula, "He");
        self.assertEqual(atoms, {"He": 1})
        self.assertEqual(idx, 2)

        formula, atoms, idx = scanFormula("MgCl2")
        self.assertEqual(formula, "MgCl2")
        self.assertEqual(atoms, {"Mg": 1, "Cl": 2})
        self.assertEqual(idx, 5)
        
        formula, atoms, idx = scanFormula("H2O+CO2")
        self.assertEqual(formula, "H2O")
        self.assertEqual(atoms, {"H": 2, "O": 1})
        self.assertEqual(idx, 3)
                         
    def testScanMixture(self):
        formulae, atoms, idx = scanMixture("")
        self.assertIsNone(formulae)
        self.assertIsNone(atoms)
        self.assertEqual(idx, 0)
        
        formulae, atoms, idx = scanMixture("H2SO4")
        self.assertEqual(formulae, [(1, "H2SO4", {"H": 2, "S": 1, "O": 4})])
        self.assertEqual(atoms, {"H": 2, "S": 1, "O": 4})
        self.assertEqual(idx, 5)

        formulae, atoms, idx = scanMixture("2H2 + O2")
        self.assertEqual(formulae, [(2, "H2", {"H": 2}), (1, "O2", {"O": 2})])
        self.assertEqual(atoms, {"H": 4, "O": 2})
        self.assertEqual(idx, 8)
