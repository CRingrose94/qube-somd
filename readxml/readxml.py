#
#
#
#
from Sire.IO import *
from Sire.Mol import *
from Sire.System import *


if __name__ == '__main__':
    # 1) Read a pdb file describing the system to simulate
    pdbfile = 'pyridine/MOL.pdb'
    p = PDB2(pdbfile)
    s = p.toSystem()
    molecules = s.molecules()
    print (molecules)
    # 2) Now we read the xml file, and store parameters for each molecule
    # TODO
    # 3) Now we create an Amberparameters object for each molecule
    molnums = molecules.molNums()

    for molnum in molnums:
        mol = molecules.at(molnum)
        print (mol)
        mol_params = AmberParameters()
        # We populate the Amberparameters object with a list of bond, angle, dihedrals
        # We look up parameters from the contents of the xml file
        # We also have to set the atomic parameters (q, sigma, epsilon)

        # We also add connectivity, bond, angles, dihedrals etc ... to mol
        # We end up adding mol_params to mol

    # By the end of this loop we have a new set of mol that looks
    # exactly like a molecules object returned by AMber().readCrdTop(...)
