#
#
#
#
from Sire.IO import *
from Sire.Mol import *
from Sire.MM import *
from Sire.System import *
from Sire.Units import *


if __name__ == '__main__':
    # 1) Read a pdb file describing the system to simulate
    pdbfile = 'pyridine/MOL.pdb'
    p = PDB2(pdbfile)
    s = p.toSystem()
    molecules = s.molecules()
    print (molecules)
    # 2) Now we read the xml file, and store parameters for each molecule
    # TODO by Sofia
    # 3) Now we create an Amberparameters object for each molecule
    molnums = molecules.molNums()

    newmolecules = Molecules()
    for molnum in molnums:
        mol = molecules.at(molnum)
        print (mol)
        mol_params = AmberParameters()
        # We populate the Amberparameters object with a list of bond, angle, dihedrals
        # We look up parameters from the contents of the xml file
        # We also have to set the atomic parameters (q, sigma, epsilon)
        editmol = mol.edit()
        atoms = editmol.atoms()
        # We update atom parameters see setAtomParameters in SireIO/amber.cpp l2122
        for atom in atoms:
            editatom = editmol.atom(atom.index())
            editatom.setProperty("charge", -1 * mod_electron) # This should be the value for that atom read from the xml file in
            editatom.setProperty("mass", 1.5 * g_per_mol) #
            editatom.setProperty("LJ", LJParameter( 0.5 * angstrom, 0.5 * kcal_per_mol))
            editatom.setProperty("ambertype", "aa")
            editmol = editatom.molecule()
        # Now we create a connectivity see setConnectivity in SireIO/amber.cpp l2144
        # XML data tells us how atoms are bonded in the molecule (Bond 'from' and 'to')

        # Now we add bond parameters to the Sire molecule. We also update amberparameters see SireIO/amber.cpp l2154

        # Now we add angle parameters to the Sire molecule. We also update amberparameters see SireIO/amber.cpp L2172

        # Now we add dihedral parameters to the Sire molecule. We also update amberparameters see SireIO/amber.cpp L2190

        # Now we work out non bonded pairs see SireIO/amber.cpp L2213
        
        molecule = editmol.commit()
        newmolecules.add(molecule)
    # By the end of this loop we have a new set of mol that looks
    # exactly like a molecules object returned by AMber().readCrdTop(...)
