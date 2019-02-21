#
#
#
#
from Sire.IO import *
from Sire.Mol import *
from Sire.MM import *
from Sire.System import *
from Sire.Units import *
from Sire.CAS import *  
from Sire.Maths import *

if __name__ == '__main__':
    # 1) Read a pdb file describing the system to simulate
    pdbfile = 'pyridine/MOL.pdb'
    p = PDB2(pdbfile)
    s = p.toSystem()
    molecules = s.molecules()
    print (molecules)


     # 2) Now we read the xml file, and store parameters for each molecule


    import xml.dom.minidom as minidom
    xmldoc = minidom.parse('pyridine/MOL_extra.xml')
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #~~~~~~~~~~~~~~~~~~~~~~ TAG NAME: TYPE ~~~~~~~~~~~~~~~~~~~~~
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    itemlist_type = xmldoc.getElementsByTagName('Type')
    dicts_type = []
    for items in itemlist_type:
        d = {}
        for a in items.attributes.values():
            d[a.name] = a.value
        dicts_type.append(d)
    dicts_tp =  str(dicts_type).split()
    print (dicts_tp)


    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #~~~~~~~~~~~~~~~~~~~~~~ TAG NAME: ATOM ~~~~~~~~~~~~~~~~~~~~~
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    itemlist_atom = xmldoc.getElementsByTagName('Atom')
    dicts_atom = []
    for items in itemlist_atom:
        d = {}
        for a in items.attributes.values():
            d[a.name] = a.value
        dicts_atom.append(d)
    dicts_at =  str(dicts_atom).split()
    print (dicts_at)

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #~~~~~~~~~~~~~~~~~~~~~~ TAG NAME: BOND ~~~~~~~~~~~~~~~~~~~~~
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    itemlist_bond = xmldoc.getElementsByTagName('Bond')
    dicts_bond = []
    for items in itemlist_bond:
        d = {}
        for a in items.attributes.values():
            d[a.name] = a.value
        dicts_bond.append(d)
    dicts_b =  str(dicts_bond).split()
    print (dicts_b)

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #~~~~~~~~~~~~~~~~~~~~~~ TAG NAME: ANGLE ~~~~~~~~~~~~~~~~~~~~
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    itemlist_angle = xmldoc.getElementsByTagName('Angle')
    dicts_angle = []
    for items in itemlist_angle:
        d = {}
        for a in items.attributes.values():
            d[a.name] = a.value
        dicts_angle.append(d)
    dicts_ang =  str(dicts_angle).split()
    print (dicts_ang)


    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #~~~~~~~~~~~~~~~~~~~~~~ TAG NAME: PROPER ~~~~~~~~~~~~~~~~~~~
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    itemlist_proper = xmldoc.getElementsByTagName('Proper')
    dicts_proper = []
    for items in itemlist_proper:
        d = {}
        for a in items.attributes.values():
            d[a.name] = a.value
        dicts_proper.append(d)
    dicts_pr =  str(dicts_proper).split()
    print (dicts_pr)



    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #~~~~~~~~~~~~~~~~~~~~~~ TAG NAME: IMPROPER ~~~~~~~~~~~~~~~~~
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    itemlist_improper = xmldoc.getElementsByTagName('Improper')
    dicts_improper = []
    for items in itemlist_improper:
        d = {}
        for a in items.attributes.values():
            d[a.name] = a.value
        dicts_improper.append(d)
    dicts_impr =  str(dicts_improper).split()
    print (dicts_impr)


    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #~~~~~~~~~~~~~~~~~~~~~~ TAG NAME: VIRTUAL SITES ~~~~~~~~~~~~
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    itemlist_VirtualSite = xmldoc.getElementsByTagName('VirtualSite')
    dicts_virtualsite = []
    for items in itemlist_VirtualSite:
        d = {}
        for a in items.attributes.values():
            d[a.name] = a.value
        dicts_virtualsite.append(d)
    dicts_vs =  str(dicts_virtualsite).split()
    #print (dicts_vs)


    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #~~~~~~~~~~~~~~~~~~~~~~ TAG NAME: RESIDUE ~~~~~~~~~~~~~~~~
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    itemlist_residue = xmldoc.getElementsByTagName('Residue')
    dicts_residue = []
    for items in itemlist_residue:
        d = {}
        for a in items.attributes.values():
            d[a.name] = a.value
        dicts_residue.append(d)
    dicts_res =  str(dicts_residue).split()
    print (dicts_res)

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #~~~~~~~~~~~~~ TAG NAME: NON BONDED FORCE ~~~~~~~~~~~~~~~~
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    itemlist_nonbond = xmldoc.getElementsByTagName('NonbondedForce')
    dicts_nonb = []
    for items in itemlist_nonbond:
        d = {}
        for a in items.attributes.values():
            d[a.name] = a.value
        dicts_nonb.append(d)
    dicts_nb =  str(dicts_nonb).split()
    print (dicts_nb)


    # 3) Now we create an Amberparameters object for each molecule
    molnums = molecules.molNums()

    newmolecules = Molecules()
    for molnum in molnums:
        mol = molecules.at(molnum)
        print (mol) #Molecule( 1 version 9 : nAtoms() = 11, nResidues() = 1 )
        mol_params = AmberParameters() #SireMol::AmberParameters()
        
        # We populate the Amberparameters object with a list of bond, angle, dihedrals
        # We look up parameters from the contents of the xml file
        # We also have to set the atomic parameters (q, sigma, epsilon)
        editmol = mol.edit()
        atoms = editmol.atoms()
        # We update atom parameters see setAtomParameters in SireIO/amber.cpp l2122 
        natoms = editmol.nAtoms()
        print("number of atoms is %s" %natoms)
        nVirtualSites = xmldoc.getElementsByTagName('VirtualSite').length #number of Virtual Sites
        
        #atoms doesn't include the virtual sites! 
        
        for atom in atoms: 
            editatom = editmol.atom(atom.index())
        
            print("The editatom object is %s"%editatom)
            
            i = int(str(atom.number()).split('(')[1].replace(")" , " "))
            editatom.setProperty("charge", float(dicts_atom[i+natoms+nVirtualSites]['charge']) * mod_electron)
            editatom.setProperty("mass", float(dicts_type[i-1]['mass']) * g_per_mol) #
            editatom.setProperty("LJ", LJParameter( float(dicts_atom[i+natoms+nVirtualSites]['sigma']) * angstrom, float(dicts_atom[i+natoms+nVirtualSites]['epsilon']) * kcal_per_mol))
            editatom.setProperty("ambertype", dicts_atom[i+natoms+nVirtualSites]['type'])
           
            editmol = editatom.molecule()
            print(editmol)

        # Now we create a connectivity see setConnectivity in SireIO/amber.cpp l2144
        # XML data tells us how atoms are bonded in the molecule (Bond 'from' and 'to')
        #bondfuncs = k_bondList
        #new_bond = editmol.setProperty("bonds_new", bondfuncs).commit()
        if natoms > 1:
            print("Set up connectivity")

            c = []
            for atom in atoms:
                i = int(str(atom.number()).split('(')[1].replace(")" , " "))
                if natoms > 1: #if numebr of atoms >1  
                    connect_prop= {}
                    connect_prop = dicts_bond[i-1]['from'], dicts_bond[i-1]['to']
                c.append(connect_prop)
            print(c)

            conn = Connectivity(editmol.info()).edit()

            for j in range(0,natoms):
                conn.connect(atoms[int(c[j][0]) ].index(), atoms[int(c[j][1]) ].index()) 
                conn.commit()   
            print(conn.commit)               
            
            for atom in atoms:
                #i = int(str(atom.number()).split('(')[1].replace(")" , " "))
                editmol.setProperty("connectivity", conn)
                print(editmol.setProperty("connectivity", conn))
            

             # Now we add bond parameters to the Sire molecule. We also update amberparameters see SireIO/amber.cpp l2154

            internalff = InternalFF()
                    
            bondfuncs = TwoAtomFunctions(mol)
            r = internalff.symbols().bond().r()

            for j in range(0,natoms):
                bf = bondfuncs.set(atoms[int(c[j][0]) ].index(), atoms[int(c[j][1]) ].index(), float(dicts_bond[j+natoms]['k'])* (float(dicts_bond[j+natoms]['length']) - r) **2  )
                
            mol = editmol.setProperty("bond", bondfuncs).commit()

        # Now we add angle parameters to the Sire molecule. We also update amberparameters see SireIO/amber.cpp L2172
        if natoms > 2:
            print("Set up angles")

            anglefuncs = ThreeAtomFunctions(mol)

            nAngles= xmldoc.getElementsByTagName('Angle').length
            at1 = []
            for i in range(0, nAngles):
                for j in range(0,natoms):
                    if dicts_angle[i]['class1']  == dicts_type[j]['class']:
                        a1 = {}
                        a1 = j
                        at1.append(a1)
            print (at1)

            at2 = []
            for i in range(0, nAngles):
                for j in range(0,natoms):
                    if dicts_angle[i]['class2']  == dicts_type[j]['class']:
                        a2 = {}
                        a2 = j
                        at2.append(a2)
            print (at2)

            at3 = []
            for i in range(0, nAngles):
                for j in range(0,natoms):
                    if dicts_angle[i]['class3']  == dicts_type[j]['class']:
                        a3 = {}
                        a3 = j
                        at3.append(a3)
            print (at3)

            theta = internalff.symbols().angle().theta()
            for j in range(0,natoms):
                anglefuncs.set( atoms[at1[j]].index(), atoms[at2[j]].index(), atoms[at3[j]].index(), float(dicts_angle[j]['k']) * ( (float(dicts_angle[j]['angle']) *degrees.value() - theta )**2 ))

            mol = editmol.setProperty("angle" , anglefuncs).commit()
             
        # Now we add dihedral parameters to the Sire molecule. We also update amberparameters see SireIO/amber.cpp L2190

        if natoms > 3:
            print("Set up dihedrals")

            amber_dihedral = AmberDihedral(dihedral.function(), Symbol("phi"))

            nProper= xmldoc.getElementsByTagName('Proper').length
            di1 = []
            for i in range(0, nProper):
                for j in range(0,natoms):
                    if dicts_proper[i]['class1']  == dicts_type[j]['class']:
                        d1 = {}
                        d1 = j
                        di1.append(d1)
            print (di1)


            di2 = []
            for i in range(0, nProper):
                for j in range(0,natoms):
                    if dicts_proper[i]['class2']  == dicts_type[j]['class']:
                        d2 = {}
                        d2 = j
                        di2.append(d1)
            print (di2)


            di3 = []
            for i in range(0, nProper):
                for j in range(0,natoms):
                    if dicts_proper[i]['class3']  == dicts_type[j]['class']:
                        d3 = {}
                        d3 = j
                        di3.append(d3)
            print (di3)


            
            di4 = []
            for i in range(0, nProper):
                for j in range(0,natoms):
                    if dicts_proper[i]['class4']  == dicts_type[j]['class']:
                        d4 = {}
                        d4 = j
                        di4.append(d4)
            print (di4)
    

            dihedralfuncs = FourAtomFunctions(mol)
    
            phi = internalff.symbols().dihedral().phi()
            for i in range(1,5):    
                for j in range(0,nProper):
                    dihedralfuncs.set(atoms[di1[j]].index(), atoms[di2[j]].index(),atoms[di3[j]].index(),atoms[di4[j]].index(),\
                     float(dicts_proper[j]['k%s'%i])*(1+ Cos(int(dicts_proper[j]['periodicity%s'%i]) * phi - float(dicts_proper[j]['phase%s'%i]))))
            mol = editmol.setProperty("dihedral" , dihedralfuncs).commit()

            print("Set up impropers")
            
            di_im1 = []
            for i in range(0, nImproper):
                for j in range(0,natoms):
                    if dicts_improper[i]['class1']  == dicts_type[j]['class']:
                        d1 = {}
                        d1 = j
                        di_im1.append(d1)
            print (di_im1)


            di_im2 = []
            for i in range(0, nImproper):
                for j in range(0,natoms):
                    if dicts_improper[i]['class2']  == dicts_type[j]['class']:
                        d2 = {}
                        d2 = j
                        di_im2.append(d2)
            print (di_im2)


            di_im3 = []
            for i in range(0, nImproper):
                for j in range(0,natoms):
                    if dicts_improper[i]['class3']  == dicts_type[j]['class']:
                        d3 = {}
                        d3 = j
                        di_im3.append(d3)
            print (di_im3)


            
            di_im4 = []
            for i in range(0, nImproper):
                for j in range(0,natoms):
                    if dicts_improper[i]['class4']  == dicts_type[j]['class']:
                        d4 = {}
                        d4 = j
                        di_im4.append(d4)
            print (di_im4)


            improperfuncs = FourAtomFunctions(mol)
    
            phi_im = internalff.symbols().improper().phi()
            for i in range(1,5):    
                for j in range(0,nImproper):
                    improperfuncs.set(atoms[di_im1[j]].index(), atoms[di_im2[j]].index(),atoms[di_im3[j]].index(),atoms[di_im4[j]].index(),\
                     float(dicts_improper[j]['k%s'%i])*(1+ Cos(int(dicts_improper[j]['periodicity%s'%i]) * phi_im - float(dicts_improper[j]['phase%s'%i]))))
            mol = editmol.setProperty("improper" , improperfuncs).commit()




        # Now we work out non bonded pairs see SireIO/amber.cpp L2213
    


        molecule = editmol.commit()
        newmolecules.add(molecule)
    # By the end of this loop we have a new set of mol that looks
    # exactly like a molecules object returned by AMber().readCrdTop(...)




