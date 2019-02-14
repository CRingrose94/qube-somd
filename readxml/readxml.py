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
    #print (dicts_tp)

    name_of_typeList =[]
    for i in range (len(dicts_tp)-1):
        l = {}
        if dicts_tp[i-1] == "{'name':" or dicts_tp[i-1] == "[{'name':":
            l= dicts_tp[i].replace("'," , " ").replace("," , " ").replace("}" , " ").replace("'" , " ")
            name_of_typeList.append(l)
    print ("The name list of the type is:")
    print (name_of_typeList)


    class_of_typeList =[]
    for i in range (len(dicts_tp)-1):
        l = {}
        if dicts_tp[i-1] == "'class':" :
            l= dicts_tp[i].replace("'," , " ").replace("," , " ").replace("}" , " ").replace("'" , " ")
            class_of_typeList.append(l)
    print ("The class list of the type is:")
    print (class_of_typeList)

    element_of_typeList =[]
    for i in range (len(dicts_tp)-1):
        l = {}
        if dicts_tp[i-1] == "'element':":
            l= dicts_tp[i].replace("'," , " ").replace("," , " ").replace("}" , " ").replace("'" , " ")
            element_of_typeList.append(l)
    print ("The element list of the type is:")
    print (element_of_typeList)

    mass_of_typeList =[]
    for i in range (len(dicts_tp)-1):
        l = {}
        if dicts_tp[i-1] == "'mass':":
            l= float(dicts_tp[i].replace("'," , " ").replace("," , " ").replace("}" , " ").replace("'" , " "))
            mass_of_typeList.append(l)
    print ("The mass list of the type is:")
    print (mass_of_typeList)
    
    
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
    #print (dicts_at)

    nameList =[]
    for i in range (len(dicts_at)-1):
        name = {}
        if dicts_at[i-1] == "{'name':" or dicts_at[i-1] == "[{'name':":
            name= dicts_at[i].replace("'," , " ").replace("," , " ").replace("}" , " ").replace("'" , " ")
            nameList.append(name)
    print ("The name list is:")
    print (nameList)

    typeList =[]
    for i in range (len(dicts_at)-1):
        t = {}
        if dicts_at[i-1] == "{'type':":
            t= dicts_at[i].replace("'," , " ").replace("," , " ").replace("}" , " ").replace("'" , " ")
            typeList.append(t)
    print ("The type list is:")
    print (typeList)


    chargeList =[]
    for i in range (len(dicts_at)-1):
        t = {}
        if dicts_at[i-1] == "'charge':":
            t= float(dicts_at[i].replace("'," , " ").replace("," , " ").replace("}" , " ").replace("'" , " "))
            chargeList.append(t)
    print ("The charge list is:")
    print (chargeList)

    sigmaList =[]
    for i in range (len(dicts_at)-1):
        s = {}
        if dicts_at[i-1] == "'sigma':":
            s= float(dicts_at[i].replace("'," , " ").replace("," , " ").replace("}" , " ").replace("'" , " "))
            sigmaList.append(s)
    print ("The sigma list is:")
    print (sigmaList)


    epsilonList =[]
    for i in range (len(dicts_at)):
        e = {}
        if dicts_at[i-1] == "'epsilon':":
            e= float(dicts_at[i].replace("'," , " ").replace("," , " ").replace("}" , " ").replace("'" , " ").replace("]" , " "))
            epsilonList.append(e)
    print ("The epsilon list is:")
    print (epsilonList)

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
    #print (dicts_b)

    fromList =[]
    for i in range (len(dicts_b)-1):
        f = {}
        if dicts_b[i-1] == "{'from':" or dicts_b[i-1] == "[{'from':":
            f= int(dicts_b[i].replace("'," , " ").replace("," , " ").replace("}" , " ").replace("'" , " "))
            fromList.append(f)
    print ("The from list is:")
    print (fromList)

    toList =[]
    for i in range (len(dicts_b)-1):
        t = {}
        if dicts_b[i-1] == "'to':":
            t= int(dicts_b[i].replace("'," , " ").replace("," , " ").replace("}" , " ").replace("'" , " "))
            toList.append(t)
    print ("The to list is:")
    print (toList)

    lengthList =[]
    for i in range (len(dicts_b)-1):
        l = {}
        if dicts_b[i-1] == "'length':":
            l= float(dicts_b[i].replace("'," , " ").replace("," , " ").replace("}" , " ").replace("'" , " "))
            lengthList.append(l)
    print ("The length list is:")
    print (lengthList)

    k_bondList =[]
    for i in range (len(dicts_b)):
        l = {}
        if dicts_b[i-1] == "'k':":
            l= float(dicts_b[i].replace("'," , " ").replace("," , " ").replace("}" , " ").replace("'" , " ").replace("]" , " "))
            k_bondList.append(l)
    print ("The k_bond list is:")
    print (k_bondList)

    cl1List =[]
    for i in range (len(dicts_b)-1):
        l = {}
        if dicts_b[i-1] == "{'class1':"or dicts_b[i-1] == "[{'class1':":
            l= dicts_b[i].replace("'," , " ").replace("," , " ").replace("}" , " ").replace("'" , " ")
            cl1List.append(l)
    print ("The class1 list is:")
    print (cl1List)

    cl2List =[]
    for i in range (len(dicts_b)-1):
        l = {}
        if dicts_b[i-1] == "'class2':":
            l= dicts_b[i].replace("'," , " ").replace("," , " ").replace("}" , " ").replace("'" , " ")
            cl2List.append(l)
    print ("The class2 list is:")
    print (cl2List)   

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
    #print (dicts_ang)



    cl1_angle_List =[]
    for i in range ( len(dicts_ang)-1):
        l = {}
        if dicts_ang[i] == "{'class1':" or dicts_ang[i] == "[{'class1':":
            l= dicts_ang[i+1].replace("'," , " ").replace("," , " ").replace("}" , " ").replace("'" , " ")
            cl1_angle_List.append(l)
    print ("The class1 angle list is:")
    print (cl1_angle_List)



    dct = {}
    for k in range(2,4):
        dct["cl%s" %k] = []
        for i in range (len(dicts_ang)-1):
            l = {}
            if dicts_ang[i] == "'class%s"%k+"':":
                d={}
                l= dicts_ang[i+1].replace("'," , " ").replace("," , " ").replace("}" , " ").replace("'" , " ")
                dct["cl"+str(k)].append(l)
        print ("The class%s "%k+"angle list is:")
        print (dct["cl"+str(k)])


    angle_List =[]
    for i in range (len(dicts_ang)-1):
        l = {}
        if dicts_ang[i] == "'angle':":
            l= float(dicts_ang[i+1].replace("'," , " ").replace("," , " ").replace("}" , " ").replace("'" , " "))
            angle_List.append(l)
    print ("The angle list is:")
    print (angle_List)


    k_angle_List =[]
    for i in range (len(dicts_ang)-1):
        l = {}
        if dicts_ang[i] == "'k':":
            l= float(dicts_ang[i+1].replace("'}]" , " ").replace("'," , " ").replace("," , " ").replace("}" , " ").replace("'" , " "))
            k_angle_List.append(l)
    print ("The k angle list is:")
    print (k_angle_List)


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
    #print (dicts_pr)

    cl1_proper_List =[]
    for i in range ( len(dicts_pr)-1):
        l = {}
        if dicts_pr[i] == "{'class1':" or dicts_pr[i] == "[{'class1':":
            l= dicts_pr[i+1].replace("'," , " ").replace("," , " ").replace("}" , " ").replace("'" , " ")
            cl1_proper_List.append(l)
    print ("The class1 proper list is:")
    print (cl1_proper_List)

    dct = {}
    for k in range(2,5):
        dct["cl%s" %k] = []
        for i in range (len(dicts_pr)-1):
            l = {}
            if dicts_pr[i] == "'class%s"%k+"':":
                d={}
                l= dicts_pr[i+1].replace("'," , " ").replace("," , " ").replace("}" , " ").replace("'" , " ")
                dct["cl"+str(k)].append(l)
        print ("The class%s "%k+"proper list is:")
        print (dct["cl"+str(k)])

    dct = {}
    for k in range(1,5):
        dct["k%s" %k] = []
        for i in range (len(dicts_pr)-1):
            l = {}
            if dicts_pr[i] == "'k%s"%k+"':":
                d={}
                l= float(dicts_pr[i+1].replace("'," , " ").replace("," , " ").replace("}" , " ").replace("'" , " "))
                dct["k"+str(k)].append(l)
        print ("The k%s "%k+"proper list is:")
        print (dct["k"+str(k)])


    dct = {}
    for k in range(1,5):
        dct["periodicity%s" %k] = []
        for i in range (len(dicts_pr)-1):
            l = {}
            if dicts_pr[i] == "'periodicity%s"%k+"':":
                d={}
                l= int(dicts_pr[i+1].replace("'," , " ").replace("," , " ").replace("}" , " ").replace("'" , " "))
                dct["periodicity"+str(k)].append(l)
        print ("The periodicity%s "%k+"proper list is:")
        print (dct["periodicity"+str(k)])



    dct = {}
    for k in range(1,5):
        dct["phase%s" %k] = []
        for i in range (len(dicts_pr)-1):
            l = {}
            if dicts_pr[i] == "'phase%s"%k+"':":
                d={}
                l= float(dicts_pr[i+1].replace("'," , " ").replace("," , " ").replace("}" , " ").replace("'" , " ").replace("]", " "))
                dct["phase"+str(k)].append(l)
        print ("The phase%s "%k+"proper list is:")
        print (dct["phase"+str(k)])


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
    #print (dicts_impr)

    cl1_improper_List =[]
    for i in range ( len(dicts_impr)-1):
        l = {}
        if dicts_impr[i] == "{'class1':" or dicts_impr[i] == "[{'class1':":
            l= dicts_impr[i+1].replace("'," , " ").replace("," , " ").replace("}" , " ").replace("'" , " ")
            cl1_improper_List.append(l)
    print ("The class1 improper list is:")
    print (cl1_improper_List)

    dct = {}
    for k in range(2,5):
        dct["cl%s" %k] = []
        for i in range (len(dicts_impr)-1):
            l = {}
            if dicts_impr[i] == "'class%s"%k+"':":
                d={}
                l= dicts_impr[i+1].replace("'," , " ").replace("," , " ").replace("}" , " ").replace("'" , " ")
                dct["cl"+str(k)].append(l)
        print ("The class%s "%k+"improper list is:")
        print (dct["cl"+str(k)])


    dct = {}
    for k in range(1,5):
        dct["k%s" %k] = []
        for i in range (len(dicts_impr)-1):
            l = {}
            if dicts_impr[i] == "'k%s"%k+"':":
                d={}
                l= float(dicts_impr[i+1].replace("'," , " ").replace("," , " ").replace("}" , " ").replace("'" , " "))
                dct["k"+str(k)].append(l)
        print ("The k%s "%k+"improper list is:")
        print (dct["k"+str(k)])


    dct = {}
    for k in range(1,5):
        dct["periodicity%s" %k] = []
        for i in range (len(dicts_impr)-1):
            l = {}
            if dicts_impr[i] == "'periodicity%s"%k+"':":
                d={}
                l= int(dicts_impr[i+1].replace("'," , " ").replace("," , " ").replace("}" , " ").replace("'" , " "))
                dct["periodicity"+str(k)].append(l)
        print ("The periodicity%s "%k+"improper list is:")
        print (dct["periodicity"+str(k)])



    dct = {}
    for k in range(1,5):
        dct["phase%s" %k] = []
        for i in range (len(dicts_impr)-1):
            l = {}
            if dicts_impr[i] == "'phase%s"%k+"':":
                d={}
                l= float(dicts_impr[i+1].replace("'," , " ").replace("," , " ").replace("}" , " ").replace("'" , " ").replace("]", " "))
                dct["phase"+str(k)].append(l)
        print ("The phase%s "%k+"improper list is:")
        print (dct["phase"+str(k)])


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


    type_vs_List =[]
    for i in range (len(dicts_vs)-1):
        t = {}
        if dicts_vs[i-1] == "[{'type':":
            t= dicts_vs[i].replace("'," , " ").replace("," , " ").replace("}" , " ").replace("'" , " ")
            type_vs_List.append(t)
    print ("The type list for virtual sites is:")
    print (type_vs_List)


    index_vs_List =[]
    for i in range (len(dicts_vs)-1):
        t = {}
        if dicts_vs[i-1] == "'index':":
            t= dicts_vs[i].replace("'," , " ").replace("," , " ").replace("}" , " ").replace("'" , " ")
            index_vs_List.append(t)
    print ("The index list for virtual sites is:")
    print (index_vs_List)

    dct = {}
    for k in range(1,4):
        dct["atom%s" %k] = []
        for i in range (len(dicts_vs)-1):
            l = {}
            if dicts_vs[i] == "'atom%s"%k+"':":
                d={}
                l= int(dicts_vs[i+1].replace("'," , " ").replace("," , " ").replace("}" , " ").replace("'" , " "))
                dct["atom"+str(k)].append(l)
        print ("The atom%s "%k+" list for the virtual sites is:")
        print (dct["atom"+str(k)])

    dct = {}
    for k in range(1,3):
        dct["wo%s" %k] = []
        for i in range (len(dicts_vs)-1):
            l = {}
            if dicts_vs[i] == "'wo%s"%k+"':":
                d={}
                l= float(dicts_vs[i+1].replace("'," , " ").replace("," , " ").replace("}" , " ").replace("'" , " "))
                dct["wo"+str(k)].append(l)
        print ("The wo%s "%k+" list for the virtual sites is:")
        print (dct["wo"+str(k)])

    dct = {}
    for k in range(1,4):
        dct["wx%s" %k] = []
        for i in range (len(dicts_vs)-1):
            l = {}
            if dicts_vs[i] == "'wx%s"%k+"':":
                d={}
                l= float(dicts_vs[i+1].replace("'," , " ").replace("," , " ").replace("}" , " ").replace("'" , " "))
                dct["wx"+str(k)].append(l)
        print ("The wx%s "%k+" list for the virtual sites is:")
        print (dct["wx"+str(k)])


    dct = {}
    for k in range(1,4):
        dct["wy%s" %k] = []
        for i in range (len(dicts_vs)-1):
            l = {}
            if dicts_vs[i] == "'wy%s"%k+"':":
                d={}
                l= float(dicts_vs[i+1].replace("'," , " ").replace("," , " ").replace("}" , " ").replace("'" , " "))
                dct["wy"+str(k)].append(l)
        print ("The wy%s "%k+" list for the virtual sites is:")
        print (dct["wy"+str(k)])

    dct = {}
    for k in range(1,4):
        dct["p%s" %k] = []
        for i in range (len(dicts_vs)-1):
            l = {}
            if dicts_vs[i] == "'p%s"%k+"':":
                d={}
                l= dicts_vs[i+1].replace("'," , " ").replace("," , " ").replace("}" , " ").replace("'" , " ")
                dct["p"+str(k)].append(l)
        print ("The p%s "%k+" list for the virtual sites is:")
        print (dct["p"+str(k)])

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
#    print (dicts_res)


    residue_List =[]
    for i in range (len(dicts_res)):
        t = {}
        if dicts_res[i-1] == "[{'name':":
            t= dicts_res[i].replace("'," , " ").replace("," , " ").replace("}]" , " ").replace("'" , " ")
            residue_List.append(t)
    print ("The residue name is:")
    print (residue_List)

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
 #   print (dicts_nb)


    coulomb_List =[]
    for i in range (len(dicts_nb)):
        t = {}
        if dicts_nb[i-1] == "[{'coulomb14scale':":
            t= float(dicts_nb[i].replace("'," , " ").replace("," , " ").replace("}]" , " ").replace("'" , " "))
            coulomb_List.append(t)
    print ("The coulomb 1-4 scale is:")
    print (coulomb_List)


    LJ_List =[]
    for i in range (len(dicts_nb)):
        t = {}
        if dicts_nb[i-1] == "'lj14scale':":
            t= float(dicts_nb[i].replace("'," , " ").replace("," , " ").replace("}]" , " ").replace("'" , " "))
            LJ_List.append(t)
    print ("The Lennard Jones 1-4 scale is:")
    print (LJ_List)


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
