import os,sys

top_files = ['inputs/toluene~methane/vacuum/input/SYSTEM.top',
             'inputs/2-cyclopentanylindole~7-cyclopentanylindole/vacuum/input/SYSTEM.top',
             'inputs/2-methylindole~methane/vacuum/input/SYSTEM.top',
             'inputs/methane~methanol/vacuum/input/SYSTEM.top',
             'inputs/methanol~ethane/vacuum/input/SYSTEM.top',
             'inputs/neopentane~methane/vacuum/input/SYSTEM.top',
             'inputs/toluene~methane/vacuum/input/SYSTEM.top'
             ]
crd_files = ['inputs/toluene~methane/vacuum/input/SYSTEM.crd',
             'inputs/2-cyclopentanylindole~7-cyclopentanylindole/vacuum/input/SYSTEM.crd',
             'inputs/2-methylindole~methane/vacuum/input/SYSTEM.crd',
             'inputs/methane~methanol/vacuum/input/SYSTEM.crd',
             'inputs/methanol~ethane/vacuum/input/SYSTEM.crd',
             'inputs/neopentane~methane/vacuum/input/SYSTEM.crd',
             'inputs/toluene~methane/vacuum/input/SYSTEM.crd' ]

pert_files = ["inputs/toluene~methane/vacuum/input/MORPH.pert",
             'inputs/2-cyclopentanylindole~7-cyclopentanylindole/vacuum/input/MORPH.pert',
             'inputs/2-methylindole~methane/vacuum/input/MORPH.pert',
             'inputs/methane~methanol/vacuum/input/MORPH.pert',
             'inputs/methanol~ethane/vacuum/input/MORPH.pert',
             'inputs/neopentane~methane/vacuum/input/MORPH.pert',
             'inputs/toluene~methane/vacuum/input/MORPH.pert' ]


# Compare single point energies between Sire and SOMD for test set using arithmetic combining rules

print ("@@@@ TESTING ARITHMETIC COMBINING RULES")

for x in range(0,len(top_files)):
    TOP = top_files[x]
    CRD = crd_files[x]
    PRT = pert_files[x]
    print ("@@@@ Single point energies for %s " % TOP)
    sys.stdout.flush()
    cmd = "~/sire.app/bin/python singlepoint.py %s %s %s %s " % (TOP, CRD, PRT, "arithmetic")
    os.system(cmd)

print ("@@@@ TESTING GEOMETRIC COMBINING RULES")

# Compare single point energies between Sire and SOMD for test set using geometric combining rules
for x in range(0,len(top_files)):
    TOP = top_files[x]
    CRD = crd_files[x]
    PRT = pert_files[x]
    print ("@@@@ Single point energies for %s " % TOP)
    sys.stdout.flush()
    cmd = "~/sire.app/bin/python singlepoint.py %s %s %s %s " % (TOP, CRD, PRT, "geometric")
    os.system(cmd)

