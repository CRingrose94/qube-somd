For a summary of the tests

python3 driver.py | grep '@@@@'

Otherwise inspect the full standard output for more detailed logs

TODO NEXT
- modify setupMovesFreeEnergy to define how to pass information about combining rules to use to OpenMM
- update interface of ::OpenMMFrEnergyST to store combining rule property
- update OpenMMFrEnergyST::initialise() to support use of arithmetic or geometric combining rules 
- verify that agreement between Sire and OpenMM energies with geometric rules is similar or better to that seen for arithmetic rules

