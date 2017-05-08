#!/usr/bin/python

# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

## @file   /GUIs/tests/Test_Module_Region.py
## @brief  PyRosetta Toolkit Test script
## @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

from __future__ import print_function

#Python Imports
import os
import sys

#Must have Toolkit added to path!
#Toolkit Imports
from app.pyrosetta_toolkit.modules.PythonPDB import PythonPDB

print "Testing General Python PDB Class"
filename = os.path.dirname(os.path.abspath(__file__))+"/data/gui/sdnNB2Pmin.pdb"
python_pdb = PythonPDB(filename)

print "Testing Cleaning a PDB"
python_pdb.clean_PDB()

print "Testing Changing Occupancy"
python_pdb.change_occupancy()

print "Testing replacing B factors with arbitrary data"
python_pdb.read_file_and_replace_b_factors(",", os.path.dirname(os.path.abspath(__file__))+"/data/gui/SASAchangeForCSV.csv", 3, 2, 5)

print "Testing removing alternate residues"
python_pdb.remove_alternate_residues()

print "Testing removing an arbitrary chain"
python_pdb.remove_chain("C")

print "Testing dumping the PDB: "
python_pdb.save_PDB(".test.output/gui_python_pdb_test.pdb")
