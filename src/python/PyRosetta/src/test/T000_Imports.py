# :noTabs=true:
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

## @author Sergey Lyskov

from __future__ import print_function

import pyrosetta
from pyrosetta import *


print( Vector1([1, 2, 3, 4]) )
print( Vector1([1, 2, 3, 4.5]) )
print( Vector1(['a', 'b', 'c', 'd']) )


print('Testing Python bindings for Rosetta...')
print('Init...')
pyrosetta.init()


print( pyrosetta.rosetta_version )
print( version() )
