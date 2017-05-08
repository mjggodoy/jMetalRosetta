# :noTabs=true:
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

## @author Sergey Lyskov

from __future__ import print_function

from rosetta import *
from pyrosetta import *

init(extra_options = "-constant_seed")  # WARNING: option '-constant_seed' is for testing only! MAKE SURE TO REMOVE IT IN PRODUCTION RUNS!!!!!
import os; os.chdir('.test.output')


print('Scoring ----------------------------------------------')
pose = core.import_pose.pose_from_file("../test/data/test_fragments.pdb")

print('Creating standard fullatom score function and scoring')
scorefxn = create_score_function('talaris2013')
scorefxn(pose)

print('Creating standard centroid score function and scoring')
scorefxn = create_score_function('score3')
scorefxn(pose)

print('Creating standard score function and scoring, again')
scorefxn = create_score_function('talaris2013', 'docking')
scorefxn(pose)
print('Creating standard score from scratch')
scorefxn = ScoreFunction()

print('Adjusting weights and scoring')
scorefxn.set_weight(core.scoring.fa_atr, 1)
scorefxn.set_weight(core.scoring.fa_pair, 1)
scorefxn.set_weight(core.scoring.fa_rep, 1)
scorefxn(pose)

print('weight for fa_atr set to: ', scorefxn.get_weight(core.scoring.fa_atr))
print(scorefxn)
print('Score break down for pose...')
scorefxn.show(pose)

print( 'all residue energies are' )
print( pose.energies().show() )
print( 'weighted energies for residue 5' )
print( pose.energies().show(5) )
#TODO: align headers in .show statements above
#TODO: add TOTAL residue energy column (weighted)

#old way:
weights = pose.energies().weights()
print( pose.energies().residue_total_energies(5).weighted_string_of( weights ) )
print( 'fa_atr of residue 5' )
print( weights[core.scoring.fa_atr] * pose.energies().residue_total_energies(5)[core.scoring.fa_atr ] )
print( "fa_atr for residue 5: ", pose.energies().residue_total_energies(5)[core.scoring.fa_atr] )

print('manually calculating 2body context-independent energies between residues 4 and 5')
rsd1 = pose.residue(4);
rsd2 = pose.residue(5);

emap =  core.scoring.EMapVector()  # TwoBodyEMapVector()
scorefxn.eval_ci_2b( rsd1, rsd2, pose, emap );

print( "fa_atr between 1 and 2: ", emap[core.scoring.fa_atr] )
print( "fa_pair between 1 and 2: ", emap[core.scoring.fa_pair] )
print( "fa_rep between 1 and 2: ", emap[core.scoring.fa_rep] )

score_types = []
for i in range(1, int(rosetta.core.scoring.end_of_score_type_enumeration) + 1):
    ii = rosetta.core.scoring.ScoreType(i)
    if weights[ii] != 0: score_types.append(ii)

for i in range(1, int(rosetta.core.scoring.end_of_score_type_enumeration) + 1):
    print( rosetta.core.scoring.ScoreType(i), end='')   # This print only number, is this correct???  (no!)
#print
#TAB-COMPLETION gives this:
#rosetta.core.scoring.ScoreType.[tab]
# TODO-lowpriority: instead of using tab-completion, accomplish the above with:
# print rosetta.core.scoring.ScoreType

print()
print('identifying hydrogen bonds in structure')

pose = core.import_pose.pose_from_file("../test/data/test_in.pdb")
scorefxn = get_fa_scorefxn() #  create_score_function('standard')
scorefxn(pose)

hbond_set = rosetta.core.scoring.hbonds.HBondSet()
pose.update_residue_neighbors();
rosetta.core.scoring.hbonds.fill_hbond_set( pose, False, hbond_set )
hbond_set.show(pose)
# TODO: can we add an hbondset lookup by residue number?
