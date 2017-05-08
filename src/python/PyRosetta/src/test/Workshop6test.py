#! /usr/bin/python
# List of commands used in PyRosetts Workshop #6

from __future__ import print_function

from rosetta import *
from pyrosetta import *
from pyrosetta.toolbox import *
from pyrosetta.teaching import *

init(extra_options = "-constant_seed")  # WARNING: option '-constant_seed' is for testing only! MAKE SURE TO REMOVE IT IN PRODUCTION RUNS!!!!!
import os; os.chdir('.test.output')


# Side Chain Conformations, the Rotamer Library, and Dunbrack Energies
pose = pose_from_file("../test/data/workshops/1YY8.clean.pdb")
scorefxn = create_score_function("talaris2013")

for i in range(1, 3):
    pose.set_chi(i, 49, 180)

# Monte Carlo Side-Chain Packing
task_pack = standard_packer_task(pose)
task_pack.restrict_to_repacking()
task_pack.temporarily_fix_everything()
task_pack.temporarily_set_pack_residue(49, True)

print( task_pack )

pack_mover = PackRotamersMover(scorefxn, task_pack)

pack_mover.apply(pose)

import os, tempfile

# Design

# work around for windows permission problem
YY8_resfile1 = tempfile.mkstemp()[1]
YY8_resfile2 = tempfile.mkstemp()[1]

generate_resfile_from_pdb("../test/data/workshops/1YY8.clean.pdb", YY8_resfile1)
generate_resfile_from_pose(pose, YY8_resfile2)

task_design = core.pack.task.TaskFactory.create_packer_task(pose)
core.pack.task.parse_resfile(pose, task_design, YY8_resfile2)

mutate_residue(pose, 49, 'E')

#os.remove(YY8_resfile1)
#os.remove(YY8_resfile2)
