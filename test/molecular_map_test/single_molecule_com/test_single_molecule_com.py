#!/usr/bin/env python2

import sys
sys.path.append('../../../src/')
import cgmap as cg
import mdtraj as md
import md_check as check

input_traj = "dppc_1_molecule.trr"
input_top  = "dppc_1_molecule.pdb"

output_traj = "dppc_1_molecule.trr"
output_top  = "dppc_1_molecule.pdb"

reference_traj = "dppc_1_molecule.trr"
reference_top  = "dppc_1_molecule.pdb"

output_dir ='./output/'
input_dir  ='./input/'
reference_dir ='./reference/'

### create new trajectories

trj = md.load(input_dir + input_traj,top=input_dir + input_top)

#actual map command
cg_trj = cg.cg_by_selection(trj,['all'],['C1'])

cg_trj.save(output_dir + output_traj) 
cg_trj[0].save(output_dir + output_top)

### check results
# reloading from scratch.

cg_traj = cg_trj.load(output_dir + output_traj,top=output_dir + output_top) 
ref_cg_traj = cg_trj.load(reference_dir + reference_traj,
        top=reference_dir + reference_top) 

result=check.md_content_equality(cg_traj,ref_cg_traj)
sys.exit(check.check_result_to_exitval(result))
