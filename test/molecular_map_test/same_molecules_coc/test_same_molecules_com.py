#!/usr/bin/env python2

import sys
sys.path.append('../../../src/')
import cgmap as cg
import mdtraj as md
import md_check as check

############################### config #####################################

input_traj = "acetonitrile.lammpstrj"
input_top  = "acetonitrile.pdb"

mapping_indices = ["index 0 or index 5",
                         "index 1 to 4"]

output_traj = "acetonitrile.lammpstrj"
output_top  = "acetonitrile.pdb"

reference_traj = "acetonitrile.lammpstrj"
reference_top  = "acetonitrile.pdb"

output_dir ='./output/'
input_dir  ='./input/'
reference_dir ='./reference/'

############################### config proc ################################

#names of cg beads created.
label_lists = ['neg', 'pos']

############################### run ########################################

### pull in trajectories
trj = md.load(input_dir + input_traj,top=input_dir + input_top)

#preprocess trajectory content by adding new masses and charges 
for a in trj.top.atoms: 
    a.mass = a.element.mass
    if a.name == 'H':
        a.charge = 0.060
    elif a.name == 'CT':
        a.charge = -0.080
    elif a.name == 'CZ':
        a.charge = 0.460
    elif a.name == 'NZ':
        a.charge = -0.560

#actual map command
cg_trj = cg.map_identical_molecules(trj, 
                                       mapping_indices, 
                                       label_lists, 
                                       mapping_function='coc')

cg_trj.save(output_dir + output_traj) 
cg_trj[0].save(output_dir + output_top)

############################### check results ###############################
# reloading results from disk.

cg_traj = cg_trj.load(output_dir + output_traj,top=output_dir + output_top) 
ref_cg_traj = cg_trj.load(reference_dir + reference_traj,
                          top=reference_dir + reference_top) 

result=check.md_content_equality(cg_traj,ref_cg_traj)
sys.exit(check.check_result_to_exitval(result))
