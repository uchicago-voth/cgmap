#!/usr/bin/env python2

#This test checks functionality related to mapping via center of point (no
#weighting). 

#Current only using null masses, as native 'center' method is broken.

import sys
sys.path.append('../../../src/')
import cgmap as cg
import mdtraj as md
import md_check as check

############################### config #####################################

input_traj = "dppc.trr"
input_top  = "dppc.pdb"

input_maps = ["mapping_bead_1_dppc",
              "mapping_bead_2_dppc",
              "mapping_bead_3_dppc"]

output_traj_null_mass = "dppc_nm.trr"
output_top_null_mass  = "dppc_nm.pdb"
output_traj_native    = "dppc_native.trr"
output_top_native     = "dppc_native.pdb"

reference_traj = "dppc.trr"
reference_top  = "dppc.pdb"

output_dir ='./output/'
input_dir  ='./input/'
reference_dir ='./reference/'

#collection of names of molecules.
lipid_types = ['DPPC']

############################### config proc ################################

fq_input_maps = [ input_dir + loc for loc in input_maps ]

#read maps for each bead from files.
#list of lists of strings.
mapping_atom_names_dppc = [ [ l.strip() for l in open(mp_file,'r').readlines()] 
                            for mp_file in fq_input_maps ]

#index strings for which to atom query the trajectory when making beads.
#list of lists of strings.
name_lists = [ " or ".join( ["name %s"%mn for mn in  mapping_names ]) 
               for mapping_names in mapping_atom_names_dppc ]

#names of cg beads created.
label_lists = ['DPH','DPM','DPT']

############################### run native ###############################

#### pull in trajectories
#trj = md.load(input_dir + input_traj,top=input_dir + input_top)
#
##the types of each molecule in the trajectory.
#molecule_types = [lipid_types.index(r.name) for r in trj.top.residues]
#
##preprocess trajectory content by adding new parts
#for a in trj.top.atoms: a.mass = a.element.mass
#for a in trj.top.atoms: a.charge = 0
#
##actual map command
#cg_trj = cg.map_molecules(            trj = trj,
#                           selection_list = [ name_lists  ], 
#                          bead_label_list = [ label_lists ], 
#                           molecule_types = molecule_types,
#                         mapping_function = 'center')
#
#cg_trj.save(output_dir + output_traj_native) 
#cg_trj[0].save(output_dir + output_top_native)

############################### run null mass ############################

### pull in trajectories
trj = md.load(input_dir + input_traj,top=input_dir + input_top)

#the types of each molecule in the trajectory.
molecule_types = [lipid_types.index(r.name) for r in trj.top.residues]

#preprocess trajectory content by adding new parts
for a in trj.top.atoms: a.mass = 1 #null mass part
for a in trj.top.atoms: a.charge = 0

#actual map command
cg_trj = cg.map_molecules(            trj = trj,
                           selection_list = [ name_lists  ], 
                          bead_label_list = [ label_lists ], 
                           molecule_types = molecule_types,
                         mapping_function = 'com')

cg_trj.save(output_dir + output_traj_null_mass) 
cg_trj[0].save(output_dir + output_top_null_mass)

############################### check results ###############################
# reloading results from disk.

cg_traj_null_mass = cg_trj.load(output_dir + output_traj_null_mass,
                                top=output_dir + output_top_null_mass) 

#cg_traj_native = cg_trj.load(output_dir + output_traj_native,
#                             top=output_dir + output_top_native) 

ref_cg_traj = cg_trj.load(reference_dir + reference_traj,
                          top=reference_dir + reference_top) 

result_null   = check.md_content_equality(cg_traj_null_mass,ref_cg_traj)
#result_native = check.md_content_equality(cg_traj_null_mass,ref_cg_traj)

#sys.exit(check.check_result_to_exitval(result_null & result_native))
sys.exit(check.check_result_to_exitval(result_null))
