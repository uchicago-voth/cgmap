#!/usr/bin/env python2

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

output_traj = "dppc.trr"
output_top  = "dppc.pdb"

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

############################### run ########################################

### pull in trajectories
trj = md.load(input_dir + input_traj,top=input_dir + input_top)

#the types of each molecule in the trajectory.
molecule_types = [lipid_types.index(r.name) for r in trj.top.residues]

#actual map command
cg_trj = cg.map_molecules(            trj = trj,
                           selection_list = [ name_lists  ], 
                          bead_label_list = [ label_lists ], 
                           molecule_types = molecule_types)

cg_trj.save(output_dir + output_traj) 
cg_trj[0].save(output_dir + output_top)

############################### check results ###############################
# reloading results from disk.

cg_traj = cg_trj.load(output_dir + output_traj,top=output_dir + output_top) 
ref_cg_traj = cg_trj.load(reference_dir + reference_traj,
                          top=reference_dir + reference_top) 

result=check.md_content_equality(cg_traj,ref_cg_traj)
sys.exit(check.check_result_to_exitval(result))
