#!/usr/bin/env python2

#This test checks functionality related to periodic boundary conditions. A
#trajectory is wrapped, then mapped, and then compared to a trajectory which is
#wrapped, then mapped. An additional mapped then wrapped case is considered.

#This test has higher numerical inconsistency than the rest. It is unclear why.
#Visual inspection shows good agreement.

import sys
sys.path.append('../../../src/')
import cgmap as cg
import mdtraj as md
import md_check as check

############################### config #####################################

xyz_tolerance = 1e-6

input_traj              = "dppc.trr"
input_top               = "dppc.pdb"
input_traj_wrapped      = "dppc_vmd_wrap.trr"
input_traj_wrapped_join = "dppc_vmd_wrap_join.trr"

input_maps = ["mapping_bead_1_dppc",
              "mapping_bead_2_dppc",
              "mapping_bead_3_dppc"]

output_traj              = "dppc.trr"
output_traj_wrapped      = "dppc_vmd_wrapped.trr"
output_traj_wrapped_join = "dppc_vmd_wrapped_join.trr"
output_top               = "dppc.pdb"

reference_traj = "dppc.trr"
reference_top  = "dppc.pdb"

output_dir    = './output/'
input_dir     = './input/'
reference_dir = './reference/'

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

############################### run ######################################

### pull in trajectories
trj              = md.load(input_dir + input_traj,              top=input_dir + input_top)
trj_wrapped      = md.load(input_dir + input_traj_wrapped,      top=input_dir + input_top)
trj_wrapped_join = md.load(input_dir + input_traj_wrapped_join, top=input_dir + input_top)

#the types of each molecule in the trajectory.
molecule_types = [lipid_types.index(r.name) for r in trj.top.residues]

#actual map command
cg_trj = cg.map_molecules(            trj = trj,
                           selection_list = [ name_lists  ], 
                          bead_label_list = [ label_lists ], 
                           molecule_types = molecule_types,
                         mapping_function = 'center',
                          center_postwrap = True)

cg_trj_wrapped = cg.map_molecules(    trj = trj_wrapped,
                           selection_list = [ name_lists  ], 
                          bead_label_list = [ label_lists ], 
                           molecule_types = molecule_types,
                         mapping_function = 'center',
                          center_postwrap = True)

cg_trj_wrapped_join = cg.map_molecules(trj = trj_wrapped_join,
                            selection_list = [ name_lists  ], 
                           bead_label_list = [ label_lists ], 
                            molecule_types = molecule_types,
                          mapping_function = 'center',
                           center_postwrap = True)

cg_trj.save(            output_dir + output_traj) 
cg_trj_wrapped.save(    output_dir + output_traj_wrapped) 
cg_trj_wrapped_join.save(output_dir + output_traj_wrapped_join) 
cg_trj[0].save(output_dir + output_top)

############################### check results ###############################
print("!!! WARNING !!!: Test performed with an absolute xyz tolerance of {}; there exist "
      "higher levels of errors here for unknown reasons.".format(xyz_tolerance))
print("!!! WARNING !!!: Look at the actual statistics, and local test output, "
      "to understand these effects.")

cg_trj               = cg_trj.load(output_dir + output_traj,             top=output_dir + output_top) 
cg_trj_wrapped       = cg_trj.load(output_dir + output_traj_wrapped,     top=output_dir + output_top) 
cg_trj_wrapped_join  = cg_trj.load(output_dir + output_traj_wrapped_join,top=output_dir + output_top) 

ref_cg_trj = cg_trj.load(reference_dir + reference_traj,
                          top=reference_dir + reference_top) 

print("Checking completely unwrapped trajectory...")
result_base    = check.md_content_equality(cg_trj,
                                            ref_cg_trj,xyz_abs_tol=xyz_tolerance)
print("Checking completely wrapped trajectory...")
result_wrapped = check.md_content_equality(cg_trj_wrapped,
                                            ref_cg_trj,xyz_abs_tol=xyz_tolerance)
print("Checking sanely unwrapped trajectory...")
result_wrapped_join = check.md_content_equality(cg_trj_wrapped_join,
                                            ref_cg_trj,xyz_abs_tol=xyz_tolerance)

sys.exit(check.check_result_to_exitval(result_base & result_wrapped & result_wrapped_join))
