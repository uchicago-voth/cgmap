#!/usr/bin/env python2

import sys
sys.path.append('../../../src/')
import cgmap as cg
import mdtraj as md
import md_check as check

############################### config #####################################

input_traj = "dppc_dopc.trr"
input_top  = "dppc_dopc.pdb"

input_maps_dppc = ["mapping_bead_1_dppc",
                   "mapping_bead_2_dppc",
                   "mapping_bead_3_dppc"]

input_maps_dopc = ["mapping_bead_1_dopc",
                   "mapping_bead_2_dopc",
                   "mapping_bead_3_dopc"]

input_maps_system = [input_maps_dppc, input_maps_dopc]

output_traj = "dppc_dopc.trr"
output_top  = "dppc_dopc.pdb"

reference_traj = "dppc_dopc.trr"
reference_top  = "dppc_dopc.pdb"

output_dir ='./output/'
input_dir  ='./input/'
reference_dir ='./reference/'

#collection of names of molecules in the fine-grained trajectory.
lipid_types = ['DPPC','DOPC']

#names of cg beads created in cg trajectory
label_lists = [['DPH','DPM','DPT'],['DOH','DOM','DOT']]

############################### config proc ################################

fq_input_maps = [[ input_dir + loc for loc in input_maps ] 
                                   for input_maps in input_maps_system ]

#creates a map (list of strings) from a filename
def read_map(filename):
    with open(filename,'r') as mapping_file:
        symbols = [ l.strip() for l in mapping_file.readlines() ]
    return symbols

mapping_atom_names_system = []
for mfiles in fq_input_maps:
    mapping_atom_names_system.append([read_map(g) for g in mfiles])

#creates a atom selection (string) from list of atom names
def map_to_atom_sel(mapping):
    return " or ".join( ["name %s"%mn for mn in mapping])

#index strings for which to atom query the trajectory when making beads.
#list of lists of strings.
name_lists = [[ map_to_atom_sel(mnames) for mnames in mapping_list] 
                                        for mapping_list in mapping_atom_names_system ]

############################### run ########################################

### pull in trajectories
trj = md.load(input_dir + input_traj,top=input_dir + input_top)

#the types of each molecule in the trajectory.
molecule_types = [lipid_types.index(r.name) for r in trj.top.residues]

#actual map command
cg_trj = cg.map_molecules(            trj = trj,
                           selection_list = name_lists, 
                          bead_label_list = label_lists, 
                           molecule_types = molecule_types,
                          transfer_labels = True)

cg_trj.save(output_dir + output_traj) 
cg_trj[0].save(output_dir + output_top)

############################### check results ###############################
# reloading results from disk.

cg_traj = cg_trj.load(output_dir + output_traj,top=output_dir + output_top) 
ref_cg_traj = cg_trj.load(reference_dir + reference_traj,
                          top=reference_dir + reference_top) 

result=check.md_content_equality(cg_traj,ref_cg_traj)

sys.exit(check.check_result_to_exitval(result))
