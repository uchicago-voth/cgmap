#!/usr/bin/env python2

import sys
sys.path.append('../../../src/')
import cgmap as cg
import mdtraj as md
import md_check as check

############################### config #####################################

input_traj = "protein.trr"
input_top  = "protein.pdb"

output_traj = "protein.trr"
output_top  = "protein.pdb"

reference_traj = "protein.trr"
reference_top  = "protein.pdb"

output_dir ='./output/'
input_dir  ='./input/'
reference_dir ='./reference/'


############################### run ########################################

### pull in trajectories
trj = md.load(input_dir + input_traj,top=input_dir + input_top)

### define mapping based on knowledge of topology
### in this instance, map every residue into a single site
for a in trj.top.atoms: a.mass = a.element.mass
for a in trj.top.atoms: a.charge = 0

# first residue is SER148 (zero index'd)
name_lists = []
label_lists = []
molecule_types = []
resREF = 148
istart = 0
iend = 0
iname = "SER"
molnum = 0


maxSize = len(list(trj.top.atoms))
stopFlag = False
tempMol = []
tempCGL = []
name_lists_key = []
for i, a in enumerate(trj.top.atoms) :
    resNAME = str(a.residue)[0:3]
    resNUM = int(str(a.residue)[3:6])
    aINDEX = a.index
        
    if resNAME not in name_lists_key :
        name_lists_key.append(resNAME)

    if (resNUM != resREF) :
        #first append name_lists and label
        iend = aINDEX - 1
        tempMol.append("index %d to %d" % (istart, iend))
        tempCGL.append(iname)

        #then update things for next residue
        iname = resNAME
        istart = aINDEX
        if resNUM < resREF :
            #stopFlag = True
            molecule_types.append(int(molnum))
            name_lists.append(tempMol)
            label_lists.append(tempCGL)
            tempMol = []
            tempCGL = []
            molnum += 1
        resREF = resNUM

    # special case if last item
    if (i == (maxSize-1)) :
        iend = aINDEX
        tempMol.append("index %d to %d" % (istart, iend))
        tempCGL.append(iname)
        molecule_types.append(int(molnum))
        name_lists.append(tempMol)
        label_lists.append(tempCGL)

#actual map command
print name_lists
print label_lists
print molecule_types

print "Lengths of all three lists should be equivalent: %d = %d = %d" % (len(name_lists), len(label_lists), len(molecule_types))

cg_trj = cg.map_unfiltered_molecules(            trj = trj,
                                      selection_list = name_lists, 
                                     bead_label_list = label_lists,
                                      molecule_types = molecule_types,
                                    mapping_function = "com")

cg_trj.save(output_dir + output_traj) 
cg_trj[0].save(output_dir + output_top)

############################### check results ###############################
# reloading results from disk.

cg_traj = cg_trj.load(output_dir + output_traj,top=output_dir + output_top) 
ref_cg_traj = cg_trj.load(reference_dir + reference_traj,
                          top=reference_dir + reference_top) 

result=check.md_content_equality(cg_traj,ref_cg_traj)

sys.exit(check.check_result_to_exitval(result))
