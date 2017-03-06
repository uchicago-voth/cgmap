#!/usr/bin/env python2
import numpy as np
import pandas as pd
from mdtraj import Trajectory
from mdtraj import Topology
from mdtraj import element
from mdtraj.geometry import distance

mapping_options = {}

def map_forces(traj,atom_indices=None,use_pbc=True):
    """Compute the center of mass for each frame.
    Parameters
    ----------
    traj : Trajectory
        Trajectory to sum forces on
    atom_indices : array-like, dtype=int, shape=(n_atoms)
            List of indices of atoms to use in computing com
    Returns
    -------
    forces : np.ndarray, shape=(n_frames, 3)
         Summed up forces on atom indices
    """

    if atom_indices is not None and len(atom_indices)>0:
        forces = traj.forces[:,atom_indices,:]
    else:
        forces = traj.forces

    mapped_forces = forces.sum(axis=1)

    return mapped_forces

def map_molecules(trj,selection_list,bead_label_list,molecule_types=None,
                  molecule_type_order=False,*args,**kwargs):
    """ This performs the mapping where each molecule has been assigned a 
    type. 
        
    Parameters
    ----------
    traj : Trajectory
        Trajectory to sum forces on
    selection_list : indexible collection of strings
    bead_label_list : indexible collection
    molecule_types : indexible collection of integers
    molecule_type_order : boolean
        Specifying molecule_type_order means that the map will be 
        reordered so that all molecules of type 0 come first, then 1, etc.
   
    -------
    traj: trajectory
        trajectory formed by applying given molecular map.
    """

    ### First, deal with optional arguments and argument validation.

    #if the array of molecule types isn't given, assume 1 molecule type.
    if molecule_types is None:
        molecule_types = [0]*trj.top.n_residues

    n_molecule_types = len(selection_list)

    if sorted(set(molecule_types)) != range(n_molecule_types):
        raise ValueError("Error in map molecules, molecule types list must "
                         "contain only and all numbers from 0 to "
                         "n_molecule_types-1.")

    if len(molecule_types) != trj.top.n_residues:
        raise ValueError("Error in map molecules, molecule types list must "
                         "have the same length as number of residues.")

    if len(selection_list) != len(bead_label_list):
        raise ValueError("Error in map molecules, must submit selection list "
                         "and bead label list of same length.")

    for i in range(n_molecule_types):
        if len(selection_list[i]) != len(bead_label_list[i]):
            raise ValueError("Error in map molecules, selection list %i and "
                             "bead label list %i must be of same length."%(i,i))

    ### generate the indices local to each molecule for mapping

    # get the first molecule index for each molecule type
    first_molecules = [molecule_types.index(i) for i in range(n_molecule_types)]

    internal_indices_list = [[] for i in range(n_molecule_types)]

    iterable = zip(selection_list,first_molecules,internal_indices_list)
    for selection,first_mol,mol_indices in iterable:
        first_index = trj.top.select("(resid == %i)"%(first_mol)).min()

        for sel in selection:
            has_index = sel.find("index")> -1
            has_name  = sel.find("name") > -1

            if has_index and has_name:
                raise ValueError("Error in map molecules, do not specify "
                                 "selection by index and by type.")
            elif has_index:
                # use atom selection language to parse selection 
                #string containing only indices on whole system, then offset later
                internal_indices = trj.top.select("%s"%(sel))

            elif has_name:
                # have to un-shift list because this will be added to current id later
                filter_string = "(resid == %i) and (%s)"%(first_mol,sel)
                internal_indices = trj.top.select(filter_string) - first_index

            if len(internal_indices)==0:
                raise ValueError("Error in map_molecules, selection string '%s'"
                                 "produced an empty list of atom indices"%sel)
            mol_indices.append(internal_indices)


    # get list of type [ (0,r0), (1,r1) etc ]
    if molecule_type_order is True:
        residue_list = sorted( enumerate(trj.top.residues),\
                               key=lambda x: molecule_types[x[0]])
    else:
        residue_list = enumerate(trj.top.residues)

    index_list  = []
    resSeq_list = []
    label_list  = []

    start_index = 0
    resSeq = 1
    for ridx,r in residue_list:
        molecule_type = molecule_types[ridx]
        for bead_idx, internal_indices in enumerate(internal_indices_list[molecule_type]):
            system_indices = internal_indices + start_index
            index_list.append(system_indices) 
            resSeq_list.append(resSeq)
            label_list.append(bead_label_list[molecule_type][bead_idx])
        resSeq = resSeq+1
        start_index = start_index + r.n_atoms

    return cg_by_index(trj, index_list, label_list, *args, **kwargs)

def map_identical_molecules(trj,selection_list,bead_label_list,*args,**kwargs):
    """ This performs the mapping assuming the entire system 
    is a set of identical molecules"""
    index_list = []
    resSeq_list = []
    label_list = []

    internal_indices_list = [] 

    if len(selection_list) != len(bead_label_list):
        raise ValueError("Error in map_identical_molecules, must submit "
                         "selection list and bead label list of same length")

    for sel in selection_list:
        internal_indices = trj.top.select("(resSeq == 1) and (%s)"%sel)
        if len(internal_indices)==0:
            raise ValueError("Error in map_identical_molecules, selection "
                             "string '%s' produced an empty list of "
                             "atom indices"%sel)
        internal_indices_list.append(internal_indices)

    start_index = 0
    resSeq = 1
    for r in trj.top.residues:
        for bead_idx, internal_indices in enumerate(internal_indices_list):
            system_indices = internal_indices + start_index
            index_list.append(system_indices) 
            resSeq_list.append(resSeq)
            label_list.append(bead_label_list[bead_idx])
        resSeq = resSeq+1
        start_index = start_index + r.n_atoms

    return cg_by_index(trj, index_list, label_list, *args, **kwargs)

def compute_center(traj,atom_indices=None,use_pbc=True):
    """Compute the center of mass for each frame.
    Parameters
    ----------
    traj : Trajectory
        Trajectory to compute center of mass for
    atom_indices : array-like, dtype=int, shape=(n_atoms)
            List of indices of atoms to use in computing center
    Returns
    -------
    center : np.ndarray, shape=(n_frames, 3)
         Coordinates of the mean position of atom_indices for each frame
    """

    if atom_indices is not None and len(atom_indices)>0:
        xyz = traj.xyz[:,atom_indices,:]
    else:
        xyz = traj.xyz

    center = np.zeros((traj.n_frames, 3))

    for i, x in enumerate(xyz):
# use periodic boundaries by centering relative to first xyz coordinate, then shift back
        if use_pbc is True:
            xyz0 = x[0,:]
            shift = traj[i].unitcell_lengths*np.floor( (x - xyz0)/traj[i].unitcell_lengths + 0.5) 
            x = x - shift
        center[i, :] = x.astype('float64').mean(axis=0)
        
    return center
mapping_options['center'] = compute_center

def compute_com(traj,atom_indices=None,use_pbc=True):
    """Compute the center of mass for each frame.
    Parameters
    ----------
    traj : Trajectory
        Trajectory to compute center of mass for
    atom_indices : array-like, dtype=int, shape=(n_atoms)
            List of indices of atoms to use in computing com
    Returns
    -------
    com : np.ndarray, shape=(n_frames, 3)
         Coordinates of the center of mass for each frame
    """

    if atom_indices is not None and len(atom_indices)>0:
        xyz = traj.xyz[:,atom_indices,:]
        masses = np.array([a.mass for a in traj.top.atoms if a.index in atom_indices])
    else:
        xyz = traj.xyz
        masses = np.array([a.mass for a in traj.top.atoms])

    com = np.zeros((traj.n_frames, 3))
    masses /= masses.sum()

    for i, x in enumerate(xyz):
    # use periodic boundaries by centering relative to 
    #first xyz coordinate, then shift back
        if use_pbc is True:
            xyz0 = x[0,:]
            shift = traj[i].unitcell_lengths*np.floor( \
                    (x - xyz0)/traj[i].unitcell_lengths + 0.5) 
            x = x - shift
        com[i, :] = x.astype('float64').T.dot(masses).flatten()
        
    return com
mapping_options['com_slow'] = compute_com
mapping_options['center_of_mass_slow'] = compute_com

def compute_com_fast(xyz_i,xyz_all,atom_indices,masses,unitcell_lengths=None):
    """Compute the center of mass for each frame.
    Parameters
    ----------
    traj : Trajectory
        Trajectory to compute center of mass for
    atom_indices : array-like, dtype=int, shape=(n_atoms)
            List of indices of atoms to use in computing com
    Returns
    -------
    com : np.ndarray, shape=(n_frames, 3)
         Coordinates of the center of mass for each frame
    """

    #com = xyz_i
    masses /= masses.sum()
    xyz = xyz_all[:,atom_indices,:]


    for i, x in enumerate(xyz):
    # use periodic boundaries by centering relative to first 
    # xyz coordinate, then shift back
        if unitcell_lengths is not None:
            xyz0 = x[0,:]
            shift = unitcell_lengths[i]*np.floor( \
                    (x - xyz0)/unitcell_lengths[i] + 0.5) 
            x = x - shift
        xyz_i[i] = x.astype('float64').T.dot(masses).flatten()
        
    return
#    return com
mapping_options['com'] = compute_com_fast
mapping_options['center_of_mass'] = compute_com_fast

def cg_by_selection(trj, selection_string_list, *args, **kwargs):
    """Create a coarse grained (CG) trajectory from list of 
    atom selections by computing centers of mass of selected 
    sets of atoms.

    Parameters
    ----------
    selection_string_list : 
        list of strings in mdtraj selection 
        language, shape=(n_beads,)
    bead_label_list : 
        list of maximum 4-letter strings to 
        label CG sites
    chain_list : 
        optional list of chain id's to split resulting 
        beads into separate chains
    resSeq_list :     
        optional list of residue sequence id's to 
        assign cg residues
    segment_id_list : 
        optional list of segment id's to assign 
        cg residues
    inplace : bool, default=False
        If ``True``, the operation is done inplace, modifying ``trj``.
        Otherwise, a copy is returned with the sliced atoms, and
        ``trj`` is not modified.
    bonds : array-like,dtype=int, shape=(n_bonds,2), default=None
        If specified, sets these bonds in new topology 
    Returns
    -------
    traj : md.Trajectory
        The return value is either ``trj``, or the new trajectory,
        depending on the value of ``inplace``.
    """

    atom_indices_list = []
    for sel_string in selection_string_list:
        atom_indices_list.append( trj.top.select(sel_string) )
        if len(atom_indices_list[-1])==0:
            print("Warning - selection string returns 0 atoms: '%s'"%sel_string)
    return cg_by_index(trj, atom_indices_list, *args, **kwargs )

def cg_by_index(trj, atom_indices_list, bead_label_list, chain_list=None, segment_id_list=None, 
                resSeq_list=None, inplace=False, bonds=None, mapping_function="com"):
    """Create a coarse grained (CG) trajectory from subsets of atoms by 
        computing centers of mass of selected sets of atoms.
    Parameters
    ----------
    atom_indices_list : 
        list of array-like, dtype=int, shape=(n_beads,n_atoms)
        List of indices of atoms to combine into CG sites
    bead_label_list : 
        list of maximum 4-letter strings to label CG sites
    chain_list : 
        optional list of chain id's to split resulting beads into separate 
        chains
    resSeq_list : 
        optional list of residue sequence id's to assign cg residues
    segment_id_list : 
        optional list of segment id's to assign cg residues
    inplace : 
        bool, default=False
        If ``True``, the operation is done inplace, modifying ``trj``.
        Otherwise, a copy is returned with the sliced atoms, and
        ``trj`` is not modified.
    bonds : array-like,dtype=int, shape=(n_bonds,2), default=None
        If specified, sets these bonds in new topology 
    mapping_function: string, default='com': how to map xyz coordinates
        options: %s

    Note - If repeated resSeq values are used, as for a repeated motiff 
        in a CG polymer, those sections most be broken into separate 
        chains or an incorrect topology will result
 
    Returns
    -------
    traj : md.Trajectory
        The return value is either ``trj``, or the new trajectory,
        depending on the value of ``inplace``.
    """%mapping_options.keys()

    if not len(atom_indices_list)==len(bead_label_list):
        raise ValueError("Must supply a list of bead labels of the "
                         "same length as a list of selected atom indices")
    for bead_label in bead_label_list:
        if not (type(bead_label) is str) or len(bead_label)>4 or len(bead_label)<1:
            raise ValueError("Specified bead label '%s' is not valid, \
                             must be a string between 1 and 4 characters"%bead_label)
    bead_label_list = [ bead_label.upper() for bead_label in bead_label_list ]

    if mapping_function not in mapping_options:
        raise ValueError("Must select a mapping function from: %s"\
                         %mapping_options.keys())

    map_coords = mapping_options[mapping_function]

    if chain_list is None:
        chain_list = np.ones(len(atom_indices_list),dtype=int)
    elif len(chain_list)!=len(atom_indices_list):
        raise ValueError("Supplied chain_list must be of the same length "
                         "as a list of selected atom indices")

    if segment_id_list is not None and len(segment_id_list)!=len(atom_indices_list):
        raise ValueError("Supplied segment_id_list must be of the same "
                         "length as a list of selected atom indices")

    if resSeq_list is not None and len(resSeq_list)!=len(atom_indices_list):
        raise ValueError("Supplied resSeq_list must be of the same "
                         "length as a list of selected atom indices")

    n_beads = len(atom_indices_list)

    xyz = np.zeros((trj.xyz.shape[0], n_beads,trj.xyz.shape[2]),
                   dtype=trj.xyz.dtype,
                   order='C')

    forces = np.zeros((trj.xyz.shape[0],n_beads,trj.xyz.shape[2]),
                      dtype=np.double,
                      order='C')
    columns = ["serial","name","element","resSeq","resName","chainID"]

    masses = np.zeros((n_beads),dtype=np.float64)
    masses_i = []
    charges = np.zeros((n_beads), dtype=np.float64)
    for ii in range(n_beads):
        atom_indices = atom_indices_list[ii]
        temp = np.array([])
        for jj in atom_indices:
            temp = np.append(temp, trj.top.atom(jj).mass)
            charges[ii] += trj.topology.atom(jj).charge
        masses_i.append(temp)
        masses[ii] = masses_i[ii].sum()

    topology_labels = []
    element_label_dict = {}

    xyz_i = np.zeros((trj.xyz.shape[0],trj.xyz.shape[2]),
                     dtype=trj.xyz.dtype,
                     order='C')

    for i in range(n_beads):
        atom_indices = atom_indices_list[i]
        bead_label = bead_label_list[i]
        #xyz_i = map_coords(trj,atom_indices)

        map_coords(xyz_i,
                   trj.xyz,
                   atom_indices,
                   masses_i[i],
                   unitcell_lengths=trj.unitcell_lengths)

        xyz[:,i,:] = xyz_i

        if "forces" in trj.__dict__ and len(trj.forces)>0:
            forces_i = map_forces(trj,atom_indices)
            forces[:,i,:] = forces_i

        if resSeq_list is not None:
            resSeq = resSeq_list[i]
        else:
            resSeq = i + 1 

        #element_label='%4s'%('B%i'%(resSeq))
        if not bead_label in element_label_dict:
            element_label='%2s'%('B%i'%(len(element_label_dict)%10))
            element_label_dict[bead_label] = element_label
        else:
            element_label = element_label_dict[bead_label]

        if element_label.strip().upper() not in element.Element._elements_by_symbol:
            element.Element(1000+resSeq, 
                            element_label, 
                            element_label, 
                            masses[i], 
                            1.0)

        topology_labels.append( [i,
                                 bead_label,
                                 element_label,
                                 resSeq,
                                 '%3s'%bead_label,
                                 chain_list[i]] )

    df = pd.DataFrame(topology_labels,columns=columns)
    topology = Topology.from_dataframe(df,bonds=bonds)
    
    if segment_id_list is not None:
        for beadidx,bead in enumerate(topology.atoms):
            bead.residue.segment_id = segment_id_list[beadidx]
        
    if inplace:
        if trj._topology is not None:
            trj._topology = topology
        trj._xyz = xyz

        return trj

    unitcell_lengths = unitcell_angles = None
    if trj._have_unitcell:
        unitcell_lengths = trj._unitcell_lengths.copy()
        unitcell_angles = trj._unitcell_angles.copy()

    time = trj._time.copy()

    new_trj = Trajectory(xyz=xyz, topology=topology, time=time,
                      unitcell_lengths=unitcell_lengths,
                      unitcell_angles=unitcell_angles)

    new_trj.forces = forces
    return new_trj
