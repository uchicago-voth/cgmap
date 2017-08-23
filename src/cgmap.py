#!/usr/bin/env python2
import numpy as np
import pandas as pd
from mdtraj import Trajectory
from mdtraj import Topology
from mdtraj import element
from mdtraj.geometry import distance
from top_manip import typed_elementwise_rep

mapping_options = {}

def mode_rows(a):
    """Efficiently returns the most common row of a 2-D array."""
    #It is unclear how this method truly works.

    #generates a numpy view
    a = np.ascontiguousarray(a)

    #Modifies the basic data type of that view to be complete rows?
    void_dt = np.dtype((np.void, a.dtype.itemsize * np.prod(a.shape[1:])))

    #Counts the occurance of rows
    _,ids, count = np.unique(a.view(void_dt).ravel(), \
                                return_index=1,return_counts=1)
    largest_count_id = ids[count.argmax()]
    most_frequent_row = a[largest_count_id]
    return most_frequent_row

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

def get_forcenorms(traj,atom_indices=None,use_pbc=True):
    """Compute a list of time-averaged (over all frames) norms of forces for all atoms in a CG site
    Parameters
    ----------
    traj : Trajectory
        Trajectory to sum forces over
    atom_indices : array-like, dtype=int, shape=(n_atoms)
        List of indices for atoms in this CG site
    Returns
    -------
    forces : np.ndarray, shape=(n_atoms)
        Time-averaged norm of forces for each atom 
    """

    if atom_indices is not None and len(atom_indices) > 0 :
        forces = traj.forces[:,atom_indices,:]
    else:
        forces = traj.forces

    mean_forces = forces.astype('float64').mean(axis=0)
    mean_forcenorms = np.square(mean_forces.astype('float64')).sum(axis=1).flatten()

    return mean_forcenorms

def map_unfiltered_molecules(trj,selection_list,bead_label_list,transfer_labels=False,
                  molecule_types=None, molecule_type_order=False,
                  return_call=False,*args,**kwargs):
    """ This performs the mapping where each molecule and its beads has been assigned 
    in the parent script. This method is used when the original residue-based logic 
    for molecule assignment does not hold (e.g., proteins).

    Parameters
    ----------
    traj : Trajectory
        Trajectory to sum forces on
    selection_list :
        Indexible collection of strings
    bead_label_list :
        Indexible collection
    transfer_labels :
        Whether to transfer over labels in @trj. Moves over resSeq, resName
        for every bead, assuming that the atoms in each bead are uniform in
        those qualities.
    molecule_types :
        Indexible collection of integers
    molecule_type_order : boolean
        Specifying molecule_type_order means that the map will be
        reordered so that all molecules of type 0 come first, then 1, etc.
    return_call: boolean
        Whether to return the arguments that cg_by_index would be called with
        instead of actually calling it. Useful for modifying the call.

    Returns
    -------
    traj: trajectory
        trajectory formed by applying given molecular map.
    -OR-
    tuple: list of arguments which would be passed to cg_by_index
    """

    # The current features do not use molecule_type_order

    n_molecule_types = len(selection_list)

    if sorted(set(molecule_types)) != range(n_molecule_types):
        raise ValueError("Error in map molecules, molecule types list must "
                         "contain only and all numbers from 0 to "
                         "n_molecule_types-1.")

    if len(selection_list) != len(bead_label_list):
        raise ValueError("Error in map molecules, must submit selection list "
                         "and bead label list of same length.")

    for i in range(n_molecule_types):
        if len(selection_list[i]) != len(bead_label_list[i]):
            raise ValueError("Error in map molecules, selection list %i and "
                             "bead label list %i must be of same length."%(i,i))

    index_list = []
    label_list = []
    for mid in range(len(selection_list)):
        selMol = selection_list[mid]
        for bid in range(len(selMol)) :
            internal_indices = trj.top.select("%s"%(selMol[bid]))
            index_list.append(internal_indices)
            label_list.append(bead_label_list[mid][bid])

    cg_trj = cg_by_index(trj, index_list, label_list, *args, **kwargs)

    return(cg_trj)


def map_molecules(trj,selection_list,bead_label_list,transfer_labels=False,
                  molecule_types=None, molecule_type_order=False,
                  return_call=False,*args,**kwargs):
    """ This performs the mapping where each molecule has been assigned a
    type.

    Parameters
    ----------
    traj : Trajectory
        Trajectory to sum forces on
    selection_list :
        Indexible collection of strings
    bead_label_list :
        Indexible collection
    transfer_labels :
        Whether to transfer over labels in @trj. Moves over resSeq, resName
        for every bead, assuming that the atoms in each bead are uniform in
        those qualities.
    molecule_types :
        Indexible collection of integers
    molecule_type_order : boolean
        Specifying molecule_type_order means that the map will be
        reordered so that all molecules of type 0 come first, then 1, etc.
    return_call: boolean
        Whether to return the arguments that cg_by_index would be called with
        instead of actually calling it. Useful for modifying the call.

    Returns
    -------
    traj: trajectory
        trajectory formed by applying given molecular map.
    -OR-
    tuple: list of arguments which would be passed to cg_by_index
    """

    ### First, deal with optional arguments and argument validation.
    if molecule_type_order is True:
        raise ValueError("molecule_type_order not currently supported.")

    #if the array of molecule types isn't given, assume 1 molecule type.
    if molecule_types is None:
        molecule_types = [0]*trj.top.n_residues

    n_molecule_types = len(selection_list)

    if sorted(set(molecule_types)) != range(n_molecule_types):
        raise ValueError("Error in map molecules, molecule types list must "
                         "contain only and all numbers from 0 to "
                         "n_molecule_types-1.")

#    if len(molecule_types) != trj.top.n_residues:
#        raise ValueError("Error in map molecules, molecule types list must "
#                         "have the same length as number of residues.")

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
            internal_indices = []
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

    if (return_call is True):
        arg_list=[trj, index_list, label_list]
        arg_list.extend(args)
        arg_list.append(kwargs)
        return(arg_list)
        #exit early.

    cg_trj = cg_by_index(trj, index_list, label_list, *args, **kwargs)

    #do a more sophisticated labeling.
    if (transfer_labels is True):

        df_aa_top = trj.top.to_dataframe()[0]
        df_cg_top = cg_trj.top.to_dataframe()[0]

        #get resSeq info.
        aa_resSeq = df_aa_top.loc[:,'resSeq']

        #find atom indices for first atoms of each residue.
        res_starting_indices = \
            np.sort(np.unique(aa_resSeq,return_index=True)[1])

        #get resids and resnames for startings atoms.
        aa_starting_resids    = df_aa_top.loc[res_starting_indices,'resSeq']
        aa_starting_resnames  = df_aa_top.loc[res_starting_indices,'resName']

        #needed for duplicating atomistic info across cg molecules
        n_sites_per_cg = [ len(desc) for desc in bead_label_list ]

        #generate and place resids
        cg_resids = typed_elementwise_rep(aa_starting_resids,molecule_types,n_sites_per_cg)
        df_cg_top.loc[:,"resSeq"] = cg_resids

        #generate and place resNames
        cg_resnames = typed_elementwise_rep(aa_starting_resnames,molecule_types,n_sites_per_cg)
        df_cg_top.loc[:,"resName"] = cg_resnames

        #convert and put back.
        cg_trj.top = Topology.from_dataframe(df_cg_top)

    return(cg_trj)

def gen_unique_overlap_mod_weights(indices_list,assume_unique=True):
    """
    Parameters
    ----------
    indices_list: list of indices to look for modifications between.

    Returns
    -------
    weights: list
        list of weights
    """

    inverse_weights = [ np.ones(len(item)) for item in indices_list ]

    sets = [ set(item) for item in indices_list ]

    #Not very pythonic, but limited by triangle loop.
    for i1 in xrange(len(sets)):
        for i2 in xrange(i1+1,len(sets)):
            elements = np.array(list(sets[i1].intersection(sets[i2])))

            locs = np.in1d(indices_list[i1],elements,assume_unique=assume_unique)
            inverse_weights[i1][locs] = inverse_weights[i1][locs] + 1

            locs = np.in1d(indices_list[i2],elements,assume_unique=assume_unique)
            inverse_weights[i2][locs] = inverse_weights[i2][locs] + 1

    #we invert the weights before we return them.
    return([ np.divide(1.0,item) for item in inverse_weights])

def map_identical_molecules(trj,selection_list,bead_label_list,transfer_labels=False,
                            return_call=False,*args,**kwargs):
    """This performs the mapping assuming the entire system
    is a set of identical molecules. The argument are almost identical to map
    molecules, but do not have to be nested in lists.

    Parameters
    ----------
    traj : Trajectory
        Trajectory to sum forces on
    selection_list :
        Collection of strings
    bead_label_list :
        Collection of strings
    transfer_labels :
        Whether to transfer over labels in @trj. Moves over resSeq, resName
        for every bead, assuming that the atoms in each bead are uniform in
        those qualities.
    return_call: boolean
        Whether to return the arguments that cg_by_index would be called with
        instead of actually calling it. Useful for modifying the call.

    Returns
    -------
    traj: trajectory
        trajectory formed by applying given molecular map.
    -OR-
    tuple: list of arguments which would be passed to cg_by_index
    """

    null_molecule_types = list(np.zeros(trj.top.n_residues,dtype=np.int32))

    out = map_molecules(            trj = trj,
                          selection_list = [selection_list],
                         bead_label_list = [bead_label_list],
                         transfer_labels = transfer_labels,
                          molecule_types = null_molecule_types,
                     molecule_type_order = False,
                             return_call = return_call,
                                   *args,
                                   **kwargs)

    return(out)

def compute_center_weighted(xyz_i,xyz_all,atom_indices,weights=None,
                            unitcell_lengths=None,center_postwrap=False):
    """Compute the weighted center over selected atoms for a coordinate matrix.

    Parameters
    ----------
    xyz_i: array-like, dtype=float, shape=(n_steps,n_sites,n_dim)
            Holds output of linear map operations.
    xyz: array-like, dtype=float, shape=(n_steps,n_atoms,n_dim)
            Holds input coordinates.
    atom_indices : array-like, dtype=int, shape=(n_sites)
            List of indices of atoms to use in computing center
    weights : array-like, dtype=float, shape=(n_sites)
            Weights used to calculate positions (normalized in function)
    unitcell_lengths: array-like, dtype=float, shape=(n_dim)
            Unitcell lengths; Positions are calculated in a minimum image.
    center_postwrap: boolean
            Whether to wrap the derived CG coordinates around the given box
            after mapping. Assumes the box is centered at L/2, where L is a
            given periodic dimension.

    Returns
    -------
    None
    """

    xyz = xyz_all[:,atom_indices,:]

    if (weights is None):
        weights = np.ones(length(atom_indices))

    weights /= weights.sum()

    if unitcell_lengths is not None and center_postwrap is False:
        #performs a consensus shift on the resulting cg mapped site.
        for i, x in enumerate(xyz):
            # use periodic boundaries by centering relative to first
            # xyz coordinate, then shift back
            xyz0 = x[0,:]
            shift = unitcell_lengths[i]*np.floor( \
                    (x - xyz0)/unitcell_lengths[i] + 0.5)

            x = x - shift

            mode_shift = mode_rows(shift)

            mapped_points = x.astype('float64').T.dot(weights).flatten()

            xyz_i[i] = mapped_points + mode_shift
    elif unitcell_lengths is not None and center_postwrap is True:
        #performs post map wrapping around box center.
        for i, x in enumerate(xyz):
            # use periodic boundaries by centering relative to first
            # xyz coordinate, then shift back
            xyz0 = x[0,:]
            shift = unitcell_lengths[i]*np.floor( \
                    (x - xyz0)/unitcell_lengths[i] + 0.5)

            x = x - shift

            mapped_points = x.astype('float64').T.dot(weights).flatten()

            shift = unitcell_lengths[i]*np.floor( \
                    (mapped_points - unitcell_lengths[i]/2)/unitcell_lengths[i] + 0.5)

            xyz_i[i] = mapped_points - shift
    else:
        for i, x in enumerate(xyz):
            xyz_i[i] = x.astype('float64').T.dot(weights).flatten()

    return None

mapping_options['center'] = compute_center_weighted
mapping_options['com'] = compute_center_weighted
mapping_options['center_of_mass'] = compute_center_weighted
mapping_options['coc'] = compute_center_weighted
mapping_options['center_of_charge'] = compute_center_weighted
mapping_options['cof'] = compute_center_weighted
mapping_options['center_of_force'] = compute_center_weighted

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
                resSeq_list=None, inplace=False, bonds=None, split_shared_atoms=False, 
                mod_weights_list=None, mapping_function="com", charge_tol=1e-5,
                center_postwrap=False):
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
    split_shared_atoms: boolean
        If specified, check to see if atoms are shared per molecule in beads. If
        so, equally divide their weight accordingly for each bead.
    mapping_function: string, default='com': how to map xyz coordinates
        options: %s
    center_postwrap: Boolean
        Whether to wrap the CG system after it is mapped. Assumes that box is
        centered at 0, and only has effect if periodic information is present.

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

    #total masse for each cg bead.
    masses = np.zeros((n_beads),dtype=np.float64)
    #list of masses for elements in cg bead.
    masses_i = []
    #masses
    for ii in range(n_beads):
        #atoms in curent cg bead.
        atom_indices = atom_indices_list[ii]
        #first, construct lists of masses in current cg bead.
        temp_masses = np.array([])
        for jj in atom_indices:
            temp_masses = np.append(temp_masses, trj.top.atom(jj).element.mass)

        masses_i.append(temp_masses)
        masses[ii] = masses_i[ii].sum()

    if hasattr(trj.top.atom(1), 'charge'):
        #total charge for each cg bead.
        charges = np.zeros((n_beads), dtype=np.float64)
        #lists of charges for in current cg bead
        charges_i = []

        #charges
        for ii in range(n_beads):

            #atoms in curent cg bead.
            atom_indices = atom_indices_list[ii]

            #first, construct lists of masses in current cg bead.
            temp_charges = np.array([])

            for jj in atom_indices:
                temp_charges = np.append(temp_charges, trj.top.atom(jj).charge)

            charges_i.append(temp_charges)
            charges[ii] = charges_i[ii].sum()
    
    forcenorm_i = []
    if mapping_function == 'cof' or mapping_function == 'center_of_force' :
        for ii in range(n_beads) :
            atom_indices = atom_indices_list[ii]
            forcenorm_i.append(get_forcenorms(trj,atom_indices))

    if mapping_function == 'coc' or mapping_function == 'center_of_charge':
        for charge in charges:
            if np.absolute(charge) < charge_tol:
                raise ValueError("Total charge on site %i is near zero" % ii)

    topology_labels = []
    element_label_dict = {}

    if (split_shared_atoms):
        mod_weights_list = gen_unique_overlap_mod_weights(atom_indices_list) 

    has_forces = False
    try:
        trj.__dict__['forces']
        test_forces = map_forces(trj, (0,))
        has_forces = True
    except TypeError:
        print("WARNING: Invalid Forces\nNo Map applied to forces")
    except KeyError:
        pass
    except:
        print("Unknown error, check your forces\nexiting...")
        raise

    for i in range(n_beads):
        atom_indices = atom_indices_list[i]
        bead_label = bead_label_list[i]
        xyz_i = xyz[:,i,:]

        if mapping_function == 'coc' or mapping_function == 'center_of_charge':
            weights = charges_i[i]
        elif mapping_function == 'com' or mapping_function == 'center_of_mass':
            weights = masses_i[i]
        elif mapping_function == 'cof' or mapping_function == 'center_of_force':
            weights = forcenorm_i[i]
        elif mapping_function == 'center':
            weights = np.ones(len(atom_indices))

        if (mod_weights_list is not None):
            weights[:] = np.multiply(weights, mod_weights_list[i])

        compute_center_weighted(xyz_i,
                                trj.xyz,
                                atom_indices,
                                weights,
                                unitcell_lengths=trj.unitcell_lengths,
                                center_postwrap=center_postwrap)

        if has_forces:    
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
