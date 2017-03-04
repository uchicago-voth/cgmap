def md_content_equality(traj_1,traj_2):

    if (not(traj_1.top == traj_2.top)):
        return False

    if (not((traj_1.xyz == traj_2.xyz).all())):
        return False

    if (not(traj_1.n_atoms == traj_2.n_atoms)):
        return False

    if (not(traj_1.n_frames == traj_2.n_frames )):
        return False

    if (not(traj_1.n_residues == traj_2.n_residues)):
        return False

    if (not(traj_1.timestep  == traj_2.timestep)):
        return False

    if (not(traj_1.unitcell_vectors == traj_2.unitcell_vectors).all()):
        return False

    if (not(traj_1.unitcell_lengths == traj_2.unitcell_lengths).all()):
        return False

    if (not(traj_1.unitcell_angles == traj_2.unitcell_angles).all()):
        return False

    return True

def check_result_to_exitval(result):
    '''Transforms boolean to command line exit value. 
    True -> 0, False -> 1. No guard logic.
    '''

    return int(not(result))
