import numpy as np

def md_content_equality(traj_1,traj_2,prefix="Traj equality: ",xyz_abs_tol=1e-8):

    result=True

    if (not(traj_1.top == traj_2.top)):
        print(prefix+"Topologies don't match.")
        result=False


    if (not((traj_1.xyz == traj_2.xyz).all())):
        print(prefix+"Warning: Coordinates don't match bit for bit.")

        sqdiff_mat = (traj_1.xyz-traj_2.xyz)**2
        min_mat    = np.minimum(abs(traj_1.xyz),abs(traj_2.xyz))

        residual       = sqdiff_mat.mean()**0.5
        max_residual   = sqdiff_mat.max()**0.5
        print(prefix+"Warning: RMS coordinate residual: {}".format(residual))
        print(prefix+"Warning: Max coordinate residual: {}".format(max_residual))

        rel_residual_mat = sqdiff_mat**0.5/min_mat

        residual         = rel_residual_mat.mean()
        max_residual     = rel_residual_mat.max()
        print(prefix+"Warning: Mean relative coordinate residual: {}".format(residual))
        print(prefix+"Warning: Max relative coordinate residual: {}".format(max_residual))

        print(prefix+"First sqdiff coordinate frame: {}".format(\
                (traj_1.xyz[:,:,1]-traj_2.xyz[:,:,1])))

        violations = np.nonzero((traj_1.xyz-traj_2.xyz) > xyz_abs_tol)

        print(prefix+"Indices violating residual ({}): {}".format(\
                xyz_abs_tol,violations))

        if (residual > xyz_abs_tol):
            result=False

    if (not(traj_1.n_atoms == traj_2.n_atoms)):
        print(prefix+"Number of atoms doesn't match.")
        print(prefix+"Traj 1 n_atoms: {}; traj 2 n_atoms\
                {}".format(traj_1.n_atoms,traj_2.n_atoms))
        result=False

    if (not(traj_1.n_frames == traj_2.n_frames )):
        print(prefix+"Number of frames doesn't match.")
        print(prefix+"Traj 1 n_frames: {}; traj 2 n_frames\
                {}".format(traj_1.n_frames,traj_2.n_frames))
        result=False

    if (not(traj_1.n_residues == traj_2.n_residues)):
        print(prefix+"Number of frames doesn't match.")
        print(prefix+"Traj 1 n_residues: {}; traj 2 n_residues\
                {}".format(traj_1.n_residues,traj_2.n_residues))
        result=False

    if (not(traj_1.timestep  == traj_2.timestep)):
        print(prefix+"Timestep value doesn't match.")
        print(prefix+"Traj 1 timestep: {}; traj 2 timestep\
                {}".format(traj_1.timestep,traj_2.timestep))
        result=False

    if (not(traj_1.unitcell_vectors == traj_2.unitcell_vectors).all()):
        print(prefix+"Unit cell vectors don't match.")
        print(prefix+"Traj 1 unitcell_vectors: {}; traj 2 unitcell_vectors\
                {}".format(traj_1.unitcell_vectors,traj_2.unitcell_vectors))
        result=False

    if (not(traj_1.unitcell_lengths == traj_2.unitcell_lengths).all()):
        print(prefix+"Unit cell lengths don't match.")
        print(prefix+"Traj 1 unitcell_lengths: {}; traj 2 unitcell_lengths\
                {}".format(traj_1.unitcell_lengths,traj_2.unitcell_lengths))
        result=False

    if (not(traj_1.unitcell_angles == traj_2.unitcell_angles).all()):
        print(prefix+"Unit cell angles don't match.")
        print(prefix+"Traj 1 unitcell_angles: {}; traj 2 unitcell_angles\
                {}".format(traj_1.unitcell_angles,traj_2.unitcell_angles))
        result=False

    return result

def check_result_to_exitval(result):
    '''Transforms boolean to command line exit value. 
    True -> 0, False -> 1. No guard logic.
    '''

    return int(not(result))
