#python3

'''Contains functions used for manipulating mdtraj topologies, 
including functions which act on numpy arrays.'''

import numpy as np
import mdtraj as md
import pandas as pd

def typed_elementwise_rep(values,types,type_repeat):
    '''Repeats values in an array a set number of times, based
    on a given type.

    Parameters
    ----------
    values: numpy array (1D)
        Contains the values that will be repeated.
    types: numpy array (1D)
        Determines the type of each element in @values.
    type_repeat:
        Given the value in @types, contains the number of times
        to repeat that element in @values

    Returns
    -------
    rep_array:
        Contains repeated values.

    Notes
    -----
    - No check on validity between arguments is done.
    - FAILS IF @values IS OF MUTABLE OBJECTS.

    Examples
    --------

    values = [ 5, 12, 7 ]
    types  = [ 1, 0,  0 ]
    type_repeat = [ 2, 4 ]

    -> [ 5, 5, 12, 12, 12, 12, 7 ,7 ]
    '''

    #assuming all arguments are valid.

    rep_array = []

    for value,vtype in zip(values,types):
        rep_array.extend([value]*type_repeat[vtype])

    return(np.array(rep_array))
