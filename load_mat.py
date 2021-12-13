#Code from https://stackoverflow.com/questions/7008608/scipy-io-loadmat-nested-structures-i-e-dictionaries
#Improve scipy.io.loadmat : Really transform .mat in dictionary

#Correct error with string array
import scipy.io as spio
import numpy as np


def loadmat(filename):
    '''
    this function should be called instead of direct spio.loadmat
    as it cures the problem of not properly recovering python dictionaries
    from mat files. It calls the function check keys to cure all entries
    which are still mat-objects
    '''
    def _check_keys(d):
        '''
        checks if entries in dictionary are mat-objects. If yes
        todict is called to change them to nested dictionaries
        '''
        for key in d:
            if isinstance(d[key],np.ndarray):
                if isinstance(d[key][0,0], spio.matlab.mio5_params.mat_struct):
                    d[key] = _todict(d[key][0,0])
        return d

    def _has_struct(elem):
        """Determine if elem is an array and if any array item is a struct"""
        struct_flag=False
        if isinstance(elem, np.ndarray) and elem.ndim>0:
            if any(isinstance(e, spio.matlab.mio5_params.mat_struct) for e in elem):
                struct_flag= True
        return struct_flag
        
        #return isinstance(elem, np.ndarray) and any(isinstance(
        #            e, spio.matlab.mio5_params.mat_struct) for e in elem)

    def _todict(matobj):
        '''
        A recursive function which constructs from matobjects nested dictionaries
        '''
        d = {}
        for strg in matobj._fieldnames:
            elem = matobj.__dict__[strg]
            if elem.shape[0]==1:
                elem=elem[0,...]
            if isinstance(elem, spio.matlab.mio5_params.mat_struct):
                d[strg] = _todict(elem)
            elif _has_struct(elem):
                d[strg] = _tolist(elem)
            else:
                d[strg] = elem
        return d

    def _tolist(ndarray):
        '''
        A recursive function which constructs lists from cellarrays
        (which are loaded as numpy ndarrays), recursing into the elements
        if they contain matobjects.
        '''
        elem_list = []
        for sub_elem in ndarray:
            if isinstance(sub_elem, spio.matlab.mio5_params.mat_struct):
                elem_list.append(_todict(sub_elem))
            elif _has_struct(sub_elem):
                elem_list.append(_tolist(sub_elem))
            else:
                elem_list.append(sub_elem)
        return elem_list
    data = spio.loadmat(filename, struct_as_record=False, squeeze_me=False)
    return _check_keys(data)