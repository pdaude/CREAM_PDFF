#!/usr/bin/env python3

import main
import numpy as np
import pydicom
import scipy.io
import time
from pathlib import Path


def getScore(case, dir, refFile):
    # Read reconstructed MATLAB-file
    file = dir / '0.mat'
    try:
        mat = scipy.io.loadmat(file)
    except:
        raise Exception('Could not read MATLAB file {}'.format(file))
    recFF = mat['ff'].flatten(order='F')/100
    recFF.shape = recFF.shape + (1,)
    # Read reference MATLAB-file
    try:
        mat = scipy.io.loadmat(refFile)
    except:
        raise Exception('Could not read MATLAB file {}'.format(refFile))
    refFF = mat['REFCASES'][0, case - 1]
    mask = mat['MASKS'][0, case - 1]
    # Calculate score
    score = 100 * \
        (1 - np.sum((np.abs(refFF - recFF) > 0.1) * mask) / np.sum(mask))
    return score


if __name__ == '__main__':

    fwqpboPath = Path(__file__).resolve().parent.absolute()
    challengePath = fwqpboPath / './challenge'

    refFile = challengePath / './refdata.mat'
    if not refFile.is_file():
        url = 'https://challenge.ismrm.org/u/5131732/ISMRM_Challenge/refdata.mat'
        raise Exception('ISMRM 2012 challenge reference data file was not found at: {}. '
                        'Please download from: {}.'.format(refFile.absolute(), url))

    modelParamsFile = fwqpboPath / './modelParams.yml'
    cases = [1, 2, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17]

    for case in cases:
        caseFile = challengePath / './{}.mat'.format(str(case).zfill(2))
        if not caseFile.is_file():
            url = 'https://challenge.ismrm.org/u/5131732/ISMRM_Challenge/data_matlab.zip'
            raise Exception('ISMRM 2012 challenge dataset {} was not found at: {}. '
                            'Please download from: {}.'.format(case, caseFile.absolute(), url))

    results = []
    for case in cases:
        dataParamsFile = challengePath / './{}.yml'.format(case)
        if case == 9:
            algoParamsFile = fwqpboPath / './algoParams2D.yml'
        else:
            algoParamsFile = fwqpboPath / './algoParams3D.yml'
        outDir = challengePath / './{}_REC'.format(str(case).zfill(2))
        t = time.time()
        main.main(dataParamsFile, algoParamsFile, modelParamsFile, outDir)
        results.append((case, getScore(case, outDir, refFile), time.time() - t))

    print()
    for case, score, recTime in results:
        print('Case {}: score {}% in {} sec'.format(
            case, round(score, 2), round(recTime, 1)))
