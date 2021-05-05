import configparser
import DICOM
import MATLAB
import numpy as np
from pathlib import Path
import yaml


# extract data parameter object representing a single slice
def getSliceDataParams(dPar, slice, z):
    sliceDataParams = dict(dPar)
    sliceDataParams['sliceList'] = [slice]
    sliceDataParams['img'] = dPar['img'].reshape(
        dPar['N'], dPar['nz'], dPar['ny']*dPar['nx'])[:, z, :].flatten()
    sliceDataParams['nz'] = 1
    return sliceDataParams


# extract dPar object representing a slab of contiguous slices starting at z
def getSlabDataParams(dPar, slices, z):
    slabDataParams = dict(dPar)
    slabDataParams['sliceList'] = slices
    slabSize = len(slices)
    slabDataParams['img'] = dPar['img'].reshape(
        dPar['N'], dPar['nz'], dPar['ny']*dPar['nx'])[:, z:z+slabSize, :].flatten()
    slabDataParams['nz'] = slabSize
    return slabDataParams


# Update algorithm parameter object aPar and set default parameters
def setupAlgoParams(aPar, N, nFAC=0):
    defaults = [
        ('nR2', 1),
        ('R2max', 100.),
        ('R2cand', [0.]),
        ('mu', 1.),
        ('nB0', 100),
        ('nICMiter', 0),
        ('multiScale', False),
        ('use3D', False),
        ('magnitudeDiscrimination', True),
        ('offresPenalty', 0.)

    ]

    for param, defval in defaults:
        if param not in aPar:
            aPar[param] = defval

    if 'graphcut' not in aPar:
        aPar['graphcut'] = 'graphcutlevel' in aPar
    
    if aPar['graphcut']:
        if 'graphcutlevel' not in aPar:
            aPar['graphcutLevel'] = 0
    else:
        aPar['graphcutLevel'] = None

    if 'realEstimates' in aPar:
        if not aPar['realEstimates'] and N==2:
            raise Exception('Real-valued estimates needed for two-point Dixon')
    elif N==2:
        aPar['realEstimates'] = True
    else:
        aPar['realEstimates'] = False

    if aPar['nR2'] > 1:
        aPar['R2step'] = aPar['R2max']/(aPar['nR2']-1)  # [sec-1]
    else:
        aPar['R2step'] = 1.0  # [sec-1]
    
    aPar['iR2cand'] = np.array(list(set([min(aPar['nR2']-1, int(R2/aPar['R2step']))
                            for R2 in aPar['R2cand']])))  # [msec]

    aPar['maxICMupdate'] = round(aPar['nB0']/10)

    # For Fatty Acid Composition, create algorithmParams for two passes: aPar and aPar['pass2']
    # First pass: use standard fat-water separation to determine B0 and R2*
    # Second pass: use B0- and R2*-maps from first pass
    if nFAC > 0:
        aPar['pass2'] = dict(aPar)  # modify algoParams for pass 2:
        aPar['pass2']['nICMiter'] = 0  # to omit ICM
        aPar['pass2']['graphcutLevel'] = None  # to omit the graphcut
        aPar['pass2']['graphcut'] = False
    
    aPar['output'] = ['wat', 'fat', 'ff', 'B0map']
    if aPar['realEstimates']:
        aPar['output'].append('phi')
    if (aPar['nR2'] > 1):
        aPar['output'].append('R2map')
    if (nFAC > 2):
        aPar['output'].append('CL')
    if (nFAC > 1):
        aPar['output'].append('PUD')
    if (nFAC > 0):
        aPar['output'].append('UD')


# Get relative weights alpha of fat resonances based on CL, UD, and PUD per UD
def getFACalphas(CL=None, P2U=None, UD=None):
    P = 11  # Expects one water and ten triglyceride resonances
    M = [CL, UD, P2U].count(None)+2
    alpha = np.zeros([M, P], dtype=np.float32)
    alpha[0, 0] = 1.  # Water component
    if M == 2:
        # F = 9A+(6(CL-4)+UD(2P2U-8))B+6C+4UD(1-P2U)D+6E+2UDP2UF+2G+2H+I+2UDJ
        alpha[1, 1:] = [9, 6*(CL-4)+UD*(2*P2U-8), 6, 4*UD*(1-P2U), 6, 2*UD*P2U,
                        2, 2, 1, UD*2]
    elif M == 3:
        # F1 = 9A+6(CL-4)B+6C+6E+2G+2H+I
        # F2 = (2P2U-8)B+4(1-P2U)D+2P2UF+2J
        alpha[1, 1:] = [9, 6*(CL-4), 6, 0, 6, 0, 2, 2, 1, 0]
        alpha[2, 1:] = [0, 2*P2U-8, 0, 4*(1-P2U), 0, 2*P2U, 0, 0, 0, 2]
    elif M == 4:
        # F1 = 9A+6(CL-4)B+6C+6E+2G+2H+I
        # F2 = -8B+4D+2J
        # F3 = 2B-4D+2F
        alpha[1, 1:] = [9, 6*(CL-4), 6, 0, 6, 0, 2, 2, 1, 0]
        alpha[2, 1:] = [0, -8, 0, 4, 0, 0, 0, 0, 0, 2]
        alpha[3, 1:] = [0, 2, 0, -4, 0, 2, 0, 0, 0, 0]
    elif M == 5:
        # F1 = 9A-24B+6C+6E+2G+2H+I
        # F2 = -8B+4D+2J
        # F3 = 2B-4D+2F
        # F4 = 6B
        alpha[1, 1:] = [9, -24, 6, 0, 6, 0, 2, 2, 1, 0]
        alpha[2, 1:] = [0, -8, 0, 4, 0, 0, 0, 0, 0, 2]
        alpha[3, 1:] = [0, 2, 0, -4, 0, 2, 0, 0, 0, 0]
        alpha[4, 1:] = [0, 6, 0, 0, 0, 0, 0, 0, 0, 0]
    return alpha


# Update model parameter object mPar and set default parameters
def setupModelParams(mPar, clockwisePrecession=False, temperature=None):

    defaults = [
        ('fatCS', [1.3]),
        ('nFAC', 0),
        ('CL', 17.4), # Derived from Lundbom 2010
        ('P2U', 0.2),  # Derived from Lundbom 2010
        ('UD', 2.6),  # Derived from Lundbom 2010
    ]

    for param, defval in defaults:
        if param not in mPar:
            mPar[param] = defval

    if 'watCS' not in mPar:
        if temperature: # Temperature dependence according to Hernando 2014
            mPar['watCS'] = 1.3 + 3.748 -.01085 * temperature # Temp in [°C]
        else:
            mPar['watCS'] = 4.7
    
    mPar['CS'] = np.array([mPar['watCS']] + mPar['fatCS'], dtype=np.float32)
    
    if clockwisePrecession:
        mPar['CS'] *= -1
    
    mPar['P'] = len(mPar['CS'])

    if mPar['nFAC'] > 0 and mPar['P'] != 11:
        raise Exception(
            'FAC excpects exactly one water and ten triglyceride resonances')
    
    mPar['M'] = 2+mPar['nFAC']

    if mPar['nFAC'] == 0:
        mPar['alpha'] = np.zeros([mPar['M'], mPar['P']], dtype=np.float32)
        mPar['alpha'][0, 0] = 1.
        if 'relAmps' in mPar:
            for (p, a) in enumerate(mPar['relAmps']):
                mPar['alpha'][1, p+1] = float(a)
        else:
            for p in range(1, mPar['P']):
                mPar['alpha'][1, p] = float(1/len(fatCS))
    elif mPar['nFAC'] == 1:
        mPar['alpha'] = getFACalphas(mPar['CL'], mPar['P2U'])
    elif mPar['nFAC'] == 2:
        mPar['alpha'] = getFACalphas(mPar['CL'])
    elif mPar['nFAC'] == 3:
        mPar['alpha'] = getFACalphas()
    else:
        raise Exception('Unknown number of FAC parameters: {}'
                        .format(mPar['nFAC']))

    # For Fatty Acid Composition, create modelParams for two passes: mPar and mPar['pass2']
    # First pass: use standard fat-water separation to determine B0 and R2*
    # Second pass: do the Fatty Acid Composition
    if mPar['nFAC'] > 0: 
        mPar['pass2'] = dict(mPar) # copy mPar into pass 2, then modify pass 1
        mPar['alpha'] = getFACalphas(mPar['CL'], mPar['P2U'], mPar['UD'])
        mPar['M'] = mPar['alpha'].shape[0]


# group slices in sliceList in slabs of reconSlab contiguous slices
def getSlabs(sliceList, reconSlab):
    slabs = []
    slices = []
    pos = 0
    for z, slice in enumerate(sliceList):
        # start a new slab
        if slices and (len(slices) == reconSlab or not slice == slices[-1]+1):
            slabs.append((slices, pos))
            slices = [slice]
            pos = z
        else:
            slices.append(slice)
    slabs.append((slices, pos))
    return slabs

    
# Update data param object, set default parameters and read data from files
def setupDataParams(dPar, outDir=None):
    if outDir:
        dPar['outDir'] = Path(outDir)
    elif 'outDir' in dPar:
        dPar['outDir'] = Path(dPar['outDir'])
    else:
        raise Exception('No outDir defined')

    defaults = [
        ('reScale', 1.0),
        ('temperature', None),
        ('clockwisePrecession', False),
        ('offresCenter', 0.),
        ('files', [])
    ]

    for param, defval in defaults:
        if param not in dPar:
            dPar[param] = defval

    if 'files' in dPar:
        dPar['files'] = [dPar['configPath'] / file for file in list(dPar['files']) if Path(dPar['configPath'] / file).is_file()]
    
    if 'dirs' in dPar:
        dPar['dirs'] = [dPar['configPath'] / dir for dir in list(dPar['dirs']) if Path(dPar['configPath'] / dir).is_dir()]
        for path in dPar['dirs']:
            dPar['files'] += [obj for obj in path.iterdir() if obj.is_file()]
    
    validFiles = DICOM.getValidFiles(dPar['files'])
    
    if validFiles:
        DICOM.updateDataParams(dPar, validFiles)
    else:
        if len(dPar['files']) == 1 and dPar['files'][0].suffix == '.mat':
            MATLAB.updateDataParams(dPar, dPar['files'][0])
        else:
            raise Exception('No valid files found')
    
    if 'reconSlab' in dPar:
        dPar['slabs'] = getSlabs(dPar['sliceList'], dPar['reconSlab'])


# Read configuration file
def readConfig(file, section):
    file = Path(file)
    with open(file, 'r') as configFile:
        try:
            config = yaml.safe_load(configFile)
        except yaml.YAMLError as exc:
            raise Exception('Error reading config file {}'.format(file)) from exc
    config['configPath'] = file.parent
    return config