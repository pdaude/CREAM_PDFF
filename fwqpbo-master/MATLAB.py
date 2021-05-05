import scipy.io
import numpy as np
from pathlib import Path


# update dPar with information retrieved from MATLAB file
# (arranged according to ISMRM fat-water toolbox)
def updateDataParams(dPar, file):
    dPar['fileType'] = 'MATLAB'
    try:
        mat = scipy.io.loadmat(file)
    except:
        raise Exception('Could not read MATLAB file {}'.format(file))
    data = mat['imDataParams'][0, 0]

    for i in range(0, 4):
        if len(data[i].shape) == 5:
            img = data[i]  # Image data (row,col,slice,coil,echo)
        elif data[i].shape[1] > 2:
            echoTimes = data[i][0]  # TEs [sec]
        else:
            if data[i][0, 0] > 1:
                dPar['B0'] = data[i][0, 0]  # Fieldstrength [T]
            else:
                clockwise = data[i][0, 0]  # Clockwiseprecession?

    if clockwise != 1:
        raise Exception('Warning: Not clockwise precession. ' +
                        'Need to write code to handle this case!')

    dPar['ny'], dPar['nx'], dPar['nz'], nCoils, dPar['N'] = img.shape
    if nCoils > 1:
        raise Exception('Warning: more than one coil. ' +
                        'Need to write code to coil combine!')

    # Get only slices in dPar['sliceList']
    if 'sliceList' not in dPar:
        dPar['sliceList'] = range(dPar['nz'])
    else:
        img = img[:, :, dPar['sliceList'], :, :]
        dPar['nz'] = len(dPar['sliceList'])
    # Get only echoes in dPar['echoes']
    dPar['totalN'] = dPar['N']
    if not 'echoes' in dPar:
        dPar['echoes'] = range(dPar['totalN'])
    else:
        img = img[:, :, :, :, dPar['echoes']]
        echoTimes = echoTimes[dPar['echoes']]
        dPar['N'] = len(dPar['echoes'])
    if 'cropFOV' in dPar:
        x0, x1 = dPar['cropFOV'][0], dPar['cropFOV'][1]
        y0, y1 = dPar['cropFOV'][2], dPar['cropFOV'][3]
        dPar['Nx'], dPar['nx'] = dPar['nx'], x1-x0
        dPar['Ny'], dPar['ny'] = dPar['ny'], y1-y0
        img = img[y0:y1, x0:x1, :, :, :]
    if dPar['N'] < 2:
        raise Exception(
            'At least 2 echoes required, only {} given'.format(dPar['N']))
    dPar['t1'] = echoTimes[0]
    dPar['dt'] = np.mean(np.diff(echoTimes))
    if np.max(np.diff(echoTimes))/dPar['dt'] > 1.05 or np.min(
      np.diff(echoTimes))/dPar['dt'] < .95:
        raise Exception('Warning: echo inter-spacing varies more than 5%')

    dPar['frameList'] = []

    dPar['dx'], dPar['dy'], dPar['dz'] = 1.5, 1.5, 5  # Ad hoc assumption on voxelsize

    # To get data as: (echo,slice,row,col)
    img.shape = (dPar['ny'], dPar['nx'], dPar['nz'], dPar['N'])
    img = np.transpose(img)
    img = np.swapaxes(img, 2, 3)

    img = img.flatten()
    dPar['img'] = img*dPar['reScale']


# Save output as MATLAB arrays
def save(output, dPar):
    dPar['outDir'].mkdir(parents=True, exist_ok=True)
    filename = dPar['outDir'] / './{}.mat'.format(dPar['sliceList'][0])
    print(r'Writing images to "{}"'.format(filename))
    scipy.io.savemat(filename, output)
