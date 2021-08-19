import scipy.io
import numpy as np
from pathlib import Path
import load_mat as lmat
import os
import os.path as op

# update dPar with information retrieved from MATLAB file
# (arranged according to ISMRM fat-water toolbox)
def updateDataParams(dPar, file):
    """

    :param dPar:
    :param file:
    List of params add in dPar:
    B0 : Fieldstrength [T]
    ny,nx,nz,totalN : x,y,z, of complex img
    sliceList : range of nz
    echoes: index of echo time
    N : number of echoes
    Nx,Ny : if cropFOV
    t1 : first echo time
    dt : echo inter spacing
    img : cmplx image (echo,slice,row,col) =(N,nz,ny,nx)
    :return:
    """
    
    dPar['fileType'] = 'MATLAB'
    try:
        mat = lmat.loadmat(file)
    except:
        raise Exception('Could not read MATLAB file {}'.format(file))


    #Get variable names in imDataParams :
    #Mandatory variables
    # images : (nx,ny,nz,ncoil,nt) complex img
    # TE : echo times [sec]
    #FieldStrength : Fieldstrength [T]
    #PrecessionIsClockwise: clockwise precession
    #Optionnal variable:
    #fieldmap : B0 fieldmap [Hz]


    Varnames=mat['imDataParams'].keys()

    if 'images' in Varnames:
        img=mat['imDataParams']['images']
    else:
        raise Exception('Warning: Not image in imDataParams' +
                        'Need to save your image (nx,ny,nz,ncoil,nt) as images')

    if 'TE' in Varnames:
        echoTimes=mat['imDataParams']['TE'] # TEs [sec]
    else:
        raise Exception('Warning: Not TE in imDataParams' +
                        'Need to save your echotimes (sec) as TE')

    if 'FieldStrength' in Varnames:
        dPar['B0']=mat['imDataParams']['FieldStrength'][0] # Fieldstrength [T]
    else:
        raise Exception('Warning: Not FieldStrength in imDataParams' +
                        'Need to save your fieldstrength (T) as FieldStrength')

    if 'PrecessionIsClockwise' in Varnames:
        clockwise=mat['imDataParams']['PrecessionIsClockwise'][0]  # Clockwiseprecession
    else:
        raise Exception('Warning: Not clockwise precession in imDataParams' +
                        'Need to save your clockwise precession as PrecessionIsClockwise')

    if 'unwrapped_phase' in Varnames:
        dPar['unwrapped_phase']=mat['imDataParams']['unwrapped_phase'][:,:,:,0,:]
        
    if 'voxelSize' in Varnames:
        dPar['dy'],dPar['dx'],dPar['dz']=mat['imDataParams']['voxelSize']

        
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



    # To get data as: (echo,slice,row,col)
    img.shape = (dPar['ny'], dPar['nx'], dPar['nz'], dPar['N'])
    img = np.transpose(img)
    img = np.swapaxes(img, 2, 3)

    img = img.flatten()
    dPar['img'] = img*dPar['reScale']


# Save output as MATLAB arrays
def save(output, dPar):
    if 'outName' in dPar:
        sse=output.pop("sse")
        algoParams=output.pop("algoParams")

        out_mat = {"params": output,"sse":sse,'algoParams':algoParams}
        print(op.abspath(dPar['outDir']))

        if not op.exists(op.abspath(dPar['outDir'])):
            os.mkdir(op.abspath(dPar['outDir']))
        filename=op.join(op.abspath(dPar['outDir']),'{}.mat'.format(dPar['outName']))
        scipy.io.savemat(filename, out_mat)
    else:
        dPar['outDir'].mkdir(parents=True, exist_ok=True)
        filename = dPar['outDir'] / './{}.mat'.format(dPar['sliceList'][0])
        print(r'Writing images to "{}"'.format(filename))
        scipy.io.savemat(filename, output)
