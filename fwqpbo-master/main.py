#!/usr/bin/env python3

import numpy as np
import sys
import optparse
import config
import fatWaterSeparation
import DICOM
import MATLAB
import cream_pdff_r2.config_spectrum as cf



# Zero pad back any cropped FOV
def padCropped(croppedImage, dPar):
    if 'cropFOV' in dPar:
        image = np.zeros((dPar['nz'], dPar['Ny'], dPar['Nx']))
        x1, x2 = dPar['cropFOV'][0], dPar['cropFOV'][1]
        y1, y2 = dPar['cropFOV'][2], dPar['cropFOV'][3]
        image[:, y1:y2, x1:x2] = croppedImage
        return image
    else:
        return croppedImage


def save(output, dPar):
    for seriesType in output: # zero pad if was cropped and reshape to row,col,slice
        if seriesType!='algoParams':
            output[seriesType] = np.moveaxis(padCropped(output[seriesType].reshape((dPar['nz'], dPar['ny'], dPar['nx'])), dPar), 0, -1)
    
    if dPar['fileType'] == 'DICOM':
        DICOM.save(output, dPar)
    elif dPar['fileType'] == 'MATLAB':
        MATLAB.save(output, dPar)
    else:
        raise Exception('Unknown filetype: {}'.format(dPar['fileType']))


# Merge output for slices reconstructed separately
def mergeOutputSlices(outputList):
    mergedOutput = outputList[0]
    for output in outputList[1:]:
        for seriesType in output:
            mergedOutput[seriesType] = np.concatenate((mergedOutput[seriesType], output[seriesType]))
    return mergedOutput


def getFattyAcidComposition(rho):
    nFAC = len(rho) - 2 # Number of Fatty Acid Composition Parameters
    eps = sys.float_info.epsilon
    CL, UD, PUD = None, None, None

    if nFAC == 1:
        # UD = F2/F1
        UD = np.abs(rho[2] / (rho[1] + eps))
    elif nFAC == 2:
        # UD = F2/F1
        # PUD = F3/F1
        UD = np.abs(rho[2] / (rho[1] + eps))
        PUD = np.abs(rho[3] / (rho[1] + eps))
    elif nFAC == 3:
        # UD = F2/F1
        # PUD = F3/F1
        # CL = F4/F1
        UD = np.abs(rho[2] / (rho[1] + eps))
        PUD = np.abs(rho[3] / (rho[1] + eps))
        CL = np.abs(rho[4] / (rho[1] + eps))
    else:
        raise Exception('Unknown number of Fatty Acid Composition parameters: {}'.format(nFAC))

    return CL, UD, PUD


# Get total fat component (for Fatty Acid Composition; trivial otherwise)
def getFat(rho, alpha):
    nVxl = np.shape(rho)[1]
    fat = np.zeros(nVxl, dtype=complex)
    for m in range(1, alpha.shape[0]):
        fat += sum(alpha[m, 1:])*rho[m]
    return fat


# Perform fat/water separation and return prescribed output
def reconstruct(dPar, aPar, mPar):

    rho, B0map, R2map = fatWaterSeparation.reconstruct(dPar, aPar, mPar)
    wat = rho[0]
    fat = getFat(rho, mPar['alpha'])

    # Conversion ppm  to Hz
    larmor = mPar['gyro'] * dPar['B0']
    B0map = B0map * larmor

    B0map = B0map[np.newaxis,:]
    R2map = R2map[np.newaxis,:]

    #phi=2*piB0 +j*R2
    #Estimated signal: (YM)
    #YM= B*AX=B*K
    #B=exp(1j*te*phi)
    #A=alpha*exp(1j*te*omega)
    #X=[wat,fat]
    #Warning X=[wat,fat]*exp(1j*t0*phi) => B=exp(1j*dte*phi)

    phi=(2*np.pi*B0map +1j*R2map)

    te=np.array([dPar["t1"]+n*dPar["dt"] for n in range(dPar["N"])])[:,np.newaxis]
    dte=np.array([n*dPar["dt"] for n in range(dPar["N"])])[:,np.newaxis]
    B = np.exp(1j * dte * phi)
    
    X = np.vstack((wat, fat))
    omega = 2 * np.pi * larmor * (np.array(mPar['CS']) - mPar['CS'][0])[np.newaxis, :]
    A =  np.exp(1j * te * omega)*np.sum(mPar['alpha'],axis=0)
    A = np.hstack((A[:, 0][:, np.newaxis], np.sum(A[:, 1:], axis=1)[:, np.newaxis]))
    
    K = np.dot(A, X)
    
    YM = B * K #

    R =(dPar["img"] - YM)/dPar['reScale']
    #Sum of square error
    SSE=np.linalg.norm(R,axis=0)



    # Prepare prescribed output
    output = {}
    output['sse'] = SSE

    #output['YM'] = YM

    if 'ff' in aPar['output']: # Calculate the fat fraction
        if aPar['magnitudeDiscrimination']:  # to avoid bias from noise
            output['FF'] = 100 * np.real(fat / (wat + fat + sys.float_info.epsilon))
        else:
            output['FF'] = 100 * np.abs(fat)/(np.abs(wat) + np.abs(fat) + sys.float_info.epsilon)

    #Correct fat and water map (divide by exp(i(w+iwR2)t1)
    Corrected_term = np.exp(1j * dPar["t1"] * phi)
    wat=wat/Corrected_term
    fat=fat/Corrected_term

    if 'wat' in aPar['output']:
        output['W'] = wat
    if 'fat' in aPar['output']:
        output['F'] = fat


    if 'B0map' in aPar['output']:
        output['B0'] = B0map

    if 'R2map' in aPar['output']:
        output['R2'] = R2map

    output['algoParams'] = aPar
    
    #output['unwrapB0']=fatWaterSeparation.unwrap_B0(B0map,dPar)

    # if 'phi' in aPar['output']:
    #     output['PH'] = np.angle(wat, deg=True) + 180
    # if 'ip' in aPar['output']: # Calculate synthetic in-phase
    #     output['ip'] = np.abs(wat+fat)
    # if 'op' in aPar['output']: # Calculate synthetic opposed-phase
    #     output['op'] = np.abs(wat-fat)


    # Do any Fatty Acid Composition in a second pass
    # if mPar['nFAC'] > 0:
    #     rho = fatWaterSeparation.reconstruct(dPar, aPar['pass2'], mPar['pass2'], B0map, R2map)[0]
    #     CL, UD, PUD = getFattyAcidComposition(rho)
    
    #     if 'CL' in aPar['output']:
    #         output['CL'] = CL
    #     if 'UD' in aPar['output']:
    #         output['UD'] = UD
    #     if 'PUD' in aPar['output']:
    #         output['PUD'] = PUD

    return output


def main(dataParamFile, algoParamFile, modelParamFile, outDir=None):
    # Read configuration files
    dPar = config.readConfig(dataParamFile, 'data parameters')
    aPar = config.readConfig(algoParamFile, 'algorithm parameters')
    mPar = cf.readConfig(modelParamFile, 'model parameters')

    # Setup configuration objects
    config.setupDataParams(dPar, outDir)
    cf.setupModelParams(mPar, dPar['clockwisePrecession'], dPar['temperature'])
    config.setupAlgoParams(aPar, dPar['N'], mPar['nFAC'])

    print('mu = {}'.format(aPar['mu']))
    print('outName = {}'.format(dPar['outName']))
    print('B0 = {}'.format(round(dPar['B0'], 2)))
    print('N = {}'.format(dPar['N']))
    print('t1/dt = {}/{} msec'.format(round(dPar['t1']*1000, 2),round(dPar['dt']*1000, 2)))
    print('nx,ny,nz = {},{},{}'.format(dPar['nx'], dPar['ny'], dPar['nz']))
    print('dx,dy,dz = {},{},{}'.format(round(dPar['dx'], 2), round(dPar['dy'], 2), round(dPar['dz'], 2)))
    print("Rescaling factor {}".format(dPar['reScale']))
    # Run fat/water processing and save output
    if aPar['use3D'] or len(dPar['sliceList']) == 1:
        if 'slabs' in dPar:
            for iSlab, (slices, z) in enumerate(dPar['slabs']):
                print('Processing slab {}/{} (slices {}-{})...'
                      .format(iSlab+1, len(dPar['slabs']), slices[0]+1, slices[-1]+1))
                slabDataParams = config.getSlabDataParams(dPar, slices, z)
                output = reconstruct(slabDataParams, aPar, mPar)
                save(output, slabDataParams) # save data slab-wise to save memory
        else:
            output = reconstruct(dPar, aPar, mPar)
            save(output, dPar)
    else:
        output = []
        for z, slice in enumerate(dPar['sliceList']):
            print('Processing slice {} ({}/{})...'
                  .format(slice+1, z+1, len(dPar['sliceList'])))
            sliceDataParams = config.getSliceDataParams(dPar, slice, z)
            output.append(reconstruct(sliceDataParams, aPar, mPar))
        save(mergeOutputSlices(output), dPar)


if __name__ == '__main__':
    # Initiate command line parser
    p = optparse.OptionParser()
    p.add_option('--dataParamFile', '-d', default='',  type="string",
                 help="File path of data parameter configuration file")
    p.add_option('--algoParamFile', '-a', default='',  type="string",
                 help="File path of algorithm parameter configuration file")
    p.add_option('--modelParamFile', '-m', default='',  type="string",
                 help="File path of model parameter configuration file")

    # Parse command line
    options, arguments = p.parse_args()

    main(options.dataParamFile, options.algoParamFile, options.modelParamFile)