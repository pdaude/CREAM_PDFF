import pydicom
import datetime
import numpy as np
from pathlib import Path


gyro = 42.58  # 1H gyromagnetic ratio


# Dictionary of DICOM tags
tagDict = {
    'Image Type': 0x00080008,
    'SOP Class UID': 0x00080016,
    'SOP Instance UID': 0x00080018,
    'Series Description': 0x0008103E,
    'Slice Thickness': 0x00180050,
    'Spacing Between Slices': 0x00180088,
    'Echo Time': 0x00180081,
    'Imaging Frequency': 0x00180084,
    'Protocol Name': 0x00181030,
    'Study Instance UID': 0x0020000D,
    'Series Instance UID': 0x0020000E,
    'Series Number': 0x00200011,
    'Slice Location': 0x00201041,
    'Image Position (Patient)': 0x00200032,
    'Rows': 0x00280010,
    'Columns': 0x00280011,
    'Pixel Spacing': 0x00280030,
    'Pixel Aspect Ratio': 0x00280034,
    'Smallest Pixel Value': 0x00280106,
    'Largest Pixel Value': 0x00280107,
    'Window Center': 0x00281050,
    'Window Width': 0x00281051,
    'Rescale Intercept': 0x00281052,
    'Rescale Slope': 0x00281053,
    'Number of frames': 0x00280008,
    'Frame sequence': 0x52009230}  # Per-frame Functional Groups Sequence


# DICOM tags for output from fat/water separation
imgTypes = {
    'phi': {'descr': 'Initial phase (degrees)', 'seriesNumber': 100},
    'wat': {'descr': 'Water-only', 'seriesNumber': 101},
    'fat': {'descr': 'Fat-only', 'seriesNumber': 102},
    'ip': {'descr': 'In-phase', 'seriesNumber': 103},
    'op': {'descr': 'Opposed-phase', 'seriesNumber': 104},
    'ff': {'descr': 'Fat Fraction (%)', 'seriesNumber': 105, 'reScaleIntercept': -100, 'reScaleSlope': 0.1},
    'R2map': {'descr': 'R2* (msec-1)', 'seriesNumber': 106},
    'B0map': {'descr': 'B0 inhomogeneity (ppm)', 'seriesNumber': 107, 'reScaleSlope': 0.001},
    'CL': {'descr': 'FAC Chain length', 'seriesNumber': 108, 'reScaleSlope': 0.01},
    'UD': {'descr': 'FAC Unsaturation degree', 'seriesNumber': 109, 'reScaleSlope': 0.01},
    'PUD': {'descr': 'FAC Polyunsaturation degree', 'seriesNumber': 110, 'reScaleSlope': 0.01}
}


# List of DICOM attributes required for the water-fat separation
reqAttributes = ['Image Type', 
                 'Echo Time', 
                 'Slice Location',
                 'Imaging Frequency', 
                 'Columns', 
                 'Rows',
                 'Pixel Spacing',
                 'Spacing Between Slices']


# Translates series description tag to M/P/R/I for magn/phase/real/imaginary
def seriesDescription2type(seriesDescription):
    for type, description in [('R', 'Real Image'), ('I', 'Imag Image')]:
        if description in seriesDescription:
            return type
    return None


# Translates image type tag to M/P/R/I for magnitude/phase/real/imaginary
def typeTag2type(tagValue):
    for type in ['M', 'P', 'R', 'I']:
        if type in tagValue:
            return type
    return None


# Retrieves DICOM element value from dataset ds at tag=key.
# Use frame for multiframe DICOM files
def getTagValue(ds, key, frame=None):
    # Philips(?) private tag containing frame tags
    if (frame is not None and
       0x2005140f in ds[tagDict['Frame sequence']].value[frame]):
        frameObject = ds[tagDict['Frame sequence']].value[frame][0x2005140f][0]
        if tagDict[key] in frameObject:
            return frameObject[tagDict[key]].value
    if tagDict[key] in ds:
        return ds[tagDict[key]].value
    return None


# Retrieve attribute from dataset ds. Use frame for multiframe DICOM files.
# Must attributes are read directly from their corresponding DICOM tag
def getAttribute(ds, attr, frame=None):
    attribute = getTagValue(ds, attr, frame)
    if attr == 'Slice Location' and attribute is None:
        attribute = getTagValue(ds, 'Image Position (Patient)', frame)
        if attribute:
            attribute = attribute[2]
    elif attr == 'Spacing Between Slices' and attribute is None:
        attribute = getTagValue(ds, 'Slice Thickness', frame)
    elif attr == 'Image Type' and attribute:  # special handling of image type
        attribute = typeTag2type(attribute)
        if not attribute:
            attribute = seriesDescription2type(
                getTagValue(ds, 'Series Description', frame))
    return attribute


# Check if attribute is in DICOM dataset ds
def AttrInDataset(ds, attr, multiframe):
    if getAttribute(ds, attr) is not None:
        return True
    elif multiframe:
        for frame in range(len(ds[tagDict['Frame sequence']].value)):
            if not getAttribute(ds, attr, frame):
                return False  # Attribute must be in all frames!
        return True
    return False


# Check if ds is a multiframe DICOM object
def isMultiFrame(ds):
    return tagDict['Number of frames'] in ds and \
        int(ds[tagDict['Number of frames']].value) > 1 and \
        tagDict['Frame sequence'] in ds


# Extract files that are readable and have all required DICOM tags
def getValidFiles(files, printOutput=False):
    validFiles = []
    for file in files:
        try:
            ds = pydicom.read_file(str(file), stop_before_pixels=True)
        except:
            if printOutput:
                print('Could not read file: {}'.format(file))
            continue
        multiframe = isMultiFrame(ds)
        hasRequiredAttrs = [AttrInDataset(ds, attr, multiframe)
                            for attr in reqAttributes]
        if not all(hasRequiredAttrs):
            if printOutput:
                print('File {} is missing required DICOM tags:'.format(file))
                for i, hasAttr in enumerate(hasRequiredAttrs):
                    if not hasAttr:
                        print(reqAttributes[i])
            continue
        else:
            validFiles.append(file)
    return validFiles


# get combination of image types for DICOM frames in frameList
def getType(frameList, printType=False):
    typeTags = [tags[2] for tags in frameList]
    numR = typeTags.count('R')
    numI = typeTags.count('I')
    numM = typeTags.count('M')
    numP = typeTags.count('P')
    if numM+numP == 0 and numR+numI > 0 and numR == numI:
        if printType:
            print('Real/Imaginary images')
        return 'RI'
    elif numM+numP > 0 and numR+numI == 0 and numM == numP:
        if printType:
            print('Magnitude/Phase images')
        return 'MP'
    elif numP == 0 and numM+numR+numI > 0 and numM == numR == numI:
        if printType:
            print('Magnitude/Real/Imaginary images')
        return 'MRI'
    elif numM+numP+numR+numI > 0 and numM == numP == numR == numI:
        if printType:
            print('Magnitude/Phase/Real/Imaginary images')
        return 'MPRI'
    else:
        raise Exception('Unknown combination of image types: ' +
                        '{} real, {} imag, {} magn, {} phase'
                        .format(numR, numI, numM, numP))


def getSOPInstanceUID():
    t = datetime.datetime.now()
    datestr = '{:04d}{:02d}{:02d}{:02d}{:02d}{:02d}{:03d}'.format(
     t.year, t.month, t.day, t.hour, t.minute, t.second, t.microsecond//1000)
    randstr = str(np.random.randint(1000, 1000000000))
    uidstr = "1.3.12.2.1107.5.2.32.35356." + datestr + randstr
    return uidstr


def getSeriesInstanceUID(dPar, seriesDescription):
    if not 'seriesInstanceUIDs' in dPar:
        dPar['seriesInstanceUIDs'] = {}
    if not seriesDescription in dPar['seriesInstanceUIDs']:
        dPar['seriesInstanceUIDs'][seriesDescription] = getSOPInstanceUID() + ".0.0.0"
    return dPar['seriesInstanceUIDs'][seriesDescription]


# update dPar with info retrieved from the DICOM files including image data
def updateDataParams(dPar, files):
    dPar['fileType'] = 'DICOM'
    frameList = []
    for file in files:
        ds = pydicom.read_file(str(file), stop_before_pixels=True)
        multiframe = isMultiFrame(ds)
        if multiframe:
            if len(files) > 1:
                raise Exception('Support for multiple multi-frame DICOM' +
                                'files not implemented yet!')
            for frame in range(len(ds[tagDict['Frame sequence']].value)):
                frameList.append([file]+[frame]+[getAttribute(ds, attr, frame)
                                 for attr in reqAttributes])
        else:  # Single frame DICOM files
            frameList.append([file]+[None]+[getAttribute(ds, attr)
                             for attr in reqAttributes])
    frameList.sort(key=lambda tags: tags[2])  # First, sort on type (M/P/R/I)
    frameList.sort(key=lambda tags: tags[3])  # Second, sort on echo time
    frameList.sort(key=lambda tags: tags[4])  # Third, sort on slice location

    type = getType(frameList, True)
    dPar['dx'] = float(frameList[0][8][1])
    dPar['dy'] = float(frameList[0][8][0])
    dPar['dz'] = float(frameList[0][9])

    dPar['B0'] = frameList[0][5]/gyro
    # [msec]->[sec]
    echoTimes = sorted(set([float(tags[3])/1000. for tags in frameList]))
    dPar['totalN'] = len(echoTimes)
    if 'echoes' not in dPar:
        dPar['echoes'] = range(dPar['totalN'])
    echoTimes = [echoTimes[echo] for echo in dPar['echoes']]
    dPar['N'] = len(dPar['echoes'])
    if dPar['N'] < 2:
        raise Exception('At least 2 echoes required, only {} given'
                        .format(dPar['N']))
    dPar['t1'] = echoTimes[0]
    dPar['dt'] = np.mean(np.diff(echoTimes))
    if np.max(np.diff(echoTimes))/dPar['dt'] > 1.05 or \
       np.min(np.diff(echoTimes))/dPar['dt'] < .95:
        print('Warning: echo inter-spacing varies more than 5%')
        print(echoTimes)
    nSlices = len(set([tags[4] for tags in frameList]))
    if 'sliceList' not in dPar:
        dPar['sliceList'] = range(nSlices)

    dPar['nx'] = frameList[0][6]
    dPar['ny'] = frameList[0][7]
    dPar['nz'] = len(dPar['sliceList'])

    if 'cropFOV' in dPar:
        x1, x2 = dPar['cropFOV'][0], dPar['cropFOV'][1]
        y1, y2 = dPar['cropFOV'][2], dPar['cropFOV'][3]
        dPar['Nx'], dPar['nx'] = dPar['nx'], x2-x1
        dPar['Ny'], dPar['ny'] = dPar['ny'], y2-y1
    else:
        x1, x2 = 0, dPar['nx']
        y1, y2 = 0, dPar['ny']
    img = []
    if multiframe:
        file = frameList[0][0]
        dcm = pydicom.read_file(str(file))
    for n in dPar['echoes']:
        for slice in dPar['sliceList']:
            i = (dPar['N']*slice+n)*len(type)
            if type == 'MP':  # Magnitude/phase images
                magnFrame = i
                phaseFrame = i+1
                if multiframe:
                    magn = dcm.pixel_array[
                        frameList[magnFrame][1]][y1:y2, x1:x2].flatten()
                    phase = dcm.pixel_array[
                        frameList[phaseFrame][1]][y1:y2, x1:x2].flatten()
                    # Abs val needed for Siemens data to get correct phase sign
                    reScaleIntercept = \
                        np.abs(getAttribute(dcm, 'Rescale Intercept',
                                            frameList[phaseFrame][1]))
                else:
                    magnFile = frameList[magnFrame][0]
                    phaseFile = frameList[phaseFrame][0]
                    mDcm = pydicom.read_file(str(magnFile))
                    pDcm = pydicom.read_file(str(phaseFile))
                    magn = mDcm.pixel_array[y1:y2, x1:x2].flatten()
                    phase = pDcm.pixel_array[y1:y2, x1:x2].flatten()
                    # Abs val needed for Siemens data to get correct phase sign
                    try:
                        reScaleIntercept = np.abs(
                            getAttribute(pDcm, 'Rescale Intercept'))
                    except:
                        print('No Rescale Intercept DICOM tag found. Using 4096.')
                        reScaleIntercept = 4096
                # For some reason, intercept is used as slope (Siemens only?)
                c = magn*np.exp(phase/float(reScaleIntercept)*2*np.pi*1j)
            # Real/imaginary images and Magnitude/real/imaginary images
            elif type in ['RI', 'MRI', 'MPRI']:
                if type == 'RI':
                    realFrame = i+1
                elif type == 'MRI':
                    realFrame = i+2
                elif type == 'MPRI':
                    realFrame = i+3
                imagFrame = i
                if multiframe:
                    realPart = dcm.pixel_array[
                        frameList[realFrame][1]][y1:y2, x1:x2].flatten()
                    imagPart = dcm.pixel_array[
                        frameList[imagFrame][1]][y1:y2, x1:x2].flatten()
                    # Assumes real and imaginary slope/intercept are equal
                    reScaleIntercept = getAttribute(
                        dcm, 'Rescale Intercept', frameList[realFrame][1])
                    reScaleSlope = getAttribute(
                        dcm, 'Rescale Slope', frameList[realFrame][1])
                else:
                    realFile = frameList[realFrame][0]
                    imagFile = frameList[imagFrame][0]
                    rDcm = pydicom.read_file(str(realFile))
                    iDcm = pydicom.read_file(str(imagFile))
                    realPart = rDcm.pixel_array[y1:y2, x1:x2].flatten()
                    imagPart = iDcm.pixel_array[y1:y2, x1:x2].flatten()
                    # Assumes real and imaginary slope/intercept are equal
                    reScaleIntercept = getAttribute(rDcm, 'Rescale Intercept')
                    reScaleSlope = getAttribute(rDcm, 'Rescale Slope')
                if reScaleIntercept and reScaleSlope:
                    offset = reScaleIntercept/reScaleSlope
                else:
                    offset = -2047.5
                c = (realPart+offset)+1.0*1j*(imagPart+offset)
            else:
                raise Exception('Unknown image types')
            img.append(c)
    dPar['frameList'] = frameList
    dPar['img'] = np.array(img)*dPar['reScale']


# Set window so that percentile % of pixels are inside
def getPercentileWindow(im, intercept, slope, percentile=95):
    lims = np.percentile(im, [(100-percentile)/2, percentile + (100-percentile)/2])
    width = lims[1]-lims[0]
    center = width/2.+lims[0]
    return center*slope+intercept, width*slope


# Sets DICOM element value in dataset ds at tag=key. Use frame for multiframe
# DICOM files. If missing tag, a new one is created if value representation VR
# is provided
def setTagValue(ds, key, val, frame=None, VR=None):
    # Philips(?) private tag containing frame tags
    if (frame is not None and
       0x2005140f in ds[tagDict['Frame sequence']].value[frame]):
        frameObject = ds[tagDict['Frame sequence']].value[frame][0x2005140f][0]
        if tagDict[key] in frameObject:
            frameObject[tagDict[key]].value = val
            return True
    if tagDict[key] in ds:
        ds[tagDict[key]].value = val
        return True
    # Else, add as new DICOM element:
    if VR:
        # Philips(?) private tag containing frame tags
        if (frame is not None and
           0x2005140f in ds[tagDict['Frame sequence']].value[frame]):
            frameObject = ds[tagDict['Frame sequence']].\
                value[frame][0x2005140f][0]
            frameObject.add_new(tagDict[key], VR, val)
            return True
        ds.add_new(tagDict[key], VR, val)
        return True
    # print('Warning: DICOM tag {} was not set'.format(key))
    return False


# Save numpy array to DICOM image.
# Based on input DICOM image if exists, else create from scratch
def saveSeries(outDir, imgType, img, dPar):
    seriesDescription = imgTypes[imgType]['descr']
    seriesNumber = imgTypes[imgType]['seriesNumber']
    if 'reScaleIntercept' in imgTypes[imgType]:
        reScaleIntercept = imgTypes[imgType]['reScaleIntercept']
    else:
        reScaleIntercept = 0.
    if 'reScaleSlope' in imgTypes[imgType]:
        reScaleSlope = imgTypes[imgType]['reScaleSlope']
    else:
        reScaleSlope = 1.
    
    seriesInstanceUID = getSeriesInstanceUID(dPar, seriesDescription)
    # Single file is interpreted as multi-frame
    multiframe = dPar['frameList'] and \
        len(set([frame[0] for frame in dPar['frameList']])) == 1
    if multiframe:
        ds = pydicom.read_file(str(dPar['frameList'][0][0]))
        imVol = np.empty([dPar['nz'], dPar['ny']*dPar['nx']], dtype='uint16')
        frames = []
    if dPar['frameList']:
        DICOMimgType = getType(dPar['frameList'])
    for z, slice in enumerate(dPar['sliceList']):
        filename = outDir / './{}.dcm'.format(slice)
        # Extract slice, scale and type cast pixel data
        pixelData = np.array([max(0, (val-reScaleIntercept)/reScaleSlope)
                      for val in img[:, :, z].flatten()])
        pixelData = pixelData.astype('uint16')
        # Set window so that 95% of pixels are inside
        windowCenter, windowWidth = getPercentileWindow(
                                            pixelData, reScaleIntercept, reScaleSlope, 95)
        if dPar['frameList']:
            # Get frame
            frame = dPar['frameList'][dPar['totalN']*slice*len(DICOMimgType)]
            iFrame = frame[1]
            if not multiframe:
                ds = pydicom.read_file(str(frame[0]))
        else:
            iFrame = None
            # Create new DICOM images from scratch
            file_meta = pydicom.dataset.Dataset()
            file_meta.MediaStorageSOPClassUID = \
                'Secondary Capture Image Storage'
            file_meta.MediaStorageSOPInstanceUID = '1.3.6.1.4.1.9590.100.' +\
                '1.1.111165684411017669021768385720736873780'
            file_meta.ImplementationClassUID = '1.3.6.1.4.1.9590.100.' + \
                '1.0.100.4.0'
            file_meta.TransferSyntaxUID = pydicom.uid.ImplicitVRLittleEndian
            ds = pydicom.dataset.FileDataset(filename, {}, file_meta=file_meta,
                                           preamble=b"\0"*128)
            # Add DICOM tags:
            ds.Modality = 'WSD'
            ds.ContentDate = str(datetime.date.today()).replace('-', '')
            ds.ContentTime = str(datetime.time())  # millisecs since the epoch
            ds.SamplesPerPixel = 1
            ds.PhotometricInterpretation = "MONOCHROME2"
            ds.PixelRepresentation = 0
            ds.HighBit = 15
            ds.BitsStored = 16
            ds.BitsAllocated = 16
            ds.SmallestImagePixelValue = '\\x00\\x00'
            ds.LargestImagePixelValue = '\\xff\\xff'
            ds.Columns = img.shape[1]
            ds.Rows = img.shape[0]
            setTagValue(ds, 'Study Instance UID',
                        getSOPInstanceUID(), iFrame, 'UI')
            setTagValue(ds, 'Pixel Spacing', [dPar['dx'], dPar['dy']], iFrame, 'DS')
            setTagValue(ds, 'Pixel Aspect Ratio', [int(dPar['dx']*100), int(dPar['dy']*100)], iFrame, 'IS')
        # Change/add DICOM tags:
        setTagValue(ds, 'SOP Instance UID', getSOPInstanceUID(), iFrame, 'UI')
        setTagValue(ds, 'Series Instance UID', seriesInstanceUID, iFrame, 'UI')
        setTagValue(ds, 'Series Number', seriesNumber, iFrame, 'IS')
        setTagValue(ds, 'Echo Time', 0., iFrame, 'DS')
        setTagValue(ds, 'Protocol Name', 'Derived Image', iFrame, 'LO')
        setTagValue(ds, 'Series Description', seriesDescription, iFrame, 'LO')
        setTagValue(ds, 'Smallest Pixel Value', np.min(pixelData), iFrame)
        setTagValue(ds, 'Largest Pixel Value', np.max(pixelData), iFrame)
        setTagValue(ds, 'Window Center', int(windowCenter), iFrame, 'DS')
        setTagValue(ds, 'Window Width', int(windowWidth), iFrame, 'DS')
        setTagValue(ds, 'Rescale Intercept', reScaleIntercept, iFrame, 'DS')
        setTagValue(ds, 'Rescale Slope', reScaleSlope, iFrame, 'DS')

        if multiframe:
            imVol[z] = pixelData
            frames.append(iFrame)
        else:
            ds.PixelData = pixelData
            ds.save_as(str(filename))

    if multiframe:
        setTagValue(ds, 'SOP Instance UID', getSOPInstanceUID())
        setTagValue(ds, 'Number of frames', len(frames))
        ds[tagDict['Frame sequence']].value = \
            [ds[tagDict['Frame sequence']].value[frame] for frame in frames]
        ds.PixelData = imVol
        filename = outDir / './0.dcm'
        ds.save_as(str(filename))


# Save all data in output as DICOM images
def save(output, dPar):
    for seriesType in output:
        outDir = dPar['outDir'] / seriesType
        outDir.mkdir(parents=True, exist_ok=True)
        print(r'Writing image{} to "{}"'.format('s'*(dPar['nz'] > 1), outDir))
        saveSeries(outDir, seriesType, output[seriesType], dPar)
