import numpy as np
from pathlib import Path
import yaml

import fwspectrum as fws
import inspect
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

def defaultModelParams():
    defaults = {
        'model':'Berglund2012',
        'watCS':4.7,
        'gyro': 42.577478
    }
    return defaults

# Update model parameter object mPar and set default parameters
def setupModelParams(mPar, clockwisePrecession=False, temperature=None):

    defaultmPar = defaultModelParams()

    for param, defval in defaultmPar.items():
        if param not in mPar:
            mPar[param] = defval

    del defaultmPar
    if temperature:  # Temperature dependence according to Hernando 2014
            mPar['watCS'] = 1.3 + 3.748 - .01085 * temperature  # Temp in [°C]

    #List of  all classes in fwspectrum.py
    FWspectrum_classes=inspect.getmembers(fws, lambda ele: inspect.isclass(ele))
    #Removing mother Class FWspectrum and simplifying class name
    FWspectrum_classes=[(FWname.split('_')[0],FWclass) for FWname,FWclass in FWspectrum_classes if '_' in FWname]
    FWspectrums=dict(FWspectrum_classes)
    if mPar['model'] in FWspectrums.keys():
        if mPar['model']!="Custom" and ('FatCS' in mPar.keys() or 'realAmps' in mPar.keys()):
            raise Exception("You can't modify FatCS and realAmps of existing models. \n "
                            "Please use the model Custom to do it")
        FatWaterspectrum=FWspectrums[mPar['model']](**mPar)
    else:
        raise Exception('Unknown spectrum model :  \n Available models are {} '.format(FWspectrums.keys()))

    mPar['fatCS']=FatWaterspectrum.getFatFreq()
    mPar['CS'] = np.array([mPar['watCS']] + mPar['fatCS'], dtype=np.float32)

    if clockwisePrecession:
        mPar['CS'] *= -1

    # Number of resonnance peaks
    mPar['P'] = len(mPar['CS'])

    #Get number of unknown fat variables
    mPar['nFAC']=FatWaterspectrum.getUFV()

    mPar['M'] = 2 + mPar['nFAC']

    mPar['alpha'] =FatWaterspectrum.getAlpha()

    del FatWaterspectrum
    #In algoparams
    species=[{"name":'water',"frequency":mPar["watCS"],"relAmps":mPar["alpha"][0,:]},{"name":'fat',"frequency":mPar["fatCS"],"relAmps":mPar["alpha"][1,:]}]
    return species


