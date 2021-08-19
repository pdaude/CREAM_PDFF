import os
import main
import os.path as op
if __name__ == '__main__':

    CREAM_PDFFPath = op.dirname(op.dirname(op.abspath(__file__)))
    algoPath=op.join(CREAM_PDFFPath,'algoParams')

    modelPath = op.join(CREAM_PDFFPath, 'modelParams')
    input_folder = '/home/pdaude/Matlab_dev/20210127_AD/romeo_yml'


    modelParamsFile = op.join(modelPath,'CustommodelParams.yml')
    algoParamsFile = op.join(algoPath,'BerglundFWQPBO_algoParams.yml')

    dataParamsFile = op.join(CREAM_PDFFPath,'simu_2020bydder.yml')
    main.main(dataParamsFile, algoParamsFile, modelParamsFile)
    # for ymlfile in os.listdir(input_folder):
    #     if ymlfile.split('.')[-1]=="yml":
    #         dataParamsFile = op.join(input_folder,ymlfile)
    #         main.main(dataParamsFile, algoParamsFile, modelParamsFile)
