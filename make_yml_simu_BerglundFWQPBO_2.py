import os
import os.path as op
import yaml
if __name__ == '__main__':
    gitPath='/home/pdaude/Projet_Python/CREAM_PDFF/'
    simu_path=op.join(gitPath,'simulation')
    simu_ImDataParams_yml='simu_yml_BerglundFWQPBO'
    simuImDataParams_dir='simuImDataParams'

    simuImDataParams_path = op.join(simu_path,simuImDataParams_dir)
    simu_ImDataParams_yml_path = op.join(simu_path,simu_ImDataParams_yml)


    if not op.exists(simu_ImDataParams_yml_path):
        os.mkdir(simu_ImDataParams_yml_path)


    for matfile in os.listdir(simuImDataParams_path):
        input_path_matfile = op.join(simuImDataParams_path , matfile)
        output_simu_dir =matfile.split('.')[0] # Remove extension .mat
        
        #Simulation with two spectrums (Good or Wrong)
        ExtensionNames=['','Wrong']
        for exN in ExtensionNames:
            output_simu_path= op.join(simu_path,'{}{}'.format(exN,output_simu_dir))
        
            output_simu_name='BerglundFWQPBO_{}{}'.format(exN,output_simu_dir)

            ### Could create directory output simu here
            #if not op.exists(output_simu_path):
            #os.mkdir(output_simu_path)
            ###
            

            ymldict={"files":[input_path_matfile],
                    "outDir":output_simu_path,
                    "outName":output_simu_name,
                    "reScale": 1,
                    }
            ymlfile=op.join(simu_ImDataParams_yml_path,'{}.yml'.format(output_simu_name))

            with open(ymlfile, 'w') as file:
                documents = yaml.dump(ymldict, file,sort_keys=False)

            del ymldict
