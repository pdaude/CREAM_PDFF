clear all
clc
close all

%% Path architecture
gitPath='/home/pdaude/Projet_Python/CREAM_PDFF/';
imDataParamsFolder='simuImDataParams';
gtFolder='simuGT';
modelParamsFolder='modelParams';


simuPath=fullfile(gitPath,'simulationR2');
imDataParamsPath=fullfile(simuPath,imDataParamsFolder);
gtPath=fullfile(simuPath,gtFolder);

if not(isfolder(simuPath))
    mkdir(simuPath)
end

if not(isfolder(imDataParamsPath))
    mkdir(imDataParamsPath)
end

if not(isfolder(gtPath))
    mkdir(gtPath)
end

%% Simulated Parameters

%Simulation R2
b0=2.895; %(T) %VIDA B0 123.260178 (MHz) -> 2.895 T
spectrumName='Hodson2008';
TEtypes={'IN/OPP','IDEAL','MINIMUM'};% TE IN/OPP ; IDEAL  
%TEtypes={'IDEAL','MINIMUM'};% TE IN/OPP ; IDEAL 
NTEs=[5,7]; %number of TE 3 5 7 9
SNRs=[50]; % SNRs 50-100-10
zNrep=100; %Number of repetition
xFF=0:1:100; %
phi=pi/6; % 30°
r2stars=[150,300]; % T2*=30ms
yB0=-300:6:300;
padding=5; % Avoid border effect
dTEmin=1.68e-3;
TEmin=0.98e-3;


%Avoiding borders effects
xFF=padarray(xFF,[0,padding],'replicate');
yB0=padarray(yB0,[0,padding],'replicate');

for TEtype=TEtypes
    TEtype=char(TEtype);
    for NTE=NTEs
        for SNR=SNRs
            for r2star=r2stars
                %Reproductibility
                rng(0,'twister')
                clear imDataParams algoParams params sse


                %Define imDataParams
                imDataParams.FieldStrength =b0;
                imDataParams.PrecessionIsClockwise = 1;
                imDataParams.voxelSize =[1 1 1];

                % Define spectrum 
                modelParamsPath=fullfile(gitPath,modelParamsFolder,join([spectrumName,modelParamsFolder,'.yml']));
                modelParams =ReadYaml(modelParamsPath);
                [species ,FWspectrum]= setupModelParams(modelParams);

                algoParams.gyro=FWspectrum.gyro;
                algoParams.species=species;

                larmor = algoParams.gyro*imDataParams.FieldStrength; % larmor freq (MHz)


                % Define TE

                %bydder TE (0.20,0.56,0.93,1.30,1.66,2.03)pi
                %TE=[0.20,0.56,0.93,1.30,1.66,2.03]/(larmor*(4.7-1.3));
                TE= simuGenericTE(NTE,imDataParams.FieldStrength,TEtype,'dTEmin',dTEmin,'TEmin',TEmin);

                switch TEtype

                case 'IN/OPP'

                    TEabrev='IO';
                case 'IDEAL'

                     TEabrev='ID';

                case 'MINIMUM'

                    TEabrev='MI';
                otherwise
                    error('type should be IDEAL and IN/OPP or MINIMUM');
                end


                TEname=sprintf('TE%d%s',NTE,TEabrev);

                imDataParams.TE=TE ;


                % Define imDataParams.images
                [XFF,YB0,ZRep]=meshgrid(xFF,yB0,ones(zNrep,1));
                N=numel(XFF);
                x=numel(xFF);
                y=numel(yB0);
                z=zNrep;

                B0=reshape(YB0,[1,N]);
                R2=r2star*ones(1,N);
                PHI=phi*ones(1,N);
                W=reshape(1-XFF/100,[1,N]);
                F=reshape(XFF/100,[1,N]);

                FreqA=(species(2).frequency-species(1).frequency);
                FatA = species(2).relAmps.*exp(2*pi*1i*larmor*TE'*FreqA);


                A = [ones(numel(TE),1), sum(FatA,2)];
                B = exp(1i*(TE'*complex(2*pi*B0,R2)+repmat(PHI,numel(TE),1)));
                X=[W;F];

                YM=B.*(A*X);
                img = reshape(YM,numel(TE),x,y,z,1);
                image=permute(img,[2,3,4,5,1]);
                %√(2−π/2) ≈ 0.66.
                %The measured SNR = S/N must then be multiplied by the 0.66 Rayleigh distribution correction factor to calculate the true SNR.
                %Create NOISE with a SNR given
                %Simulation 1 2 3
                Smean=mean(abs(image(:,:,:,1,1)),'all');
                AmpN=Smean/SNR/sqrt(2-pi/2);
                Noise_signal=AmpN*complex(randn(size(image)),randn(size(image)));
                SNRmeas=Smean/std(abs(Noise_signal(:)));
                disp(SNRmeas)



                % image+Noise
                imDataParams.images=image+Noise_signal;

                %Normalise over 1 / 99% instead of max to avoid outliers.

                img_r=prctile(abs(imDataParams.images(:,:,:,1,1)),99,'all');

                if img_r==0;img_r=1;end % If img_r is null;

                imDataParams.images=imDataParams.images/img_r;

                %Saving ImdataParams & Ground truth (GT) params/sse


                simuImDataParamsPath=fullfile(imDataParamsPath,sprintf('%s_SNR%d_%s_%d_%s',imDataParamsFolder,SNR,TEname,r2star,spectrumName));

                save(simuImDataParamsPath,'imDataParams');

                % Ground truth params and sse;
                params.B0=YB0;
                params.R2=reshape(R2,[x,y,z]);
                params.F=reshape(F,[x,y,z])/img_r;
                params.W=reshape(W,[x,y,z])/img_r;
                params.PH=reshape(PHI,[x,y,z]);
                params.FF=XFF;
                [sse,~]=calculate_residual(params,algoParams,imDataParams);

                simuGTPath=fullfile(gtPath,sprintf('%s_SNR%d_%s_%d_%s',gtFolder,SNR,TEname,r2star,spectrumName));

                save(simuGTPath,'params','sse','algoParams');
            end
            
        end
            
    end
    
end


