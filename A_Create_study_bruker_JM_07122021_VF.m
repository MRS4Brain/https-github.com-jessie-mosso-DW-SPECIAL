%% From a PV 360 v1 Bruker exp. folder, generates MATLAB structures with necessary 
%params + data (fid or ser file) + fid.refscan (fidrefscan) if it exists. 
% Jessie Mosso @ EPFL - 17/02/2021
%includes: normalization, reading the files not opened in topspin, grpdly, offset 

clear; close all;

%% DO ADAPT THIS SECTION: Enter exp nb. + directory and run the program 
folderexp="G:\41-DEMO_DWSPECIAL\";
expnb=18;
timeacq='';%'_BDL727';%'' or am or pm
dateexp='2024'; 

%% DO NOT ADAPT THE FOLLOWING
%% Parameters - method file 
study.path=char(folderexp+num2str(expnb));
methodfile = fileread(folderexp +num2str(expnb)+ "\method");

%date and time
startind_date=strfind(methodfile,['$$ ' dateexp]); %adapt the data here is needed
date=methodfile(startind_date+3:startind_date+21);
study.time=date;

%literal names = acquisition time (sec), spectral width (Hz), sequence, np
%(complex nb *2), nucleus, resonance frequency B0, ppmoffset, tr(sec), RG,
% voxel dimensions (mm), acq type (MRS), nrep, nav
BRUKERparamnames=["PVM_SpecAcquisitionTime=","##$PVM_SpecSWH=( 1 )", ...
    "##$Method=","##$PVM_SpecMatrix=( 1 )","##$PVM_Nucleus1Enum=", ...
    ["##$PVM_FrqRef=( 8 )" + sprintf('\n')], ...
    ["##$PVM_FrqWorkOffsetPpm=( 8 )" + sprintf('\n')], ...
    "##$PVM_RepetitionTime=","##$PVM_RgValue=", ...
    ["##$PVM_VoxArrSize=( 1, 3 )" + sprintf('\n')],"##$PVM_EncSpectroscopy=", ...
    "PVM_NRepetitions=","PVM_NAverages=","PVM_EchoTime="];

GUIparamnames=["acq_time","spectralwidth","sequence","np","nucleus",...
    "resfreq","ppm_ref","tr","gain","voxs","acqtype","nrep","nav","te"];
numerical=[true, true, false, true, false, true, true, true, true,false, ...
    false, true, true,true];
arraylistbol=[false, false, false, false, false, false, false, false, ...
    false, true, false, false, false,false];

for par=1:length(BRUKERparamnames)
    startind=strfind(methodfile,BRUKERparamnames(par));
    if isempty(startind)
    else
    startind=startind(1);
    param=methodfile(startind+strlength(BRUKERparamnames(par)));
    
    if arraylistbol(par)
        kparam=1;
        while methodfile(startind+strlength(BRUKERparamnames(par))+kparam+1)~='#'
            param=[param methodfile(startind+strlength(BRUKERparamnames(par))+kparam)];
            kparam=kparam+1;
        end
        param=strsplit(param);
    else
        kparam=1;
%         disp(methodfile(startind+strlength(BRUKERparamnames(par))+kparam))
        while methodfile(startind+strlength(BRUKERparamnames(par))+kparam)~=['#',' ','$',sprintf('\n')]
            param=[param methodfile(startind+strlength(BRUKERparamnames(par))+kparam)];
            kparam=kparam+1;
        end
    end
    if numerical(par)
        param=str2double(param);
    end
    varname=matlab.lang.makeValidName(GUIparamnames(par));
    study.(varname{1})=param; 
    end 
end 

%adjustements
study.acq_time=study.acq_time*10^-3; % in seconds
study.np=study.np*2; %real/imag
study.tr=study.tr/1000; %in seconds 
if study.acqtype=='Yes'
    study.acqtype='MRS';

else 
    disp("can't open this type of data")
end
study.format='Matlab';
study.params.sfrq=study.resfreq;
study.params.arraydim=1; %MRS
study.params.np=study.np;
study.params.sw=study.spectralwidth;
study.params.tr=study.tr;
study=rmfield(study,'tr');
study.params.te=study.te; 
stuy=rmfield(study,'te');
study.params.gain=study.gain;
study=rmfield(study,'gain');
study.params.vox1=str2double(study.voxs(1));
study.params.vox2=str2double(study.voxs(2));
study.params.vox3=str2double(study.voxs(3));
study=rmfield(study,'voxs');

%% Group Delay - acqus file
if isfile(folderexp +num2str(expnb)+ "\acqus")
    acqusfile = fileread(folderexp +num2str(expnb)+ "\acqus");
    startind_grpdly=strfind(acqusfile,"##$GRPDLY= ");
    startind_grpdly=startind_grpdly(1);
    grpdly=acqusfile(startind_grpdly+strlength("##$GRPDLY= "));
    kgrp=1;
    while acqusfile(startind_grpdly+strlength("##$GRPDLY= ")+kgrp+1)~='#'
        grpdly=[grpdly acqusfile(startind_grpdly+strlength("##$GRPDLY= ")+kgrp)];
        kgrp=kgrp+1;
    end
    grpdly=str2double(grpdly);

elseif isfile(folderexp +num2str(expnb)+ "\acqp")
    acqpfile = fileread(folderexp +num2str(expnb)+ "\acqp");
    startind_grpdly=strfind(acqpfile,"##$ACQ_RxFilterInfo=( 2 )" + sprintf('\n') + "(");
    startind_grpdly=startind_grpdly(1);
    grpdly=acqpfile(startind_grpdly+strlength("##$ACQ_RxFilterInfo=( 2 )" + sprintf('\n') + "("));
    kgrp=1;
    while acqpfile(startind_grpdly+strlength("##$ACQ_RxFilterInfo=( 2 )" + sprintf('\n') + "(")+kgrp+1)~='#'
        grpdly=[grpdly acqpfile(startind_grpdly+strlength("##$ACQ_RxFilterInfo=( 2 )" + sprintf('\n') + "(")+kgrp)];
        kgrp=kgrp+1;
    end
    grpdly=strsplit(grpdly);
    grpdly=str2double(grpdly{1});
end
study.params.grpdly=grpdly;

%% ADCoverflow? - acqp file

acqpfile = fileread(folderexp +num2str(expnb)+ "\acqp");
startind_overflow=strfind(acqpfile,"##$ACQ_adc_overflow=( 2 )" + sprintf('\n')); %2 for 2 receiver channels
if isempty(startind_overflow)==1
    study.params.adcoverflow='No';
else
startind_overflow=startind_overflow(1);
overflow=acqpfile(startind_overflow+strlength("##$ACQ_adc_overflow=( 2 )" + sprintf('\n')));
kov=1;
while acqpfile(startind_overflow+strlength("##$ACQ_adc_overflow=( 2 )" + sprintf('\n')) +kov+1)~='#'
    overflow=[overflow acqpfile(startind_overflow+strlength("##$ACQ_adc_overflow=( 2 )" + sprintf('\n'))+kov)];
    kov=kov+1;
end
overflow=strsplit(overflow);
if strcmp(overflow{1},'Yes') | strcmp(overflow{2},'Yes')
    study.params.adcoverflow='Yes';
    f=msgbox(['ADC overflow during acquisition E' num2str(expnb)]);
else
    study.params.adcoverflow='No';
end 
end 


%% Read data

disp(study.sequence)

%job0
if isfile(folderexp +num2str(expnb)+ "\ser")
    disp ("2D ser, " + "nav=" + num2str(study.nav) + " nrep=" + num2str(study.nrep) + " saved")
    fileid=fopen(folderexp +num2str(expnb)+ "\ser",'r','ieee-le'); %read binary format
    if fileid == -1
        disp('Cannot open file');
        return
    end
else
    if study.nrep==1
        dimtype="1D";
    else
        dimtype="2D";
    end 
    disp (dimtype + " fid, " + "nav=" + num2str(study.nav) + " nrep=" + num2str(study.nrep)+ " saved")
    fileid=fopen(folderexp +num2str(expnb)+ "\fid",'r','ieee-le'); %read binary format
    if fileid == -1
        disp('Cannot open file');
        return
    end
end

buffer=fread(fileid,'int32'); 
nbptsfid=length(buffer)/(2*study.nrep);

for rep=1:study.nrep
    buffer_ser(rep,:)=buffer((rep-1)*(nbptsfid*2)+1:rep*(nbptsfid*2))';
end
ser_c=buffer_ser(:,1:2:end)+1i*buffer_ser(:,2:2:end);
fclose(fileid);

%offset 
if study.ppm_ref ~= 0
    offset_hz=study.ppm_ref*study.resfreq;
    dw=1/study.spectralwidth;
    t=[0:dw:(study.np/2-1)*dw];
    tmat=repmat(t,study.nrep,1);
    ser_c_shift=ser_c.*exp(1i.*(2*pi*offset_hz).*tmat);
    study.data.real(:,1,:)=real(ser_c_shift);
    study.data.imag(:,1,:)=-imag(ser_c_shift); %flips the spectrum
else
    study.data.real(:,1,:)=real(ser_c);
    study.data.imag(:,1,:)=-imag(ser_c); %flips the spectrum
end 

%filename and liststring
if study.nrep>1
    filename=['Bruker_' date(1:10)  '_'  num2str(expnb) timeacq '_ser.mat'];
    %nt
    study.params.nt=study.nrep*study.nav;

else
    filename=['Bruker_' date(1:10)  '_'  num2str(expnb) timeacq '_fid.mat'];
    %nt
    study.params.nt=study.nav;
end 

%normalize: xNA/RG/voxelsize
voxvol=study.params.vox1.*study.params.vox2.*study.params.vox3;
study.data.real=study.data.real.*study.nav./study.params.gain./voxvol;
study.data.imag=study.data.imag.*study.nav./study.params.gain./voxvol;

study.filename=filename;
study.liststring=char(folderexp + study.filename);

%multiplicity
study.multiplicity=study.nrep;

%process
study.process.lsfid=round(grpdly)+1;
study.process.apodizefct='exponential'; %default
study.process.apodparam1= zeros(1,study.nrep);%default
study.process.apodparam2=zeros(1,study.nrep);%default
study.process.transfsize=0; %default
study.process.appltoarray1=0;%default
study.process.appltoarray2=0;%default
study.process.phasecorr0=zeros(1,study.nrep);%default
study.process.phasecorr1=zeros(1,study.nrep);%default
study.process.DCoffset=0; %default


study=rmfield(study,'nav');
study=rmfield(study,'nrep');

save(folderexp + filename,'study')

%simple display
figure()
if study.multiplicity>1
    fid=squeeze(study.data.real)+1i*squeeze(study.data.imag); 
elseif study.multiplicity==1
    fid=squeeze(study.data.real)'+1i*squeeze(study.data.imag)'; 
end 
fid_sum=sum(fid(:,:),1);
fid_sum=[fid_sum(round(grpdly)+1:end), zeros(1,round(grpdly))];
ft=fftshift(fft(fid_sum./study.params.nt,[],2),2);
plot(real(ft),'DisplayName','job0')
hold on %in case of presence of refscan
legend()
title(num2str(expnb)) 


%% Refscan
if isfile(folderexp + num2str(expnb)+ "\fid.refscan")

    BRUKERparamnames_refscan=["PVM_RefScanNA=","##$PVM_RefScanRG="];
    localparamnames_refscan=["NArefscan","RGrefscan"];
    numerical=[true, true];

    for par=1:length(BRUKERparamnames_refscan)
        startind=strfind(methodfile,BRUKERparamnames_refscan(par));
        startind=startind(1);
        param=methodfile(startind+strlength(BRUKERparamnames_refscan(par)));
        kparam=1;
        while methodfile(startind+strlength(BRUKERparamnames_refscan(par))+kparam+1)~=['#',' ','$']
            param=[param methodfile(startind+strlength(BRUKERparamnames_refscan(par))+kparam)];
            kparam=kparam+1;
        end
        param=str2double(param);
        varname=matlab.lang.makeValidName(localparamnames_refscan(par));
        refscan.(varname{1})=param;    
    end 

    %store data
    disp('refscan saved')
    fileid=fopen(folderexp + num2str(expnb)+ "\fid.refscan",'r','ieee-le'); %read binary format
    if fileid == -1
        disp('Cannot open file');
        return
    end
    buffer=fread(fileid,'int32'); 
    nbptsfid=length(buffer)/2;
    buffer_c_ref=buffer(1:2:end)+1i*buffer(2:2:end); 
    fclose(fileid);

    study.data.real=zeros(1,1,study.np/2);
    study.data.imag=zeros(1,1,study.np/2);
    study.data.real(1,1,:)=real(buffer_c_ref);
    study.data.imag(1,1,:)=-imag(buffer_c_ref); %flips the spectrum

    %filename and liststring
    filename=['Bruker_' date(1:10)  '_'  num2str(expnb)  timeacq '_fidrefscan.mat'];
    study.filename=filename;
    study.liststring=char(folderexp + study.filename);


    %normalize: xNA/RG/voxelsize 
    study.data.real=study.data.real.*refscan.NArefscan./refscan.RGrefscan./voxvol;
    study.data.imag=study.data.imag.*refscan.NArefscan./refscan.RGrefscan./voxvol;

    %multiplicity
    study.multiplicity=1; %1D

    %process
    study.process.lsfid=round(grpdly)+1;
    study.process.apodizefct='exponential'; %default
    study.process.apodparam1= 0;%default
    study.process.apodparam2=0;%default
    study.process.transfsize=0; %default
    study.process.appltoarray1=0;%default
    study.process.appltoarray2=0;%default
    study.process.phasecorr0=0;%default
    study.process.phasecorr1=0;%default
    study.process.DCoffset=0; %default
    study.process.B0=zeros(1,study.np/2);

    %nt
    study.params.nt=refscan.NArefscan;

    save(folderexp + filename,'study');

    %simple display
    fidrefscan=squeeze(study.data.real)'+1i*squeeze(study.data.imag)'; 
    fidrefscan_sum=sum(fidrefscan,1);
    fidrefscan_sum=[fidrefscan_sum(round(grpdly)+1:end), zeros(1,round(grpdly))];
    ftrefscan=fftshift(fft(fidrefscan_sum./study.params.nt,[],2),2);
    figure; 
    plot(real(ftrefscan),'DisplayName','refscan')
    legend()
    title(num2str(expnb))

end
