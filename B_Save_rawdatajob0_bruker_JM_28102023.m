%% From a PV 360 v1 Bruker exp. folder, adds rawdatajob0 to the study structure from A_
% Jessie Mosso @ EPFL - 28/10/2023
%has been used only for quadrature surface coil - has to be adapted for
%other coil configurations 

clear; close all;

%% DO ADAPT THIS SECTION: Enter exp nb. + directory and run the program 
folder="G:\41-DEMO_DWSPECIAL\";
expnb=18;
timeacq='';%'' or am or pm
dateacq='2024-04-19';
dosave=false;
phaseauto=true; %alternatively - enter manually the phase between the loops "phaseval", for complex summation 
phaseval=90;
nloops=2; 
scalingloops=[100,100]; %can be adapted depending on the desired coil combination

%% DO NOT ADAPT THE FOLLOWING
%% get params: na, nrep, rxarrayphases, 
methodfile = fileread(folder +num2str(expnb)+ "\method");

BRUKERparamnames=["PVM_NRepetitions=","PVM_NAverages=",["##$PVM_ArrayPhase=( 2 )" + sprintf('\n')]];

GUIparamnames=["nrep","nav","rxarrayphases"];
numerical=[true, true, true];
arraylistbol=[false, false, true];

for par=1:length(BRUKERparamnames)
    startind=strfind(methodfile,BRUKERparamnames(par));
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

na=study.nav;
nrep=study.nrep;
if phaseauto
    rxarrayphases=study.rxarrayphases
else
    rxarrayphases=[0,phaseval];
    study.rxarrayphases=rxarrayphases
end 

clear study; 

%% load data rawdata job0

pathrawdata=folder +num2str(expnb)+ "\rawdata.job0";

fileid=fopen(pathrawdata,'r','ieee-le'); %read binary format
if fileid == -1
    disp('Cannot open file');
    return
end

bufferrawdata=fread(fileid,'int32'); 
nbptsfid=length(bufferrawdata)/(2*na*nrep*nloops); %2 is for real, imag
fclose(fileid);

fid_reorganized=zeros(nloops,na*nrep,nbptsfid);

for avrep=1:na*nrep
    for coil=1:nloops
        fid_reorganized(coil,avrep,:)=bufferrawdata((avrep-1)*nbptsfid*2*nloops+(coil-1)*nbptsfid*2+1:2:(avrep-1)*nbptsfid*2*nloops+(coil-1)*nbptsfid*2+nbptsfid*2)+1i*bufferrawdata((avrep-1)*nbptsfid*2*nloops+(coil-1)*nbptsfid*2+2:2:(avrep-1)*nbptsfid*2*nloops+(coil-1)*nbptsfid*2+nbptsfid*2);
    end 
end 

%
fid_comb_re=zeros(na*nrep,nbptsfid);
fid_comb_im=zeros(na*nrep,nbptsfid);

fid_reorganized_rephased_re=zeros(nloops,na*nrep,nbptsfid);
fid_reorganized_rephased_im=zeros(nloops,na*nrep,nbptsfid);

for coil=1:nloops
    fid_comb_re=fid_comb_re + ...
                    (squeeze(real(fid_reorganized(coil,:,:))).*cos(rxarrayphases(coil)/180*pi) ... 
                    -squeeze(imag(fid_reorganized(coil,:,:))).*sin(rxarrayphases(coil)/180*pi)).*scalingloops(coil);

    fid_reorganized_rephased_re(coil,:,:)=(squeeze(real(fid_reorganized(coil,:,:))).*cos(rxarrayphases(coil)/180*pi) ... 
                    -squeeze(imag(fid_reorganized(coil,:,:))).*sin(rxarrayphases(coil)/180*pi)).*scalingloops(coil);

    fid_comb_im= fid_comb_im+ ...
                    (squeeze(real(fid_reorganized(coil,:,:))).*sin(rxarrayphases(coil)/180*pi) ... 
                    +squeeze(imag(fid_reorganized(coil,:,:))).*cos(rxarrayphases(coil)/180*pi)).*scalingloops(coil);

    fid_reorganized_rephased_im(coil,:,:)=(squeeze(real(fid_reorganized(coil,:,:))).*sin(rxarrayphases(coil)/180*pi) ... 
                    +squeeze(imag(fid_reorganized(coil,:,:))).*cos(rxarrayphases(coil)/180*pi)).*scalingloops(coil);
end 


%fid before coil combination but with rephasing
fid_reorganized_rephased=fid_reorganized_rephased_re-1i*fid_reorganized_rephased_im;

%fid after coil combination
fid_comb=fid_comb_re-1i*fid_comb_im; 

grpdly=77;
%fid after coil combination and grpdly shift
fid_comb_grp=[fid_comb(:,grpdly:end),fid_comb(:,1:grpdly-1)];

ft_comb_grp=fftshift(fft(fid_comb_grp,[],2),2);

if dosave

    load([convertStringsToChars(folder) 'Bruker_' dateacq '_' num2str(expnb) timeacq '_ser.mat'])
    
    study.datajob0.real=zeros(size(fid_comb,1),1,size(fid_comb,2));
    study.datajob0.real(:,1,:)=real(fid_comb);
    
    study.datajob0.imag=zeros(size(fid_comb,1),1,size(fid_comb,2));
    study.datajob0.imag(:,1,:)=imag(fid_comb);
    voxvol=study.params.vox1.*study.params.vox2.*study.params.vox3;
    
    study.datajob0.real=study.datajob0.real./study.params.gain./voxvol.*na/(nloops*scalingloops(1));
    study.datajob0.imag=study.datajob0.imag./study.params.gain./voxvol.*na/(nloops*scalingloops(1));
    
    save([convertStringsToChars(folder) 'Bruker_' dateacq '_' num2str(expnb) timeacq '_ser.mat'],'study')
end 

figure; 
plot(real(sum(ft_comb_grp(1:1:end,:))))
hold on 
