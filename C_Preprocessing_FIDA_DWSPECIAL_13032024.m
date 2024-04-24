%% Preprocessing of data after reading with A_ and B_
% Jessie Mosso @ EPFL - 13/03/2024
% Requires FID-A from Jamie Near  
% https://github.com/CIC-methods/FID-A
clear;clc;close all; 

%% add FID-A master to Matlab path
%uses custom FID-A functions found in custom_FID-A_functions


%% DO ADAPT THIS SECTION: 
LBall=10; %LB applied for processing only
doprocessing=true;
dosavesum=true;
dosaveprocessing=true; 
spectralreginf=0.5; %ppm %limits for spectral registration
spectralregsup=4.3; %ppm
filelist = dir([pwd '/raw/*ser.mat']);


%% DO NOT ADAPT THE FOLLOWING
%% step 1 - do preprocessing with FIDA
if doprocessing
for file=1:length(filelist)
    disp(filelist(file).name)
    f1=figure;
    
    %% 1-Convert Bruker study to FID A structure
    out0a=convertNicoStudyFidAformat_isison([pwd '/raw'], filelist(file).name);
    out0b=convertNicoStudyFidAformat_isisoff([pwd '/raw'], filelist(file).name);
    %apply LB
    dw=out0a.dwelltime;
    tt=[0:dw:dw*(out0a.n-1)];
    out0alb=out0a;
    fids0alb=out0alb.fids.*repmat(exp(-tt*pi*LBall).',1,out0a.averages);
    out0alb.fids=fids0alb;
    out0alb.specs=fftshift(fft(out0alb.fids.',[],2),2).';
    out0blb=out0b;
    fids0blb=out0blb.fids.*repmat(exp(-tt*pi*LBall).',1,out0b.averages);
    out0blb.fids=fids0blb;
    out0blb.specs=fftshift(fft(out0blb.fids.',[],2),2).';
    
    figure(f1);
    subplot(2,2,1)
    for k=1:out0a.averages
        plot(real(out0a.specs(:,k)))
        hold on 
        xlim([2900,3050])
    end
    for k=1:out0b.averages
        plot(real(out0b.specs(:,k)))
        hold on 
        xlim([2900,3050])
    end 
    title(['0-ini'])
    
    figure(f1);
    subplot(2,2,2)
    for k=1:out0alb.averages
        plot(real(out0alb.specs(:,k)))
        hold on 
        xlim([2900,3050])
    end 
    for k=1:out0blb.averages
        plot(real(out0blb.specs(:,k)))
        hold on 
        xlim([2900,3050])
    end 
    title(['1-ini - LB' num2str(LBall)])
    
    %% 2-align av.
    [out1alb,fsa,phsa]=op_alignAverages_fd_jm_conj(out0alb,4.7+spectralreginf,4.7+spectralregsup,0.5,'y');%5.2,9,0.5,'y');
    [out1blb,fsb,phsb]=op_alignAverages_fd_jm_conj(out0blb,4.7+spectralreginf,4.7+spectralregsup,0.5,'y');%5.2,9,0.5,'y');
    % fs        = Vector of frequency shifts (in Hz) used for alignment.
    % phs       = Vector of phase shifts (in degrees) used for alignment.

    figure(f1);
    subplot(2,2,3)
    for k=1:out1alb.averages
        plot(real(out1alb.specs(:,k)))
        hold on 
        xlim([2900,3050])
    end
    for k=1:out1blb.averages
        plot(real(out1blb.specs(:,k)))
        hold on 
        xlim([2900,3050])
    end
    title(['2-spec. reg. JNear - LB' num2str(LBall)])
    
    %remove LB
    out1a=out1alb;
    fids1a=out1alb.fids.*repmat(exp(tt*pi*LBall).',1,out1alb.averages);
    out1a.fids=fids1a;
    out1a.specs=fftshift(fft(out1a.fids.',[],2),2).';
    out1b=out1blb;
    fids1b=out1blb.fids.*repmat(exp(tt*pi*LBall).',1,out1blb.averages);
    out1b.fids=fids1b;
    out1b.specs=fftshift(fft(out1b.fids.',[],2),2).'; 

    %% 3-outlier removal 
    [out2a,metrica,badAveragesa]=op_rmbadaverages_jm(out1a,1.5,'f'); %performs 10Hz LB inside
    [out2b,metricb,badAveragesb]=op_rmbadaverages_jm(out1b,1.5,'f'); %performs 10Hz LB inside
    
    %apply LB
    out2alb=out2a;
    fids2alb=out2a.fids.*repmat(exp(-tt*pi*LBall).',1,out2a.averages);
    out2alb.fids=fids2alb;
    out2alb.specs=fftshift(fft(out2alb.fids.',[],2),2).';
    out2blb=out2b;
    fids2blb=out2b.fids.*repmat(exp(-tt*pi*LBall).',1,out2b.averages);
    out2blb.fids=fids2blb;
    out2blb.specs=fftshift(fft(out2blb.fids.',[],2),2).';
    
    figure(f1);
    subplot(2,2,4)
    for k=1:out2alb.averages
        plot(real(out2alb.specs(:,k)))
        hold on 
        xlim([2900,3050])
    end
    for k=1:out2blb.averages
        plot(real(out2blb.specs(:,k)))
        hold on 
        xlim([2900,3050])
    end
    title(['3-rmv motion corrupted av. - LB' num2str(LBall)])
    title(filelist(file).name)

    %combine on/off
    clear fidtot;
    fida=out1a.fids.'; 
    fidb=out1b.fids.'; 
    fidtot(1:2:size(fida,1)*2,:)=fida; 
    fidtot(2:2:size(fida,1)*2,:)=fidb;
    
    ind=1;
    clear fid2sum;
    for k=1:size(fidtot,1)/2
        fid2sum(ind,:)=sum(fidtot((k-1)*2+1:k*2,:)); 
        ind=ind+1;
    end 
    
    ind=1;
    clear fidmocor;
    for k=1:size(fidtot,1)/2 %from 1 to 80 pairs
        if ismember(k,badAveragesa) % if in pair k, the odd is an outlier (= present in the 1st outlier list) - out
        else % if in pair k, the odd is an not an outlier 
            if ismember(k,badAveragesb) % if in pair k, the even is an outlier (= present in the 2nd outlier list) - out
            else % if in pair k, neither the odd nor the even are outliers
                fidmocor(ind,:)=fid2sum(k,:);
                ind=ind+1; 
            end 
        end 
    end 
    fidmocor=conj(fidmocor);
    
    %% 4-add all the info to the Matlab study structure
    if dosaveprocessing
        load([pwd '/raw/' filelist(file).name]);
        study.fidaprocess.phsa=phsa;
        study.fidaprocess.fsa=fsa;
        study.fidaprocess.metrica=metrica;
        study.fidaprocess.badAveragesa=badAveragesa; 
        study.fidaprocess.phsb=phsb;
        study.fidaprocess.fsb=fsb;
        study.fidaprocess.metricb=metricb;
        study.fidaprocess.badAveragesb=badAveragesb;

        study.params.nt=size(fidmocor,1)*2; 
        study.multiplicity=size(fidmocor,1);

        study.process.apodparam1=zeros(1,size(fidmocor,1));
        study.process.apodparam2=zeros(1,size(fidmocor,1));
        study.process.phasecorr0=zeros(1,size(fidmocor,1));
        study.process.phasecorr1=zeros(1,size(fidmocor,1));

        study.data.real=zeros(size(fidmocor,1),1,size(fidmocor,2));
        study.data.real(:,1,:)=real(fidmocor);
        study.data.imag=zeros(size(fidmocor,1),1,size(fidmocor,2));
        study.data.imag(:,1,:)=imag(fidmocor);

        if ~exist([pwd '/processed/'], 'dir')
           mkdir([pwd '/processed/'])
        end
        save([pwd '/processed/' filelist(file).name(1:end-4) '_processed.mat'],'study')   
    end 
    end 
end 
    
%% STEP 2 - APPLY PREPROCESSING, SUM AND PHASE THE SUM 

if dosavesum
for file=1:length(filelist)
    load([pwd '/processed/' filelist(file).name(1:end-4) '_processed.mat'],'study')   

    fidmocor=squeeze(study.data.real)+1i*squeeze(study.data.imag);
    sumfid=sum(fidmocor)./(size(fidmocor,1).*2); %normalization after outlier removal - should be done either here or at quantification
   
    %% save
    study.data.real=zeros(1,1,study.np/2);
    study.data.imag=zeros(1,1,study.np/2);

    study.data.real(1,1,:)=real(sumfid); 
    study.data.imag(1,1,:)=imag(sumfid);

    study.multiplicity=1;
    study.process.lsfid=0;
    study.process.apodparam1=0;
    study.process.apodparam2=0;
    study.process.phasecorr0=0;
    study.process.phasecorr1=0;
    study.process.B0=zeros(1,study.np/2);

    filename=filelist(file).name(1:end-4);
    if ~exist([pwd '/processed/sum/'], 'dir')
        mkdir([pwd '/processed/sum/'])
    end
    save([pwd '/processed/sum/SUM_' filename '_processed.mat'],'study');
    
    %plot
    figure; 
    ftcorr=fftshift(fft(sumfid./study.params.nt,[],2),2); 
    plot(real(ftcorr))
    hold on 
    load([pwd '/raw/' filelist(file).name(1:end-4) '.mat'],'study')   
    fidini=squeeze(study.data.real)+1i*squeeze(study.data.imag);
    sumfidini=sum(fidini)./study.params.nt;
    sumfidini=[sumfidini(77:end),zeros(1,76)]; 
    ftini=fftshift(fft(sumfidini,[],2),2); 
    plot(real(ftini))
    legend('corr','ini')

end
end 
