%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%                                                               %%%%%%
%%%%%%                       SIGNAL ANALYSIS                         %%%%%%
%%%%%%                                                               %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% NOTE


%% PULIZIA WORKSPACE
clc 
clear all 
close all

%% CARICAMENTO SEGNALI

% load('DataRaw\2021-01-02_15.23.28_ExpGusma_PercorsoCasa_70\ExpGusma_Session1_Number3_Calibrated_SD.mat')
% 
load('DataRaw\BASELINE\2020-12-07_22.49.42_ExpGusma_MultiSession\ExpGusma_Session3_Number3_Calibrated_SD.mat')
accDataRaw =[Number3_Accel_LN_X_CAL , Number3_Accel_LN_Y_CAL ,Number3_Accel_LN_Z_CAL];
gyroDataRaw = [ Number3_Gyro_X_CAL , Number3_Gyro_Y_CAL , Number3_Gyro_Z_CAL] ;
magDataRaw = [Number3_Mag_X_CAL  Number3_Mag_Y_CAL  Number3_Mag_Z_CAL];



%% TIME LINE

%frequenza & linespace
N = size(accDataRaw,1);
Fs = 51.2;
Ts=1/Fs;
t = 0:Ts:((N-1)*(1/Fs));  


%% PRINT DATARAW 

printSignal3D(gyroDataRaw(:,1)  ,  gyroDataRaw(:,2)  , gyroDataRaw(:,3)   ,1:N,'Gyro Data Raw','sample','Gyro(deg/s^2)',1);

printSignal3D(accDataRaw(:,1)  ,  accDataRaw(:,2) ,  accDataRaw(:,3)  ,1:N,'Acceleration Data RAw','Sample','Accelaration(m/s^2)',2);
  
printSignal3D(magDataRaw(:,1)  ,  magDataRaw(:,2) ,  magDataRaw(:,3)  ,1:N,'magnetometro Data Raw ' , 'Sample' , 'Cambio campo magnetico ( boh )',3);
  
%% COSTANT BIAS Acc Gyro

meanBiasAccx=mean(accDataRaw(:,1));
meanBiasAccx
meanBiasAccy=mean(accDataRaw(:,2));
meanBiasAccy
meanBiasAccz=mean(accDataRaw(:,3));
meanBiasAccz


meanBiasGyrox=mean(gyroDataRaw(:,1));
meanBiasGyrox
meanBiasGyroy=mean(gyroDataRaw(:,2));
meanBiasGyroy
meanBiasGyroz=mean(gyroDataRaw(:,3));
meanBiasGyroz


%% velocity random walk

accNavFrame=zeros(N,3);
% accDataRaw(:,3)=accDataRaw(:,3)-9.8.*ones(N,1);
for i=1:N
    Rrpy=rotationMatrix(gyroDataRaw(i,1),gyroDataRaw(i,2),gyroDataRaw(i,3));
    accNavFrame(i,:)=Rrpy*accDataRaw(i,:)';
end
vNavFrame=zeros(N,3);
for i=2:N
    vNavFrame(i,:)=vNavFrame(i-1,:)+accNavFrame(i,:).*Ts;
end


printSignal3D(vNavFrame(:,1)  ,  vNavFrame(:,2)  , vNavFrame(:,3)   ,1:N,'Velocity Random Walk','sample','Velocity(m/s)',5);

%% HISTOGRAM ANALISYS (FREQUENCY DISTRIBUTION)


    figure(4)
    subplot(2,3,1);
    histogramAnalysis(accDataRaw(:,1),'accXData histogram','Frequency(Hz)','Magnitude');

    subplot(2,3,2);
    histogramAnalysis(accDataRaw(:,2),'accYData histogram','Frequency(Hz)','Magnitude');

    subplot(2,3,3);
    histogramAnalysis(accDataRaw(:,3),'accZData histogram','Frequency(Hz)','Magnitude');

    subplot(2,3,4);
    histogramAnalysis(gyroDataRaw(:,1),'gyroXData histogram','Frequency(Hz)','Magnitude');

    subplot(2,3,5);
    histogramAnalysis(gyroDataRaw(:,2),'gyroYData histogram','Frequency(Hz)','Magnitude');

    subplot(2,3,6);
    histogramAnalysis(gyroDataRaw(:,3),'gyroZData histogram','Frequency(Hz)','Magnitude');
%     
%     subplot(3,3,7);
%     histogramAnalysis(magDataRaw(:,1),'Magnetometer_X_Data histogram','Frequency(Hz)','Magnitude');
% 
%     subplot(3,3,8);
%     histogramAnalysis(magDataRaw(:,2),'Magnetometer_Y_Data histogram','Frequency(Hz)','Magnitude');
% 
%     subplot(3,3,9);
%     histogramAnalysis(magDataRaw(:,3),'Magnetometer_Z_Data histogram','Frequency(Hz)','Magnitude');



%% 

% Welch's method, named after Peter D. Welch, is an approach for spectral density estimation. It is used in physics, 
% engineering, and applied mathematics for estimating the power of a signal at different frequencies. 
% 
% signal=[variance ,a];
% 
% figure(6)
% pwelch(signal,[],[],[],Fs,'centered','power')
% legend('variance','magnitude accelaration')
% hold on
% pwelch(a)

% figure(3)
% [pxx,f] = pwelch(signal,120,50,120,Fs,'centered','power');
% plot(f,10*log10(pxx))
% xlabel('Frequency (Hz)')
% ylabel('Magnitude (dB)')
% grid

% 
% figure(4)
% [pxx,f] = pwelch(signal(:,1),[],[],[],Fs,'centered','power');
% pmax = pwelch(signal(:,1),[],[],[],Fs,'maxhold','centered','power');
% pmin = pwelch(signal(:,1),[],[],[],Fs,'minhold','centered','power');
% 
% plot(f,pow2db(pxx))
% hold on
% plot(f,pow2db([pmax pmin]),':')
% hold off
% xlabel('Frequency (kHz)')
% ylabel('Power (dB)')
% legend('pwelch','maxhold','minhold')
% title('Centered Power Spectrum Estimates')
% grid


% figure(5)
% [pxx,f] = pwelch(signal(:,2),[],[],[],Fs,'centered','power');
% pmax = pwelch(signal(:,2),[],[],[],Fs,'maxhold','centered','power');
% pmin = pwelch(signal(:,2),[],[],[],Fs,'minhold','centered','power');
% 
% plot(f,pow2db(pxx))
% hold on
% plot(f,pow2db([pmax pmin]),':')
% hold off
% xlabel('Frequency (kHz)')
% ylabel('Power (dB)')
% legend('pwelch','maxhold','minhold')
% title('Centered Power Spectrum Estimates')
% grid



% signal è il*segnale da analizzare.
% fs è il*numero*di*elementi*su*cui*calcolare*gli*spettri*da*mediare,*ovvero*la*
% durata della* finestra* temporale* da* considerare.* Inserire* il* valore* della*
% frequenza*di*campionamento*fs*equivale*a*impostare*una*finestra*di*1s.*
% [* ]: questo* simbolo* si* usa* per* l’omissione* dei* parametri,* che* significa*
% accettare* le* impostazioni* di* default.* In* particolare* questo* parametro*
% riguarda* il* numero* di* campioni* che* vengono* sovrapposti* tra* le* varie*
% finestre.*Di*default*la*sovrapposizione*è*fissata*al*50%.
% nfft è*il* numero* di* campioni* del* segnale* su* cui* calcolare*la* trasformata* di*
% Fourier* per* la* stima* dello* spettro* della* finestra* temporale.* Il* valore* di*
% default* è* pari* alla* potenza* di* 2* successiva* alla* lunghezza* della* finestra,*
% questo*serve*per*migliorare*il*costo*computazionale*dell’algoritmo*fft*(Fast*
% Fourier* Transform* è* infatti* ottimizzata* se* il* numero* di* elementi* da*
% analizzare*è*una*potenza*di*2).*Gli*elementi*che*eccedono*la*lunghezza*della*
% finestra*vengono*utilizzati*per*lo*zero*padding. E’*stato*scelto*di*modificare*
% questa*impostazione* di* default* perché,*essendo*fs=1000,* nfft* sarebbe* stato*
% 1024*e*lo*zero*padding*non*sarebbe*stato*efficace*con*soli*24*elementi.*Alla*
% riga* 6* del* codice* è* stato* quindi* impostato* nfft* pari* a* due* potenze* di* 2*
% successive* al* valore* di* fs* (nfft=2048),* utilizzando* la* funzione* di* MATLAB
% nextpow2.m,*che*arrotonda*al*numero*intero*successivo*il*logaritmo*in*base*
% 2*del*numero*dato*in*ingresso.
% fs è*la*frequenza*di*campionamento*del*segnale*da*analizzare*(xA).
% ‘onesided’ significa*che*xA*in*ingresso*è*un*segnale*reale.
% 
% 
% figure(7)
% nfft=2^(nextpow2(Fs)+1);
% [pxx,f] = pwelch(signal(:,1),Fs,[],nfft,Fs,'centered','power');
% pmax = pwelch(signal(:,1),Fs,[],nfft,Fs,'maxhold','centered','power');
% pmin = pwelch(signal(:,1),Fs,[],nfft,Fs,'minhold','centered','power');
% 
% plot(f,pow2db(pxx))
% hold on
% plot(f,pow2db([pmax pmin]),':')
% hold off
% xlabel('Frequency (kHz)')
% ylabel('Power (dB)')
% legend('pwelch','maxhold','minhold')
% title('Centered Power Spectrum Estimates')
% grid
% 
% 
% figure(8)
% nfft=2^(nextpow2(Fs)+1);
% [pxx,f] = pwelch(signal(:,2),Fs,[],nfft,Fs,'centered','power');
% pmax = pwelch(signal(:,2),Fs,[],nfft,Fs,'maxhold','centered','power');
% pmin = pwelch(signal(:,2),Fs,[],nfft,Fs,'minhold','centered','power');
% 
% plot(f,pow2db(pxx))
% hold on
% plot(f,pow2db([pmax pmin]),':')
% hold off
% xlabel('Frequency (kHz)')
% ylabel('Power (dB)')
% legend('pwelch','maxhold','minhold')
% title('Centered Power Spectrum Estimates')
% grid
% 
%    
%% PRINT METHOD 

function printSignal(signal,t,Title,Xax,Yax,indexFigure)
    figure(indexFigure)
    plot(t, signal)
    xlabel(Xax)
    ylabel(Yax)
    title(Title)
end

function printSignalSub(signal,t,Title,Xax,Yax)
    plot(t, signal)
    xlabel(Xax)
    ylabel(Yax)
    title(Title)
end
function printSignal3D(signalX,signalY,signalZ,t,Title,Xax,Yax,indexFigure)
    figure(indexFigure)
    plot(t, signalX)
    hold on
    plot(t,signalY)
    hold on
    plot(t,signalZ)
    xlabel(Xax)
    ylabel(Yax)
    title(Title)
    legend('Xax','Yax','Zax')
end

function printQuaternions(signalW,signalX,signalY,signalZ,t,Title,Xax,Yax,indexFigure)
    figure(indexFigure)
    plot(t,signalW)
    hold on
    plot(t, signalX)
    hold on
    plot(t,signalY)
    hold on
    plot(t,signalZ)
    xlabel(Xax)
    ylabel(Yax)
    title(Title)
    legend('Xax','Yax','Zax','Wax')
end


%% FUNCTIONS

function histogramAnalysis(data,Title,Xax,Yax)
    FData=fft(data);
    assex=linspace(-1/2,1/2,numel(data));
    printSignalSub(abs(fftshift(FData))/numel(FData),assex,Title,Xax,Yax);
end

function [Rrpy]= rotationMatrix(phi,theta,psi)
    cphi=cos(deg2rad(phi));
    ctheta=cos(deg2rad(theta));
    cpsi=cos(deg2rad(psi));
    sphi=sin(deg2rad(phi));
    stheta=sin(deg2rad(theta));
    spsi=sin(deg2rad(psi));
    
    Rrpy=[ ctheta*cphi   stheta*spsi*cphi-sphi*cpsi      stheta*cphi*cpsi + spsi*sphi 
           ctheta*sphi    stheta*sphi*spsi + cpsi*cphi   stheta*cpsi*sphi-cphi*spsi 
           -stheta       ctheta*spsi                       ctheta*cpsi];
end

function [accDataRaw , gyroDataRaw ,quaternions]= acquisizione_dati(WalkInertial,soggetto,stringa,i)
    
    %3-axes acceleration
    AccX=getfield(WalkInertial(soggetto),'AccX',stringa);
    AccY=getfield(WalkInertial(soggetto),'AccY',stringa);
    AccZ=getfield(WalkInertial(soggetto),'AccZ',stringa);
    %The unit of the calibrated accelerometer is meters per square second(m/s2)
    accXDataRaw=cell2mat(AccX(i));
    accYDataRaw=cell2mat(AccY(i));
    accZDataRaw=cell2mat(AccZ(i));
    accDataRaw=[accXDataRaw(1:23821) accYDataRaw accZDataRaw];
    %primo soggetto accex 1:6180
    %secondo soggetto accex 1:6183
    %3-axis gyroscope
    GyroX=getfield(WalkInertial(soggetto),'GyroX',stringa);
    GyroY=getfield(WalkInertial(soggetto),'GyroY',stringa);
    GyroZ=getfield(WalkInertial(soggetto),'GyroZ',stringa);
    %The unit of the calibrated accelerometer is meters per square second(m/s2)
    gyroXDataRaw=cell2mat(GyroX(i));
    gyroYDataRaw=cell2mat(GyroY(i));
    gyroZDataRaw=cell2mat(GyroZ(i));
    gyroDataRaw=[gyroXDataRaw gyroYDataRaw gyroZDataRaw];
    

    %quaternion obatained from angvel
    QuatW=getfield(WalkInertial(soggetto),'QuatX',stringa);
    QuatX=getfield(WalkInertial(soggetto),'QuatX',stringa);
    QuatY=getfield(WalkInertial(soggetto),'QuatX',stringa);
    QuatZ=getfield(WalkInertial(soggetto),'QuatX',stringa);
    QuatW=cell2mat(QuatW(i));
    QuatX=cell2mat(QuatX(i));
    QuatY=cell2mat(QuatY(i));
    QuatZ=cell2mat(QuatZ(i));
    %quaternion
    quaternions=quaternion(QuatW  , QuatX , QuatY  ,  QuatZ );
end

function [data]=lowPassAndSmothing(data)
    %provare a tagliare a 15 hz per movimento del braccio e a grado 10 

    %low-pass filtering
    filteredgyroX=butterworth(data(:,1),6,0.3,'low');
    filteredgyroY=butterworth(data(:,2),6,0.3,'low');
    filteredgyroZ=butterworth(data(:,3),6,0.3,'low');
    %smoothing signal
    data(:,1)=movmean(filteredgyroX,70);
    data(:,2)=movmean(filteredgyroY,70);
    data(:,3)=movmean(filteredgyroZ,70);
end


function [newData] = butterworth(dataRaw , grade , cutoff,type)
    [a,b]=butter(grade,cutoff,type);
    %figure(10),freqz(a,b)
    newData=filtfilt(a,b,dataRaw);
end

function [ roll , pitch , yaw]= inverseSolution(rotMatrix)
    pitch=atan2(-rotMatrix(3,1) , sqrt(rotMatrix(3,2)^2 +rotMatrix(3,3)^2));
    if -pi/2 < pitch && pitch < pi/2
        yaw=atan2(rotMatrix(3,2),rotMatrix(3,3));
        pitch=atan2(-rotMatrix(3,1) , sqrt(rotMatrix(3,2)^2 +rotMatrix(3,3)^2));
        roll=atan2(rotMatrix(2,1),rotMatrix(1,1));
    else if pi/2 < pitch && pitch <3*pi/2
            yaw=atan2(-rotMatrix(3,2),-rotMatrix(3,3));
            pitch=atan2(-rotMatrix(3,1),-sqrt(rotMatrix(3,2)^2 + rotMatrix(3,3)^2));
            roll=atan2(-rotMatrix(2,1),-rotMatrix(1,1));
        end
    end
end