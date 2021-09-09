%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%                                                               %%%%%%
%%%%%%                       StepDetection                          %%%%%%
%%%%%%                                                               %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% NOTE
% IN QUESTA VERSIONE UTILIZZEREMO LA TECNICA DI INTEGRAZIONE DI EULERO DI 

function [s]=stepDetection(accDataRaw,gyroDataRaw)


%% INIZIALIZZAZIONI




%frequenza & linespace
N = size(accDataRaw,1);
Fs = 51.2;
Ts=1/Fs;
t = 0:Ts:((N-1)*(1/Fs)); 
sampleLine=1:N;


%%

printSignal3D(accDataRaw(:,1)  ,  accDataRaw(:,2) ,  accDataRaw(:,3)  ,1:N,'Acceleration Data RAw','Sample','Accelaration(m/s^2)',1);
  
printSignal(accDataRaw(:,3),1:N,'acc y','sample','ms2',8)
% printSignal3D(gyroDataRaw(:,1)  ,  gyroDataRaw(:,2)  , gyroDataRaw(:,3)   ,1:N,'Gyro Data Raw','sample','Gyro(deg/s^2)',1);

%% SIGNAL FILTERING
%VALUTARE SE TENERE O MENO SOLO PER GYRO RAW

% % %low-pass filtering
accDataRaw(:,1)=butterworth(accDataRaw(:,1),2,0.4,'low');
accDataRaw(:,2)=butterworth(accDataRaw(:,2),2,0.4,'low');
accDataRaw(:,3)=butterworth(accDataRaw(:,3),2,0.4,'low');
% 
%smoothing signal
accDataRaw(:,1)=movmean(accDataRaw(:,1),10);
accDataRaw(:,2)=movmean(accDataRaw(:,2),10);
accDataRaw(:,3)=movmean(accDataRaw(:,3),10);
 

%% COMPUTE MAGNITUDE ACCELERATION 

a=sqrt(accDataRaw(:,1).^2+accDataRaw(:,2).^2+accDataRaw(:,3).^2);

%% COMPUTE LOCAL ACCELERATION VARIANCE

variance=movvar(a,20);
printSignal(a, 1:N , ' varianza dell accelerazione' , 'Samples' ,'m/s^2' , 2  )

%% FIND PEAKS

ppm=70;
tempoPasso=60/ppm;
tempoStep=tempoPasso*2;
samplesStep=(tempoStep-0.3)*128;

MinPeakHeight=6;
MinPeakDistance=samplesStep;
%60=MinPeakDistance/2;

figure(3)
findpeaks(variance,sampleLine,'MinPeakDistance',60,'MinPeakHeight',6)
[peaks,locs,w,p]=findpeaks(variance,sampleLine,'MinPeakDistance',60,'MinPeakHeight',6);

% col tempo
% findpeaks(variance,Fs,'MinPeakDistance',1.30,'MinPeakHeight',6)
% [peaks,locs,w,p]=findpeaks(variance,sampleLine,'MinPeakDistance',166.4,'MinPeakHeight',6);

%% COMPUTE THRESHOLD VARIABLES

diffPeaks=diff(locs);
meandiffPeaks=mean(diffPeaks);
meanPeaks=mean(peaks);
if meandiffPeaks < 81
    T1=5;
    T2=4;
else
    T1=2;
    T2=0.35;
end


%% COMPUTE THRESHOLD LINE

% threshold 1 si pùò considerare PROVVISSORIA STANCE FASE almeno quando la varianza scende sotto T1 (restandoci per una finestra limitata)
% threshold 2 si considera DEFINITIVAMENTE STANCE FASE quando in unintervallo predefinito la varianza scende sotto T2
% si può considerare PROVVISSORIO SWING FASE quando B1 sale rispetto il sample precedente (rimanendoci per una finestra)e quando il massimo di B2 è 0 in una finestra

B1=zeros(N,1);
B2=zeros(N,1);

for i=1:N
    if(variance(i)>T1)
        B1(i)=T1;
    end
    if(variance(i)<T2)
        B2(i)=T2;
    end
end


%% STEP SAMPLE DETECT



startStanceFase=[1 0];
endStanceFase=[];
for i=2:(N-100)
    %stance fase
    if (B1(i-1)>B1(i))  && (max(B1(i:i+(MinPeakDistance*1/30)))== 0) && (max(B2(i:i+(MinPeakDistance*1/7)))== T2)  && (i-(endStanceFase(end,1))> (MinPeakDistance*1/30)) && size(startStanceFase,1) == size(endStanceFase,1)
               startStanceFase=[startStanceFase ; i,0];
    % swing phase
    
    elseif (B1(i-1)<B1(i)) && (min(B2(i:i+(MinPeakDistance*1/30)))==0) && (i-(startStanceFase(end,1)) > (MinPeakDistance*1/30) && size(startStanceFase,1)-1==size(endStanceFase,1) )
            
                endStanceFase=[endStanceFase ;i, 0];

    
    end
end

%% STEP FREQUENCY

samplesF=diff(startStanceFase(2:length(startStanceFase)-1,1));
media=mean(samplesF);
strideDuration=media/Fs;
stepDuration=strideDuration/2;
stepFrequency=60/stepDuration;

%% STEP DETECTION PLOTTING

% printStepDetection(a,variance,B1,B2,startStanceFase,endStanceFase,1:N,'Step Detection analysis','Time(s)', 'acceleration(m/s^2)',2);

printStepDetection2(variance,B1,B2,endStanceFase,startStanceFase,1:N,'Step Detection analysis','Sample', 'acceleration(m/s^2)',4);

%% ANGULAR DISPLACEMENT

vAng=gyroDataRaw;
angDi=zeros(N,3);
for i=2:N
    angDi(i,:)=angDi(i-1,:)+vAng(i,:).*Ts;
end

% % %low-pass filtering
 
% 
%smoothing signal
%  angDi(:,1)=movmean(angDi(:,1),50);
%  angDi(:,2)=movmean(angDi(:,2),50);
%  angDi(:,3)=movmean(angDi(:,3),50);

angDi=movmean(angDi,128);

printSignal3D(angDi(:,1)  ,  angDi(:,2) ,  angDi(:,3)  ,1:N,'Angular Displacement','Sample','Degree(deg)',5);

%% REFERENCE FRAME TRASFORMATION FROM SENSOR TO GLOBAL  RPY

accNavFrame=zeros(N,3);
% accDataRaw(:,3)=accDataRaw(:,3)-9.8.*ones(N,1);
for i=1:N
    Rrpy=rotationMatrix(angDi(i,1),angDi(i,2),angDi(i,3));
    accNavFrame(i,:)=Rrpy*accDataRaw(i,:)';
end

%%   VELOCITY repesct to NavFrame WITH ZUPT
%la velocità per ogni passo non va accumulata , perchè ognni passo ha una
%velocità a sè


vNavFrame=zeros(startStanceFase(end,1),3);
vBias=zeros(startStanceFase(end,1),3);
k=1;
for i=2:startStanceFase(end,1)
    if i >= endStanceFase(k,1) && i < startStanceFase(k+1,1)
            vNavFrame(i,:)=vNavFrame(i-1,:)+ accNavFrame(i,:).*Ts*1/3;
        elseif i == startStanceFase(k+1,1)
            velocity(k)=sqrt( mean(vNavFrame(endStanceFase(k,1):startStanceFase(k+1,1),1)).^2  +  mean(vNavFrame(endStanceFase(k,1):startStanceFase(k+1,1),2)).^2) ;
            k=k+1;
    end
end  

figure(8)
plot( velocity)
title('velocity')
% hold on
% legend('x','y','z')
strideVelocity=mean(velocity);
%% DISPLACEMENT

% lo spostamento non va accumulato tra i diversi passi , ogni passo ha il suo
% spostamento. Lo spostamento per ogni campione del passo va accumulata.

displacement=zeros(startStanceFase(end,1),3);
desplacementAtStep=zeros(size(startStanceFase,1),3);
k=1;
for i=2:startStanceFase(end,1)
    if i >= endStanceFase(k,1) && i < startStanceFase(k+1,1)
        displacement(i,:)=displacement(i-1,:)+vNavFrame(i,:).*Ts;
    elseif i == startStanceFase(k+1,1)
            desplacementAtStep(k,:)=displacement(i-1,:);
            k=k+1;
            %displacement(i,:)=desplacementAtStep(k-1);
    else 
        %displacement(i,:)=desplacementAtStep(k-1);
    end
end   


%% STRIDE LENGTH

strideLength=zeros(size(startStanceFase,1)-1,1);
for i=2:size(startStanceFase,1)-2
        strideLength(i)=sqrt(desplacementAtStep(i,1).^2+desplacementAtStep(i,2).^2);
end

figure(9)
plot( strideLength)
title('strideLength')
%% POSITION
position=zeros(size(startStanceFase,1)-1,size(startStanceFase,1)-1);
for i=2:size(startStanceFase,1)-2
        position(i,1)=position(i-1,1)+cos(deg2rad(angDi(startStanceFase(i,1),3)))*strideLength(i);
        position(i,2)=position(i-1,2)+sin(deg2rad(angDi(startStanceFase(i,1),3)))*strideLength(i);
end
position=position(1:end-1,:);
figure(6)
plot(position(:,1),position(:,2))
title('tracking')

%% SHOW RESULTS

disp('stride lentgh (m)')
mean(strideLength)

% meandiffPeaks
% strideVelocity


%% CREATE STRUCTURE

s=struct('strideVelocity',strideVelocity,'stepFrequency',stepFrequency,'stepDuration',stepDuration,'strideLength',mean(strideLength));

%% PRINT METHOD 

function printSignal(signal,t,Title,Xax,Yax,indexFigure)
    figure(indexFigure)
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

function printStepDetection(signal1,signal2,signal3,signal4,endStanceFase,startStanceFase,t,Title,Xax,Yax,indexFigure)
    figure(indexFigure)
    plot(t, signal1)
    hold on
    plot(t,signal2)
    hold on
    plot(t,signal3)
    hold on
    plot(t,signal4)
    hold on 
    plot(startStanceFase(:,1),startStanceFase(:,2),'o');
    hold on
    plot(endStanceFase(:,1),endStanceFase(:,2),'*');
    xlabel(Xax)
    ylabel(Yax)
    title(Title)
    legend('Magnitude','Variance','startStanceFase','endStanceFase')
end

function printStepDetection2(signal2,signal3,signal4,endStanceFase,startStanceFase,t,Title,Xax,Yax,indexFigure)
    figure(indexFigure)
    plot(t,signal2)
    hold on
    plot(t,signal3)
    hold on
    plot(t,signal4)
    hold on 
    plot(startStanceFase(:,1),startStanceFase(:,2),'o');
    hold on
    plot(endStanceFase(:,1),endStanceFase(:,2),'*');
    xlabel(Xax)
    ylabel(Yax)
    title(Title)
    legend('Variance','startStanceFase','endStanceFase')
end


function printStep(signal1,signal2,signal3,signal4,signal5,endStanceFase,startStanceFase,t,Title,Xax,Yax,indexFigure)
    figure(indexFigure)
    plot(t, signal1)
    hold on
    plot(t,signal2)
    hold on
    plot(t,signal3)
    hold on
    plot(t,signal4)
    hold on 
    plot(t,signal5)
    hold on 
    plot(startStanceFase(:,1),startStanceFase(:,2),'o');
    hold on
    plot(endStanceFase(:,1),endStanceFase(:,2),'*');
    xlabel(Xax)
    ylabel(Yax)
    title(Title)
    legend('Magnitude','Variance')
end

%% FUNCTIONS

function [accDataRaw , gyroDataRaw ]= acquisizione_dati(WalkInertial,soggetto,stringa,i)
    
    limite= 15000;
    %3-axes acceleration
    AccX=getfield(WalkInertial(soggetto),'AccX',stringa);
    AccY=getfield(WalkInertial(soggetto),'AccY',stringa);
    AccZ=getfield(WalkInertial(soggetto),'AccZ',stringa);
    %The unit of the calibrated accelerometer is meters per square second(m/s2)
    accXDataRaw=cell2mat(AccX(i));
    accYDataRaw=cell2mat(AccY(i));
    accZDataRaw=cell2mat(AccZ(i));
    accDataRaw=[accXDataRaw(1: limite ,1) accYDataRaw(1: limite ,1) accZDataRaw(1: limite ,1)];
     %primo soggetto accex 24000   21000  20000
     %esperimento non valido 1  
    %terzo soggetto accex 23700 23000   23100
    %esperimento 3 non valido
    %quinto soggetto acce 19000  19200   19250
    %esperimento 2 ha una frequenza alta
    %ottavo soggetto acc 22450  22600 29000
    %terzo e primo esperimento sbagliatissimo secondo almeno funzionano i segnali
    %nono soggetto  17950  18950 19950
    % 1 esp accellerazione nulla su yz 2 esp x risulta nulla 
    %3 esp nulla yz fino a che si annula anche x
    %undicesimo sogg 18400
    % tredicesimo soggetto 15000
    %quindicesimo soggetto  17000 si salva il terzo esperimento
    %diciasettesimo soggetto 22000 non si salva niente
    %diciannovesimo soggetto 15000
    %3-axis gyroscope
    GyroX=getfield(WalkInertial(soggetto),'GyroX',stringa);
    GyroY=getfield(WalkInertial(soggetto),'GyroY',stringa);
    GyroZ=getfield(WalkInertial(soggetto),'GyroZ',stringa);
    %The unit of the calibrated accelerometer is meters per square second(m/s2)
    gyroXDataRaw=cell2mat(GyroX(i));
    gyroYDataRaw=cell2mat(GyroY(i));
    gyroZDataRaw=cell2mat(GyroZ(i));
    gyroDataRaw=[gyroXDataRaw(1: limite ,1) gyroYDataRaw(1: limite ,1) gyroZDataRaw(1: limite ,1)];
    

%     %quaternion obatained from angvel
%     QuatW=getfield(WalkInertial(soggetto),'QuatX',stringa);
%     QuatX=getfield(WalkInertial(soggetto),'QuatX',stringa);
%     QuatY=getfield(WalkInertial(soggetto),'QuatX',stringa);
%     QuatZ=getfield(WalkInertial(soggetto),'QuatX',stringa);
%     QuatW=cell2mat(QuatW(i));
%     QuatX=cell2mat(QuatX(i));
%     QuatY=cell2mat(QuatY(i));
%     QuatZ=cell2mat(QuatZ(i));
%     %quaternion
%     quaternions=quaternion(QuatW  , QuatX , QuatY  ,  QuatZ );
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

function [newData] = butterworth(dataRaw , grade , cutoff,type)
    [a,b]=butter(grade,cutoff,type);
    %figure(10),freqz(a,b)
    newData=filtfilt(a,b,dataRaw);
end
end