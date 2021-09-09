%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%                                                               %%%%%%
%%%%%%                      Kalman Filter                            %%%%%%
%%%%%%                                                               %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% NOTE
% IN QUESTA VERSIONE UTILIZZEREMO LA TECNICA DI INTEGRAZIONE DI EULERO DI 



function[X]= LinearKalmanFilter(accDataRaw,gyroDataRaw,accDataRawInit,gyroDataRawInit)

%% CARICAMENTO BASELINE

gyroDataRawInit=[];
accDataRawInit=[];

load('DataRaw\BASELINE\2020-12-07_22.49.42_ExpGusma_MultiSession\ExpGusma_Session1_Number3_Calibrated_SD.mat');
accDataRawInit =[accDataRawInit ;Number3_Accel_LN_X_CAL , Number3_Accel_LN_Y_CAL ,Number3_Accel_LN_Z_CAL];
gyroDataRawInit = [gyroDataRawInit ; Number3_Gyro_X_CAL , Number3_Gyro_Y_CAL , Number3_Gyro_Z_CAL] ;

load('DataRaw\BASELINE\2020-12-07_22.49.42_ExpGusma_MultiSession\ExpGusma_Session2_Number3_Calibrated_SD.mat');
accDataRawInit =[accDataRawInit ;Number3_Accel_LN_X_CAL , Number3_Accel_LN_Y_CAL ,Number3_Accel_LN_Z_CAL];
gyroDataRawInit = [gyroDataRawInit ; Number3_Gyro_X_CAL , Number3_Gyro_Y_CAL , Number3_Gyro_Z_CAL] ;

load('DataRaw\BASELINE\2020-12-07_22.49.42_ExpGusma_MultiSession\ExpGusma_Session3_Number3_Calibrated_SD.mat');
accDataRawInit =[accDataRawInit ;Number3_Accel_LN_X_CAL , Number3_Accel_LN_Y_CAL ,Number3_Accel_LN_Z_CAL];
gyroDataRawInit = [gyroDataRawInit ; Number3_Gyro_X_CAL , Number3_Gyro_Y_CAL , Number3_Gyro_Z_CAL] ;

load('DataRaw\BASELINE\2020-12-07_22.49.42_ExpGusma_MultiSession\ExpGusma_Session4_Number3_Calibrated_SD.mat');
accDataRawInit =[accDataRawInit ;Number3_Accel_LN_X_CAL , Number3_Accel_LN_Y_CAL ,Number3_Accel_LN_Z_CAL];
gyroDataRawInit = [gyroDataRawInit ; Number3_Gyro_X_CAL , Number3_Gyro_Y_CAL , Number3_Gyro_Z_CAL] ;




%% TIME LINE

%frequenza & linespace
N = size(accDataRaw,1);
Fs = 51.2;
Ts=1/Fs;
t = 0:Ts:((N-1)*(1/Fs));  


%% INITIALIZATION 

%NOISE MESUREMENT MATRIX
% Baseline orientation 
vAng=gyroDataRawInit;
angDiInit=zeros(size(gyroDataRawInit,1),3);
for i=2:size(gyroDataRawInit,1)
    angDiInit(i,:)=angDiInit(i-1,:)+vAng(i,:).*Ts;
end
vAng=gyroDataRaw;
angDi=zeros(N,3);
for i=2:N
    angDi(i,:)=angDi(i-1,:)+vAng(i,:).*Ts;
end

%variance
varianceMesurex=var(gyroDataRawInit(:,1));
varianceMesurey=var(gyroDataRawInit(:,2));
varianceMesurez=var(gyroDataRawInit(:,3));
R=diag([varianceMesurex , varianceMesurey , varianceMesurez]);
% R=diag([1,1,1]).*10^(-8);
Z=[varianceMesurex,varianceMesurey,varianceMesurez];

%STATE VECTOR
meanBiasGyrox=mean(gyroDataRawInit(:,1));
meanBiasGyroy=mean(gyroDataRawInit(:,2));
meanBiasGyroz=mean(gyroDataRawInit(:,3));
X=zeros(6,N);
X(:,1)=[0; 0; 0; 0; 0; 0];

% %IMPUT SYSTEM
% U=gyroDataRaw';

%ERROR COVERIANCE MATRIX
varianceBiasGyrox=var(gyroDataRawInit(:,1));
varianceBiasGyroy=var(gyroDataRawInit(:,2));
varianceBiasGyroz=var(gyroDataRawInit(:,3));
varianceAnglex=var(gyroDataRawInit(:,1));
varianceAngley=var(gyroDataRawInit(:,2));
varianceAnglez=var(gyroDataRawInit(:,3));
Q=diag([varianceAnglex,varianceAngley,varianceAnglez,varianceBiasGyrox,varianceBiasGyroy,varianceBiasGyroz]);

%NOISE MODEL COVARIANCE
% P=diag([10^(-8),10^(-8),10^(-8),10^(-8),10^(-8),10^(-8)]);
P=diag([1,1,1,1,1,1]).*10^(-8);

% STATE TRANSIOZIONE MATRIX
A=diag([1,1,1,1,1,1])+ [zeros(3 ,3) diag([-Ts,-Ts,-Ts]); zeros(3 ,6)];

%MESUREMENT TRANSITION MATRIX
H=[diag([1,1,1]) zeros(3,3)];

%MESUREMENT EXIT VECTOR
Y=zeros(3,N);

%MESUREMENT MATRIX
B=[diag([Ts,Ts,Ts]) zeros(3,3) ; zeros(3,6)];

for k=1:N-1
    %%  PREDICTION pHASE
    U=[gyroDataRaw(k,:)'; 0;0;0];    
    
    %PREDICTION OF THE STATE VECTOR
    Xpredict=A*X(:,k)+B*U;
    
    %PREDICTION OF ERROR COVARIANCE MATRIX
    Ppredict=A*P*A'+Q;

    %MESUREMENT EXIT VECTORs
    Y(:,k)=H*X(:,k)+ Z';
    
    %% UPDATE PHASE
 
    % KALMAN GAIN
    K= Ppredict*H'*(H*Ppredict*H'+R)^(-1);
    
    %non sono sicuro della corretteza di questa parte
    % INNOVATION MATRIX
    InnMat=Y(:,i)-H*Xpredict;
    
    %STATE VECTOR UPDATE
    X(:,k+1)=Xpredict + K *InnMat;
    
    %ERROR COVARIANCE MATRIX UPDATE
    P= (eye(6) - K*H)*Ppredict;
    
end

%% SIGNAL ANALYSIS

% filtered with kalman filter
printSignal3D(X(1,:),X(2,:),X(3,:),1:N,'Linear Kalman filter performance for orientation estimatinon','sample','degre',4);

%unfiltered orientation signal
printSignal3D(angDi(:,1),angDi(:,2),angDi(:,3),1:N,'Unfilter performance for orientation','sample','degre',5);

%comparison
print2Signal(X(1,:),angDi(:,1),1:N,'comparison betwen filtered and unfiltered','sample','degre',6)


%% FUNCTION

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

function print2Signal(signalX,signalY,t,Title,Xax,Yax,indexFigure)
    figure(indexFigure)
    plot(t, signalX)
    hold on
    plot(t,signalY)
    xlabel(Xax)
    ylabel(Yax)
    title(Title)
    legend('kalman filter','unfiltered')
end

function printSignal(signal,t,Title,Xax,Yax,indexFigure)
    figure(indexFigure)
    plot(t, signal)
    xlabel(Xax)
    ylabel(Yax)
    title(Title)
end

function [gyroDataRaw]= acquisizione_dati_gyro(WalkInertial,soggetto,stringa,i)
   
    %3-axis gyroscope
    GyroX=getfield(WalkInertial(soggetto),'GyroX',stringa);
    GyroY=getfield(WalkInertial(soggetto),'GyroY',stringa);
    GyroZ=getfield(WalkInertial(soggetto),'GyroZ',stringa);
    %The unit of the calibrated accelerometer is meters per square second(m/s2)
    gyroXDataRaw=cell2mat(GyroX(i));
    gyroYDataRaw=cell2mat(GyroY(i));
    gyroZDataRaw=cell2mat(GyroZ(i));
    gyroDataRaw=[gyroXDataRaw(1:30000,1) gyroYDataRaw(1:30000,1) gyroZDataRaw(1:30000,1)];
   
end
function [accDataRaw]=acquisizione_dati_acc(WalkInertial,soggetto,stringa,i)
    AccX=getfield(WalkInertial(soggetto),'AccX',stringa);
    AccY=getfield(WalkInertial(soggetto),'AccY',stringa);
    AccZ=getfield(WalkInertial(soggetto),'AccZ',stringa);
    %The unit of the calibrated accelerometer is meters per square second(m/s2)
    accXDataRaw=cell2mat(AccX(i));
    accYDataRaw=cell2mat(AccY(i));
    accZDataRaw=cell2mat(AccZ(i));
    accDataRaw=[accXDataRaw(1:30000,1) accYDataRaw(1:30000,1) accZDataRaw(1:30000,1)];
end

function [accDataRaw , gyroDataRaw ]= acquisizione_dati(WalkInertial,soggetto,stringa,i)
    
    %3-axes acceleration
    AccX=getfield(WalkInertial(soggetto),'AccX',stringa);
    AccY=getfield(WalkInertial(soggetto),'AccY',stringa);
    AccZ=getfield(WalkInertial(soggetto),'AccZ',stringa);
    %The unit of the calibrated accelerometer is meters per square second(m/s2)
    accXDataRaw=cell2mat(AccX(i));
    accYDataRaw=cell2mat(AccY(i));
    accZDataRaw=cell2mat(AccZ(i));
    accDataRaw=[accXDataRaw(1:23821,1) accYDataRaw accZDataRaw];
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
    

   
end

function [newData] = butterworth(dataRaw , grade , cutoff,type)
    [a,b]=butter(grade,cutoff,type);
    %figure(10),freqz(a,b)
    newData=filtfilt(a,b,dataRaw);
end


end