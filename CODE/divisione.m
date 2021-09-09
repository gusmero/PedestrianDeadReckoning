% load('DataRaw\2021-01-02_15.23.28_ExpGusma_PercorsoCasa_70\ExpGusma_Session1_Number3_Calibrated_SD.mat')
% 


clear all
close all

load('DataRaw\2021-01-11_15.02.55_ExpGusma_SD_Session7\ExpGusma_Session7_Number3_Calibrated_SD.mat')
accDataRaw =[Number3_Accel_LN_X_CAL , Number3_Accel_LN_Y_CAL ,Number3_Accel_LN_Z_CAL];
gyroDataRaw = [ Number3_Gyro_X_CAL , Number3_Gyro_Y_CAL , Number3_Gyro_Z_CAL] ;
magDataRaw = [Number3_Mag_X_CAL  Number3_Mag_Y_CAL  Number3_Mag_Z_CAL];


%1:4612
%5131:9229
%9623:13780

%14540:18380
%18900:22340
%22840:26310

%26830:30170
%30600:33360
%33950:36840

int=33950:36840;

Number3_Accel_LN_X_CAL=accDataRaw(int  ,1);
Number3_Accel_LN_Y_CAL=accDataRaw(int ,2);
Number3_Accel_LN_Z_CAL=accDataRaw(int,3);

Number3_Gyro_X_CAL=gyroDataRaw(int ,1);
Number3_Gyro_Y_CAL=gyroDataRaw(int ,2);
Number3_Gyro_Z_CAL=gyroDataRaw(int,3);

Number3_Mag_X_CAL=magDataRaw(int,1);
Number3_Mag_Y_CAL=magDataRaw(int  ,2);
Number3_Mag_Z_CAL=magDataRaw(int ,3);



save('ExpGusma_Session9_Number3_Calibrated_SD.mat','Number3_Accel_LN_X_CAL','Number3_Accel_LN_Y_CAL','Number3_Accel_LN_Z_CAL','Number3_Gyro_X_CAL','Number3_Gyro_Y_CAL','Number3_Gyro_Z_CAL','Number3_Mag_X_CAL','Number3_Mag_Y_CAL','Number3_Mag_Z_CAL');



