

s1='DataRaw\';
s0={'70ppm\','85ppm\','100ppm\'};
s3='ExpGusma_Session';

s5='_Number3_Calibrated_SD.mat';

for i=1:9

if i==1 || i==2 ||i==3
    s2=s0(1);
elseif i==4 || i==5 || i==6
    s2=s0(2);
else
    s2=s0(3);
end
s4=string(i);
stringa=strcat(s1,s2(1),s3,s4,s5);

load(stringa)
accDataRaw =[Number3_Accel_LN_X_CAL , Number3_Accel_LN_Y_CAL ,Number3_Accel_LN_Z_CAL];
gyroDataRaw = [ Number3_Gyro_X_CAL , Number3_Gyro_Y_CAL , Number3_Gyro_Z_CAL] ;

s(i)=stepDetection(accDataRaw,gyroDataRaw);

end


% %% verifica della step length 
% for i=1:9
%     strideVelocity(i)=getfield(s(i),'strideVelocity')
%     strepFrequency(i)=getfield(s(i),'stepFrequency')
% end
% 
% step=strideVelocity.* stepFrequency;
