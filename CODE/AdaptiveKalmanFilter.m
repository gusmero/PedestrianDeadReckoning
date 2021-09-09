

function[X]= AdaptiveKalmanFilter(accDataRaw,magDataRaw,gyroDataRaw,accDataRawInit,gyroDataRawInit,magDataRawInit)
%% TIME LINE

%frequenza & linespace
N = size(accDataRaw,1);
Fs = 128;
Ts=1/Fs;
t = 0:Ts:((N-1)*(1/Fs));  

%% INIZIALIZZAZIONE

c=1.3;
c0=15;
 

[q_am_Init , q_w_Init]=quaternion(accDataRawInit,magDataRawInit,gyroDataRawInit);

[q_am,q_w]=quaternion(accDataRaw,magDataRaw,gyroDataRaw);
% inizializzazione variabili del filtro di kalman  

X(1,:)=[1 0 0 0];

%
H=eye(4);

%NOISE MODEL COVARIANCE
P=10^(-6)*eye(4);

%NOISE MESUREMENT MATRIX
R=diag([var(q_am_Init(:,1)), var(q_am_Init(:,2)) , var(q_am_Init(:,3)) , var(q_am_Init(:,4)) ]);
% R=10^(-6)*eye(4);

%NOISE STATE MATRIX
Q=diag([var(q_w_Init(:,1)), var(q_w_Init(:,2)) , var(q_w_Init(:,3)) , var(q_w_Init(:,4)) ]);
% Q=10^(-10)*eye(4);

%MESUREMENT VECTOR
Z=zeros(N,4);


%% Adaptive Kalman  Filter

for k=2:N


%TRASFORMATION MATRIX
F=(eye(4)+(Ts/2.*quaternionMatrix([0 gyroDataRaw(k,:)])));

% STATE PREDICTION
X_hat=F*X(k-1,:)';

%NOISE MODEL COVARIANCE PREDICTION
P_hat=F*P*F'+Q;

%to avoid discontinues MESUREMENT VECTOR
if q_am(k,:)*X_hat > 0
    Z(k,:)=q_am(k,:)/norm(q_am(k,:));
else
    Z(k,:)=-q_am(k,:)/norm(-q_am(k,:));
end
    
% NORMALIZATION OF MESUREMENT VECTOR
Z(k,:)=Z(k,:)/norm(Z(k,:)); 


%Z_hat PREDICTED MESUREMENT VECTOR
Z_hat=eye(4)*X_hat;

% RESIDUAL CALCULATION
for i=1:size(Z,2)
    residual(i)=Z(i)-Z_hat(i);
end

% STANDARD RESIDUAL
r_primo=abs(residual./mad(residual));

% P_WEIGH

for i=1:4
    if r_primo <= c
        P_wgh(i,i)=1/R(i,i);
    else
        P_wgh(i,i)=(c/r_primo(1,i))*(1/R(i,i));
    end
end
%Calculation of P_wgh non-digonal
for i=1:4
    for j=4:1
        if r_primo(1,i) <= c && r_primo(1,j)<=c
            P_wgh(i,i)=1/R(i,j);
        elseif r_primo(1,i) > c || r_primo(1,j) > c
                P_wgh(i,i)=(c/max(r_primo(1,i),r_primo(1,j)))*(1/R(i,j));
        end
    end
end

% LIST-SQUARE ESTIMATOR
X_tilde=inv(F'*P*F)*F'*P*Z';
%X_tilde=(F'*P_wgh*F)/(F'*P_wgh*Z);

%STATE DISCREPANCY STATISTICS
delta_X_tilde=norm(X_tilde-X_hat)/sqrt(trace(P_hat));

% ADAPTIVE FACTOR

if delta_X_tilde <= c0
    iota=1;
else
    iota=c0/delta_X_tilde;
end   
    
% GAIN MATRIX
K=1/iota.*P_hat*H'*inv(1/iota.*H*P_hat*H'+inv(P));

% UPDATING TEH STATE VECTORE
X(k,:)=(X_hat+K*residual')';

%NORMALIZATION OF STATE VECTOR
X(k,:)=X(k,:)/norm(X(k,:));

%updatingstate varriance matrix
P=(eye(4)-K*H)*P_hat;

end
  


 %% PRINT 
 
 eulafter=attitudeHeading(X);
 
 eulbefore=attitudeHeading(q_w);
 
 printSignal3D(gyroDataRaw(:,1)  ,  gyroDataRaw(:,2)  , gyroDataRaw(:,3)   ,1:N,'Gyro Data Raw','sample','Gyro(deg/s^2)',2);

printSignal3D(eulbefore(1,:),eulbefore(2,:),eulbefore(3,:),1:N,'Orientamento from giroscope','sample','degre',3);

printSignal3D(eulafter(1,:),eulafter(2,:),eulafter(3,:),1:N,'Orientamento filtrato','sample','degre',4);

 printSignal(eulafter(3,:),1:N ,'Orientamento from giroscope','sample','degre',6)
end

%% FUNCTION

function [eul]=attitudeHeading(X)
    eul=zeros(3,size(X,1));
    for i=1:size(X,1)
        eul(:,i)=quat2eul(X(i,:),"XYZ");
    end
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

function print2Signal(signalX,signalY,t,Title,Xax,Yax,indexFigure)
    figure(indexFigure)
    plot(t, signalX)
    hold on
    plot(t,signalY)
    xlabel(Xax)
    ylabel(Yax)
    title(Title)
    legend('Xax','Yax','Zax')
end

function printSignal(signal,t,Title,Xax,Yax,indexFigure)
    figure(indexFigure)
    plot(t, signal)
    xlabel(Xax)
    ylabel(Yax)
    title(Title)
end




function [q_am,q_w]=quaternion(accDataRaw,magDataRaw,gyroDataRaw)
    
   N = size(accDataRaw,1);
Fs = 128;
Ts=1/Fs;

    % HEADING CALCULATION
    a_x=accDataRaw(:,1);
    a_y=accDataRaw(:,2);
    a_z=accDataRaw(:,3);


    for i=1:N
        if accDataRaw(i,3) >= 0
            q_a(i,:)=[ sqrt((a_z(i) + 1)/2 )  -a_y(i)/(2*sqrt((a_z(i) + 1) / 2))     a_x(i)/(2*sqrt((a_z(i) + 1) / 2))  0  ];
        else
            q_a(i,:)=[ -a_y(i)/(2*sqrt((1-a_z(i)) / 2))  sqrt((1-a_z(i)) /2)    0    a_x(i)/(2*sqrt((1-a_z(i)) / 2))];
        end
    end

    %from magnetometer

    % calcolo il vettore del campo magnetico ruotato
    l=zeros(N,3);
    for i=1:N
        l(i,:)=quat2rotm(q_a(i,:))*magDataRaw(i,:)';
    end
    l_x=l(:,1);
    l_y=l(:,2);
    l_z=l(:,3);

    for i=1:N
        if l(i,3) >= 0
            q_m(i,:)= [(sqrt(l_x(i)^2 + l_y(i)^2 +  l_x(i) * sqrt(l_x(i)^2 + l_y(i)^2)))/sqrt(2*(l_x(i)^2 + l_y(i)^2))   0   0  l_y(i)/(sqrt(2)*sqrt(l_x(i)^2 + l_y(i)^2 +  l_x(i) * sqrt(l_x(i)^2 + l_y(i)^2))) ];
        else
            q_m(i,:)= [l_y(i)/(sqrt(2)*sqrt(l_x(i)^2 + l_y(i)^2 -  l_x(i) * sqrt(l_x(i)^2 + l_y(i)^2)))  0    0    (sqrt(l_x(i)^2 + l_y(i)^2 -  l_x(i) * sqrt(l_x(i)^2 + l_y(i)^2)))/sqrt(2*(l_x(i)^2 + l_y(i)^2))];
        end
    end
   

    % FUSION
    for i=1:N
        
        q_am(i,:)=quatmultiply(q_a(i,:),q_m(i,:));
    end

    % HEADING ESTIMATION USING ANGULAR RATE
    q_w(1,:)=[1 0 0 0];
    for i=2:N
        q_w(i,:)=(eye(4)+(Ts/2.*quaternionMatrix([0 gyroDataRaw(i,:)])))*q_w(i-1,:)';
    end

end

function [quatMet]=quaternionMatrix(q)
    quatMet =[ q(1)  -q(2)  -q(3)  -q(4)
               q(2)  q(1)   -q(4)   q(3)
               q(3)   q(4)   q(1)    -q(2)
              q(4)   -q(3)   q(2)    q(1)];
end


function [Rrpy]= rotationMatrix(phi,theta,psi)
    cphi=cos(deg2rad(phi));
    ctheta=cos(deg2rad(theta));
    cpsi=cos(deg2rad(psi));
    sphi=sin(deg2rad(phi));
    stheta=sin(deg2rad(theta));
    spsi=sin(deg2rad(psi));
  
    Ry= [ cphi 0 sphi; 0 1 0 ; -sphi 0  cphi];
    Rz = [ ctheta -stheta 0 ; stheta ctheta 0 ; 0 0 1 ];
    Rx = [  1 0 0 ; 0  cpsi -spsi ; 0 spsi cpsi];
    
    Rrpy=Ry*Rz*Rx;
end


function [Rrpy]= rotationMatrixQuaternion(q)
    
    Rrpy=[q(1)^2+q(2)^2-q(3)^2-q(4)^2   2*(q(2)*q(3)+q(1)*q(4))   2*(q(2)*q(4)-q(1)*q(3))
          2*(q(2)*q(3)-q(1)*q(4))    q(1)^2-q(2)^2+q(3)^2-q(4)^2  2*(q(3)*q(4)+q(1)*q(2))
          2*(q(2)*q(4)+q(1)*q(3))     2*(q(3)*q(4)-q(1)*q(2))      q(1)^2-q(2)^2-q(3)^2+q(4)^2];  
    
end


% function [phi , theta,psi]= InverseSolutionQuaternio(q)
%     phi=arctan
% end

