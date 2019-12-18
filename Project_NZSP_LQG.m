% * * * * *  * * * * *  * * * * *  * * * * *  * * * * *  * * * * *  * %
% * *  * * * * *  DC MOTOR Control Using Buck Converter  * * * * *  * %
% * * * * *  * * * * *  AERO-625 Fall 2010  * * * * *  * * * * *  * %
% * * * * *  * * * * *  * * * * *  * * * * *  * * * * *  * * * * *  * %

% Time History Settings
h=5e-5;% Integration Time Step
T=5e-4; % Sampling Time For System
t_final=0.2; % Final Time

% Linearized System Model For DC-DC Step Down Converter Fed DC Motor
A=[0,-50,0,0;2500,0,-2500,0;0,380.23,-760.46,-17.49;0,0,651.56,-11.93];
B=[600,0;0,0;0,0;0,-14164.3];
C=[0 0 1 0;0 0 0 1];
D=zeros(2);
Isys=eye(4);
% %States
%x1: Inductor Current
%x2: Capacitor Voltage
%x3: Motor Current
%x4: Motor Speed

% %Controls
%D: Duty Cycle
%TL: Load Torque


% %Outputs
%x3: Motor Current
%x4: Motor Speed

% %Initial Conditions
x0=[0;0;0;0];
ym=[1.82;40];
xe0=[0;0;0;0];%Xe0=E[X0]
P0= [10 0 0 0;0 10 0 0;0 0 1 0;0 0 0 1];%P0= E{(X0-Xe0)(X0-Xe0)'} ;P0!=0

% Weighing Matrices
Qnz = [100,0,0,0;0,100,0,0;0,0,100,0;0,0,0,10];
Rnz= [100 0;0 100];
[K,Qc,Rc,M,P,E]=LQRDJV(A,B,Qnz,Rnz,T);
%time vector
t= 0:h:(t_final+T-h);

% total no of sample points
nframes= t_final/h;
nsamples= t_final/T;

%hold interval
thold= T/h;

% Memory Initialization
x= zeros(nsamples+1,4);
u= zeros(nsamples+1,2);
xe= zeros(nsamples+1,4);
xes_n_hold= zeros(nframes+1,4);
us_n_hold= zeros(nframes+1,2);

% Discretization Block
[phi,gamma] = c2d (A,B,T);
[Pi12,Pi22]=QPMCALC((phi-Isys),gamma,C,D);
x_ss=Pi12*ym;
u_ss=Pi22*ym;

% Process Noise Co-variance Matrices due to Modeling Uncertainities

  %sys noise for inductor current (ESR) =0.1 A2, capacitor voltage (ESR)= 0.1 V2
  %sys noise for motor current=0.05 A2, motor speed = 0.044rad2/s2
    Q=[0.1 0 0 0;0 0.1 0 0;0 0 0.05 0;0 0 0 0.044];
    G=0.1*eye(4);

  % Measurement noise covariance matrix
    R = [0.024 0;0 0.02];
    V=0.01*eye(2);

% Noise vector Generation
  randn('state',0);
  % Process Noise
    w= [sqrt(0.1)*randn(1,nsamples+1);sqrt(0.1)*randn(1,nsamples+1);sqrt(0.05)*randn(1,nsamples+1);sqrt(0.044)*randn(1,nsamples+1)] ;
  % Measurement Noise
    v= [sqrt(0.024)*randn(1,nsamples+1);sqrt(0.02)*randn(1,nsamples+1)] ;

% Discrete time State Estimation
  [Md,Pd,Zd,Ed] = DLQE(phi,G,C,Q,R);

% Initial Condition Settings
x(1,:) = (x0)';
xe(1,:)= (xe0)';

% Holding 1st Sample
for k=1:thold
       xes_n_hold(k,:)=xe(1,:);
       us_n_hold(k,:)=u(1,:);
end
 
xold= x0+(G*w(:,1)) ;

%Control Law
uold= u_ss+ (K*(x_ss-xe0));

% Control Saturation
if(uold(1,1)>0.8)
    uold(1,1)=0.8;
end
if(uold(2,1)<0)
    uold(2,1)=0;
end   

yold= C*x0 + (V*v(:,1));

u(1,:) = (uold)';

% Discrete time Kalman filter:
   
 % Prediction for state vector and covariance:
   xe_old = phi*xe0 + gamma*uold; % Time Update
   P_old = (phi * P0 * phi') + (G*Q*G'); % Time Update
   
%a posteriori recursion
 for d=1:nsamples
 %Plant
   xnew= phi*xold + gamma*uold + G*w(:,d+1);
   ynew= C*xnew + (V*v(:,d+1));
   x(d+1,:) = (xnew)';   
   
   % Compute Kalman gain factor:
   L = P_old * C'* inv(( C * P_old * C')+ R);

   % Correction based on observation:
   xe_new = xe_old + L*(ynew-C*xe_old); %Measurement Update
   P_new = P_old - (L * C * P_old); %Measurement Update
   xe(d+1,:) = (xe_new)';

   % NZSP Controller
   u_new= u_ss+ (K*(x_ss-xe_new));
   
   % Control Saturation
    if(u_new(1,1)>0.8)
        u_new(1,1)=0.8;
    end
    if(u_new(2,1)<0)
        u_new(2,1)=0;
    end       
    
   u(d+1,:) = (u_new)';
   
   % Estimator
   xe_old = phi*xe_new + gamma*u_new; % Time Update
   P_old = (phi * P_new * phi') + (G*Q*G'); % Time Update  
   
   for k=1:thold
       xes_n_hold((d*thold)+k,:)=xe(d+1,:);
       us_n_hold((d*thold)+k,:)=u(d+1,:);
   end
   
   uold = u_new;
   xold = xnew;
 end

xf=zeros(nframes+1,4);
tmp0=zeros(4,1);
t1=0:h:t_final;

% Sampled Data Block- Discreter Inputs to Continuous Plant
[Ac,Bc] = c2d (A,B,h);
xf(1,:)=x0';

for n=1:nframes
tmp0= (Ac*(xf(n,:))') + (Bc*(us_n_hold(n,:))');
xf(n+1,:)=tmp0';
end

Lk=[0.1022    0.0703;
    0.3384    0.1659;
    0.1525    0.0941;
    0.0784    0.2483];

% Closed Loop Tracking Robustness

A_cl=[phi, -gamma*K; (phi -(Lk*C)-(gamma*K)),(Lk*C)];
B_cl=[gamma * (Pi22 + ( K * Pi12 ));gamma * (Pi22 + ( K * Pi12 ))];
C_cl=[C C];
sys2 = ss(A_cl,B_cl,C_cl,0,T);

% High Frequency m(w) Specifications
num_mw=[1 2];
den_mw=20;
mw_i = tf(den_mw,num_mw);

% Low Frequency Specificationa
low_fr = tf(10,[1 0]);

% Closed Loop Frequency Bounds
 figure(2);
 hold all
 sigma(low_fr);
 sigma(sys2);
 sigma(mw_i);
 
% Time Response Plots

figure(1);
subplot(2,2,1);
plot(t1,xf(:,1),'-',t,xes_n_hold(:,1),'-',t1,xf(:,2),'-',t,xes_n_hold(:,2),'-');
title('Buck Converter Time Response');
ylabel('Current (A), Voltage (V)');
xlabel('Time (Sec)');
legend('IL','ILe','Va','Vae');
AXIS([0 t_final 0.0 8]);
GRID ON;

subplot(2,2,2);
plot(t1,xf(:,3),'-',t,xes_n_hold(:,3),'-',t1,xf(:,4),'-',t,xes_n_hold(:,4),'-');
title('Motor Time Response');
ylabel('Current (A), Speed (rad/sec)');
xlabel('Time (Sec)');
legend('Ia','Iae','Wm','Wme');
AXIS([0 t_final 0.0 50]);
GRID ON;

subplot(2,2,3);
plot(t,us_n_hold(:,1),'-');
title('Buck Converter Control Time Response');
ylabel('Duty Cycle D');
xlabel('Time (Sec)');
AXIS([0 t_final 0 1]);
GRID ON;

subplot(2,2,4);
plot(t,us_n_hold(:,2),'-');
title('Motor Control Time Response');
xlabel('Time (Sec)');
ylabel('Load Torque N-m');
AXIS([0 t_final 0 0.2]);
GRID ON;
