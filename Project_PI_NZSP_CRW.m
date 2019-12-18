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
yi0=[0;0];
u0=[0;0.01];
ym=[1.82;40];

% Weighing Matrices
Q1 = [100,0,0,0;0,100,0,0;0,0,100,0;0,0,0,100];
R= [100 0;0 100];
Q2= [100 0;0 100];
S= [10 0;0 10];
Q_aug=[Q1,zeros(4,2),zeros(4,2);zeros(2,4),R,zeros(2);zeros(2,4),zeros(2),Q2];

% Discretization Block
[phi,gamma] = c2d (A,B,T);

% Augmented Matrices
A_aug=[phi,gamma,zeros(4,2);zeros(2,4),eye(2),zeros(2);T*C,T*D,eye(2)];
B_aug=[zeros(4,2);T*eye(2);zeros(2)];

% DLQR Cost Function
[K,S,E] = DLQR(A_aug,B_aug,Q_aug,S);

K1=K(1:2,1:4);% States Feed Back Gains
K2=K(1:2,5:6);% Control Rate Feed Back Gains
K3=K(1:2,7:8);% PI Control Gains

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
yi= zeros(nsamples+1,2);
ui= zeros(nsamples+1,2);

xs_n_hold= zeros(nframes+1,4);
us_n_hold= zeros(nframes+1,2);
yis_n_hold= zeros(nframes+1,2);
uis_n_hold= zeros(nframes+1,2);

% Exogenous Disturbance
q=10;
W= q*T*([1;1;1;1]);

[Pi12,Pi22]=QPMCALC((phi-Isys),gamma,C,D);
x_ss=Pi12*ym;
u_ss=Pi22*ym;

%Full state feedback
ui0= (K2*u_ss) + (K1*x_ss) - (K1*x0) - (K2*u0) -(K3*yi0);

% Initial Condition Settings for Loop Operations
xold= x0+W;
uold= u0;
yiold=yi0;
uiold=ui0;

x(1,:)  = (x0+W)';
u(1,:)  = (u0)';
yi(1,:) = (yi0)';
ui(1,:) = (ui0)';

% Initial Sample Holding
for k=1:thold
       xs_n_hold(k,:)=x(1,:);
       us_n_hold(k,:)=u(1,:);
       yis_n_hold(k,:)=yi(1,:);
       uis_n_hold(k,:)=ui(1,:);
 end

%sample and hold operation
for d=1:nsamples
    
    xnew = (phi*xold + gamma*uold) + W;
    unew = uold+ (T*uiold);
    
 % Control Saturation
    if(unew(1,1)>0.8)
        unew(1,1)=0.8;
    end
    if(unew(2,1)<0)
        unew(2,1)=0;
    end    
    
    yinew=(T*C*xold) + yiold -(T*ym) ;
    uinew= (K2*u_ss) + (K1*x_ss) - (K1*xnew) -(K2*unew) - (K3*yinew);
        
    % stacking and storing new varibles is memory
    x(d+1,:)  = xnew';
    u(d+1,:)  = unew';
    yi(d+1,:) = yinew';
    ui(d+1,:) = uinew';
    
    uold  = unew;
    xold  = xnew;
    yiold = yinew;
    uiold = uinew;
    
    for k=1:thold
       xs_n_hold((d*thold)+k,:)=x(d+1,:);
       us_n_hold((d*thold)+k,:)=u(d+1,:);
       yis_n_hold((d*thold)+k,:)=yi(d+1,:);
       uis_n_hold((d*thold)+k,:)=ui(d+1,:);
    end
end

xf=zeros(nframes+1,4);
tmp0=zeros(4,1);
x_area=zeros(4,1);
x_c=zeros(4,1);

t1=0:h:t_final;

% Sampled Data Block
[Ac,Bc] = c2d (A,B,h);
xf(1,:)=x0';

for n=1:nframes
tmp0= (Ac*(xf(n,:))') + (Bc*(us_n_hold(n,:))');
xf(n+1,:)=tmp0';
end 

% Closed Loop Tracking
A_cl=[phi,gamma,zeros(4,2);-T*K1,eye(2)-(K2*T),-K3*T;T*C,T*D,eye(2)];
B_cl=[zeros(4,2);T*(Pi22 + ( K1 * Pi12 ));zeros(2)];
C_cl=[C,zeros(2),zeros(2)];
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

% Plots Time Domain
figure(1);
 
subplot(2,2,1);
plot(t1,xf(:,1),'-',t1,xf(:,2),'-');
title('Buck Converter Time Response');
ylabel('Current (A), Voltage (V)');
xlabel('Time (Sec)');
%legend('IL','Va');
AXIS([0 t_final 0.0 6]);
GRID ON;

subplot(2,2,2);
plot(t1,xf(:,3),'-',t1,xf(:,4),'-');
title('Motor Time Response');
ylabel('Current (A), Speed (rad/sec)');
xlabel('Time (Sec)');
%legend('Ia','Wm');
AXIS([0 t_final 0.0 44]);
GRID ON;

subplot(2,2,3);
plot(t,us_n_hold(:,1),'-');
title('Buck Converter Control Time Response');
ylabel('Duty Cycle D');
xlabel('Time (Sec)');
AXIS([0 t_final 0 0.9]);
GRID ON;

subplot(2,2,4);
plot(t,us_n_hold(:,2),'-');
title('Motor Control Time Response');
xlabel('Time (Sec)');
ylabel('Load Torque N-m');
AXIS([0 t_final 0 0.06]);
GRID ON;
