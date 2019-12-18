%  Transition rate matrices
% format long;
R1 = [-33.65,15,15, 3.65,0, 0, 0, 0;...
     1095,-1113.65, 0, 0, 15, 3.65, 0, 0;...
     1095,0,-1113.65,0,15,0,3.65,0;...
     730,0,0,-760,0,15,15,0;...
     0,1095,1095,0,-2193.65,0,0,3.65;...
     0,730,0,1095,0,-1840,0,15;...
     0,0,730,1095,0,0,-1840,15;...
     0,0,0,0,730,1095,1095,-2920];

R2 = [-123.65,60,60, 3.65,0, 0, 0, 0;...
     1095,-1158.65, 0, 0, 60, 3.65, 0, 0;...
     1095,0,-1158.65,0,60,0,3.65,0;...
     730,0,0,-850,0,60,60,0;...
     0,1095,1095,0,-2193.65,0,0,3.65;...
     0,730,0,1095,0,-1885,0,60;...
     0,0,730,1095,0,0,-1885,60;...
     0,0,0,0,730,1095,1095,-2920];
 
 Rs=[-43.8,43.8;2190,-2190];
 Is=eye(2);
 R3=[-47.851,43.8,4.051,0;2190,-2199.582,0,9.582;...
     774.711,0,-818.511,43.8;0,1234.792,2190,-3424.792];
 I3=eye(4);
Iden = eye(8);
unk=[1;1;1;1;1;1;1;1];
p0=[1 0 0 0 0 0 0 0];
ps0=[1 0];
p3_0=[1 0 0 0];
 Qe=[-47.85,2190,774.711,0;43.8,-2199.582,0,1234.792;...
     4.05,0,-818.511,2190;1,1,1,1];
p=Qe\[0;0;0;1];
% step size in hours
delta= 1;

% sample time in hours
t_sample= 1;
% no .of steps in a sample
nsamples = t_sample/delta;
% total time in hours
t_final =8;
% total samples
nframes =t_final/t_sample;

%Matrix multiplication
R1e= Iden + ((R1*delta)/8760);
R2e= Iden + ((R2*delta)/8760);
Rse= Is + ((Rs*delta)/8760);
R3e= I3 + ((R3*delta)/8760);
% Steady State Probability and Transition Rates for loads in Normal Weather
R1_inf= R1e^nsamples;
p1_inf=p0*R1_inf;

% Steady State Probability and Transition Rates for loads in Rough Weather
R2_inf= R2e^nsamples;
p2_inf=p0*R2_inf;

% Steady State Probability and Transition Rates for weather
Rs_inf= Rse^nsamples;
ps_inf=ps0*Rs_inf;

% Memory Initialization
x= zeros(nframes+1,1);
x(1,:) = p3_0(1,3)+p3_0(1,4);

for n=1:nframes
    R3_t= R3e^(nsamples*n);
    p3_t= p3_0*R3_t;
    x(n+1,:)= p3_t(1,3)+p3_t(1,4);
end
R3_inf=R3e^nsamples;
t=0:t_sample:t_final;
plot(t,x(:,1),'-');
title('Probability Loss of load for Delta=1 hrs');
ylabel('P_L');
xlabel('time (hrs)');