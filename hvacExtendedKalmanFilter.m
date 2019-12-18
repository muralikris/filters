function hvacModelParameters= hvacExtendedKalmanFilter(inputTrainingData,outputTrainingData, optionsArg)

% Discrete-time Extended Kalman Filter for realtime adaptive system identification for nonlinear systems

% For a discrete system at time-sample k, with states(x_k, contains control inputs) and observed outputs (y_k)
% Model Representation is as follows:
% x_k+1= f_k(x_k) + e_k
% y_k= g_k(x_k) + n_k
% F_k= df_k/dx_k.(x_k), Jacobian
% G_k= dg_k/dx_k.(x_k), Jacobian
% Q_k =E(e_k.e_k') process noise matrix
% R_k=E(n_k.n_k')) observation noise matrix
% P_k= covariance matrix 

if nargin < 3
    optionsArg = {};
end

options = ParseOptions (optionsArg);
epsilon=  options.epsilon ; % Convergence threshold for adaptive model learning
delta= options.delta;% initial value for weighing matrices, 
ts=options.resampleInMin *60; %  Sampling Time In Seconds 

hvacModelParameters= [];
if isempty(inputTrainingData) || isempty(outputTrainingData)
    return
end

% Smart Thermostat Hvac Data Analytics

% ARMAX Model : y(k)= a1*y(k-1) + a2*y(k-2)+ a3*Tout(k-1) + a4*u(k-1)
%yk: Inside Air Temperature at sample k
N= 2;% N: Number of internal states of the model, 2nd order hvac model 

% Inputs/ Controls
%u: Hvac On Off State
%Tout: Outside Weather Temperature

% Outputs
%yk: Inside Air Temperature 

% Training Data Initialization
T_out= inputTrainingData(1,:); %Outside Weather Temperature 
hvac_state= inputTrainingData(2,:);% Hvac On Off State
% u: can be fractional for multi-stage
T_in= outputTrainingData(1,:); %Inside Air Temperature

% Parameters used in Model: a1-a4
W= 4;% Total No.of parameters
M= N+W; % Total number of states for EKF

L= length(T_in);% Total no. of samples in data

% Run Recursion for EKF
if L>10 % sufficient number of samples
    
    % Initializing Weighing Matrices
    yEst= T_in(1);
    zEst= T_in(1);
    %estimatedParameters= [65;68;0.2; 0.01;0.02;0.03];
    estimatedParameters= zeros(M,1);   

    P= [100 0 0 0 0 0;
        0 100 0 0 0 0;
        0 0 100 0 0 0;
        0 0 0 100 0 0;
        0 0 0 0 1 0;
        0 0 0 0 0 100];%MxM P0= E{(X0-Xe0)(X0-Xe0)'} ;P0!=0
    
    % Process Noise
    Q= [10 0 0 0 0 0;
        0 10 0 0 0 0;
        0 0 10 0 0 0;
        0 0 0 10 0 0;
        0 0 0 0 0.1 0;
        0 0 0 0 0 10];%MxM
    
    % Measurement Noise-low for Dreamstat2
    R= 10;%1x1        

    % 2nd order HVAC Model                 
    for timInd= 2:L
        prevTimeInd= timInd-1;
        Tout= T_out(prevTimeInd);
        u= hvac_state(prevTimeInd);
        Tin= estimatedParameters(1);
        Tin_k_1= estimatedParameters(2);
        a1= estimatedParameters(3);
        a2= estimatedParameters(4); 
        a3= estimatedParameters(5); 
        a4= estimatedParameters(6); 
        
        % State space representation for the above ARMAX model
        %[yk; yk-1]= [a1, a2; 1, 0]* [yk-1; yk-2]+ [a3 a4; 0; 0]*[uk-1 ; Tout_k-1];
        % yk-1= [1,0]* [yk-1; yk-2];                
        f_x = [(a1*Tin) + (a2*Tin_k_1) + (a3*u) + (a4*Tout);
            Tin;
            a1;
            a2;
            a3;
            a4];
        
        g_x= Tin;
        
        %states = [Tin; Tin_k_1; a1; a2; a3; a4]
        % F- Jacobian of f_x wrt states: 6x6
        % G- Jacobian of g_x wrt states: 6x6
        F= [ a1, a2, Tin, Tin_k_1, u, Tout ;
            1, 0, 0, 0, 0, 0 ;
            0, 0, 1, 0, 0, 0 ;
            0, 0, 0, 1, 0, 0 ;
            0, 0, 0, 0, 1, 0;
            0, 0, 0, 0, 0, 1];
        
        G= [1, 0, 0, 0, 0, 0] ;       
        
        % EKF Implementation
        x_hat_k_m = f_x;
        P_k_m = (F*P*F') + Q;        
        K_k =(P_k_m*G')/( (G*P_k_m*G') + R);
        yCurr= T_in(prevTimeInd);        
        err_k= yCurr - g_x;
        x_hat_k = x_hat_k_m + (K_k*err_k) ;                                    

%         if any(isnan(x_hat_k))
%             debug=1;
%         end
        P = (eye(M) - K_k*G)*P_k_m;
        estimatedParameters= x_hat_k;
        yEst=[yEst, Tin]; %#ok<AGROW>
        zEst=[zEst, Tin_k_1]; %#ok<AGROW>
        
    end
    
    if options.enablePlot
        figure
        plot(T_in);
        hold on
        plot(yEst);
        hold on
        plot (zEst);
        hold off
        legend('airTempMeas','airTempEst','airDiffTempEst');
    end
end

hvacModelParameters= estimatedParameters(3:end);
end


function options = ParseOptions (optionsArg)
	options.delta = GetNumericOption (optionsArg, 'delta', 0.005); % Initialization
    options.lambda  = GetNumericOption(optionsArg, 'lambda', 1); % <=1 
    options.epsilon  = GetNumericOption(optionsArg, 'epsilon', 0.001); % < 0.01 
    options.onSwing  = GetNumericOption(optionsArg, 'onSwing',0.5);% In deg fahrenheit
    options.offSwing  = GetNumericOption(optionsArg, 'offSwing',0);% In deg fahrenheit
    options.enablePlot = GetLogicalOption(optionsArg, 'enablePlot', false);   
    options.resampleInMin  = GetNumericOption(optionsArg, 'resampleInMin', 1);
end