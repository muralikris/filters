function [modelWeights,yEstimationAbsErr]= recursiveLeastSquaresWithForget(X,y,optionsArg)

% Recursive Least Squares with Forgetting Factor for adaptive system identification aka model learning

% M: Number of coefficients to learn in the model
% N: Number of Time periods in the data
% Ts: Data Sampling Interval (for information only)  

% System Model Representation:
% y = X'w + noise

% Input arguments:
% X = Model state matrix, dim MxN
% y = output vector, dim Nx1

% Output arguments:
% err = a priori estimation error, dim Nx1
% w = model weights/ coefficients, dim Mx1

if nargin < 3
    optionsArg = {};
end
options = ParseOptions (optionsArg);

epsilon=  options.epsilon ; % Convergence threshold for adaptive model learning
lambda= options.lambda ; %forgetting factor depends on sampling rate and avg on time
% memory= 1\ (1-lambda)
% high sampling rate means lot of redundant(fake) values so dont forget data fast lest you miss real useful data
% long on-time  means fewer cycles and we need to capture more cycles to identify system, hence longer memory 
if lambda>1
    lambda=1 ;
end

% Total no. of samples in data
N= length(y);

if options.useIntercept    
    interceptRow= ones(1,N); %Row vector
    X=[interceptRow; X];
end

M= size(X,1); % no. of model coefficients

% Initialization
delta= options.delta;% initial value, (1/delta) >0 a large number
P= eye(M)/delta;
modelWeights= zeros(M,1);

% Run Recursive Least Squares
convergeWeights= []; % Debug purpose
convergeTimeLengths= [];% Debug purpose

if N>2 % 2nd order system
    % error vector
    yEstimationAbsErr= [0 ; 0];
    
    % Loop
    convergeTimeLen= 0;
    for timInd= 3:N
        xVec= X(:,timInd);
        currErr= y(timInd)- xVec'*modelWeights; % Scalar
        
        tmp1=lambda + ( xVec' * P* xVec); % Scalar
        invTmp1= 1/tmp1; % Scalar
        k= invTmp1*P*xVec;
        
        tmpMat1= P* (xVec*xVec')*P;
        P=lambda^(-1)*(P-  invTmp1*tmpMat1);
        modelWeights= modelWeights+k*(currErr);
        absErr= abs(currErr);
        yEstimationAbsErr=[yEstimationAbsErr ; absErr]; %#ok<AGROW>
        
        if absErr < epsilon
            convergeTimeLen= convergeTimeLen +1;
        else
            if convergeTimeLen >60 %1hr
                convergeWeights=[convergeWeights, modelWeights ]; %#ok<AGROW>
                convergeTimeLengths=[convergeTimeLengths, convergeTimeLen]; %#ok<AGROW>
            end
            convergeTimeLen= 0;
        end
    end
        
else
    yEstimationAbsErr =[];
    modelWeights= [];
end

end

function options = ParseOptions (optionsArg)
	options.delta = GetNumericOption (optionsArg, 'delta', 0.005); % Initialization
    options.lambda  = GetNumericOption(optionsArg, 'lambda', 1); % <=1 
    options.epsilon  = GetNumericOption(optionsArg, 'epsilon', 0.001); % < 0.01 
    options.enablePlot = GetLogicalOption(optionsArg, 'enablePlot', false);
    options.useIntercept = GetLogicalOption(optionsArg, 'useIntercept', true); %learn intercept weight of model 
end
