clear;

% Add search path
addpath(genpath(pwd));

UK = readtable("UK_yields.csv");
date = table2array(UK(:, 1));
date = datetime(date, 'InputFormat', 'dd/MM/yyyy');
UK = table2array(UK(:, 2: end));

US_factors = readtable("US_factors.csv");
US_factors = table2array(US_factors);

maturity = [12, 24, 36, 60, 84, 120, 240]';

%% DNS model
yield = UK; % in-sample yields 
n_period = 12; % number of forecasting period
n_contract = size(yield, 2); % number of contracts
n_factor = 3; % number of NS factors
n_obs = size(yield, 1); % number of observations 
mdl1 = ssm(@(params) DNS(params, yield, maturity));

A0 = [3; 3; 3];   
B0 = [0.1; 0.1; 0.1];   
mu0 = [1; 1; 1];
%D0 = repelem(0.1, n_contract)'; % Diagonal covariance matrix for measurement equation
D0 = [repelem(0.1, n_contract)'; 0.1]; % Band covariance matrix
%D0 = [0.1; 0.1]; % Full covariance matrix 
param0 = [A0; B0; mu0; D0]; % initial parameters 

options = optimoptions('fminunc','MaxFunEvals',25000,'algorithm','quasi-newton', ...
    'TolFun',1e-8,'TolX',1e-8,'MaxIter',1000,'Display','off');

[est_mdl1, params1] = estimate(mdl1, yield, param0, 'Display', 'off', ...
    'options', options, 'Univariate', true);

mu = params1(2*n_factor+1: 3*n_factor)';
[~, ~, ~, ~, ~, ~, ~, deflated_yields1] = DNS(params1, yield, maturity);
deflated_states1 = filter(est_mdl1, deflated_yields1);
estimated_states1 = deflated_states1 + mu;
estimated_yields1 = estimated_states1 * est_mdl1.C';
residual1 = yield - estimated_yields1;

rmse_DNS = sqrt( mean(residual1 .^ 2) ); 

%% DNS_FR model
yield = UK; % in-sample yields
factor = US_factors; 

n_factor = 3; % number of NS factors
n_contract = size(yield, 2); % number of contracts
n_exp = size(factor, 2); % number of explanatory variables 
n_obs = size(yield, 1); % number of observations 

A0 = [3; 3; 3];    
B0 = [0.1; 0.1; 0.1];   
mu0 = [1; 1; 1];
coe0 = repelem(1, n_contract*n_exp)'; 
%D0 = repelem(0.1, n_contract)'; % Diagonal covariance matrix for measurement equation
D0 = [repelem(0.1, n_contract)'; 0.1]; % Band covariance matrix
%D0 = [0.1; 0.1]; % Full covariance matrix
param0 = [A0; B0; mu0; coe0; D0]; % initial parameters

mdl2 = ssm(@(params) DNS_FR(params, yield, factor, maturity));    

options = optimoptions('fminunc','MaxFunEvals',50000,'algorithm','quasi-newton', ...
    'TolFun',1e-8,'TolX',1e-8,'MaxIter',1000,'Display','off');

[est_mdl2, params2] = estimate(mdl2, [yield, factor], param0, 'Display', 'off', ...
    'options', options, 'Univariate', true);

mu = params2(2*n_factor+1: 3*n_factor)';
coe = zeros(n_contract, n_exp);
coe(:) = params2(3*n_factor+1: 3*n_factor+n_contract*n_exp);
[~, ~, ~, ~, ~, ~, ~, deflated_yields2] = DNS_FR(params2, yield, factor, maturity);
[deflated_states2, ~, output2] = filter(est_mdl2, deflated_yields2);
estimated_states2 = deflated_states2 + mu;
estimated_yields2 = estimated_states2 * est_mdl2.C' + factor*coe';
residual2 = yield - estimated_yields2;

rmse2 = sqrt( mean(residual2 .^ 2) ); 

