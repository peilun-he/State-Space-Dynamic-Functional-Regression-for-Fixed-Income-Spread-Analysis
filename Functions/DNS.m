function [A, B, C, D, mean0, cov0, state_type, deflated_yield] = DNS(parameter, yield, maturity)
% Parameter-to-matrix map for Dynamic Nelson Siegel model 
%
% For observation vector y(t) and state vector x(t), create a state-space 
% model (SSM) of the form:
% State equation:       x(t) - mu     = A * (x(t-1) - mu) + B * u(t)
% Observation equation: y(t) - C * mu = C * (x(t)  -  mu) + D * e(t)
% Both u(t) and e(t) and multivariate standard normally distributed
%
% 9+T parameters in total: (T = number of contracts)
%   3 for A
%   3 for B
%   3 for mu
%   T for D

n_factor = 3;
n_contract = size(yield, 2);

% Check the length of maturity
if n_contract ~= length(maturity)
    error("The number of contracts and the length of maturity does not match! ");
end

lambda = 0.0609;

A = diag(2 ./ ( 1+exp(-parameter(1: n_factor)) ) - 1);  
B = diag(parameter(n_factor+1: 2*n_factor)); 
C = [ones(size(maturity)), (1 - exp(-lambda*maturity))./(lambda*maturity), ...
    ( (1 - exp(-lambda*maturity))./(lambda*maturity) - exp(-lambda*maturity) ) ];

if length(parameter) == 3*n_factor + n_contract % diagonal covariance matrix
    D = diag(parameter(3*n_factor+1: end));
elseif length(parameter) == 3*n_factor + n_contract + 1 % band covariance matrix
    var_mat = diag(parameter(3*n_factor+1: end-1).^2); % variance matrix
    rho = sqrt(1 + pi^2 / (1 + 4*n_contract^2)) * 1 / (1+exp( -parameter(end) )) - 1/2 * sqrt(1 + pi^2 / (1 + 4*n_contract^2));
    for i = 1: n_contract - 1
        var_mat(i, i+1) = rho * sqrt(var_mat(i, i)) * sqrt(var_mat(i+1, i+1));
        var_mat(i+1, i) = var_mat(i, i+1); 
    end
    D = chol(var_mat)';
elseif length(parameter) == 3*n_factor + 2 % full covariance matrix 
    var_mat = zeros(n_contract, n_contract); 
    rho = 2 / ( 1+exp(-parameter(end)) ) - 1;
    for i = 1: n_contract
        for j = 1: n_contract
            var_mat(i, j) = parameter(end-1)^2 * rho^abs(j-i);
        end
    end
    D = chol(var_mat)';
else
    error("The number of parameters is incorrect! "); 
end


mu = parameter(2*n_factor+1: 3*n_factor);
intercept = C * mu; 

deflated_yield = bsxfun(@minus, yield, intercept'); 

mean0      = [];
cov0       = [];
state_type = [];

end



