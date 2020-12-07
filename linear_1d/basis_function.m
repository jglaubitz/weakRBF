%% basis_function
% Author: Jan Glaubitz 
% Date: Nov 17, 2020
%
% Set up the RBF (basis function )
%
%  INPUT:
%  basis : basis function (string)
%
%  OUTPUT:
%  rbf : basis function (function 
%%
function rbf = basis_function( basis )
    
    % radial basis functions 
    rbf_G = @(ep,r) exp(-(ep*r).^2); % Gaussians
    rbf_MQ = @(ep,r) sqrt(1 + (ep*r).^2); % multiquadrics
    rbf_IQ = @(ep,r) 1./(1 + (ep*r).^2); % inverse quadrics 
    rbf_cubic = @(ep,r) (ep*r).^3; % cubic 
    rbf_quintic = @(ep,r) (ep*r).^5; % quintic
    rbf_TPS = @(ep,r) (ep*r).^2 .* log(ep*r); % thin plate spline (TPS)

    % Choose the RBF
    if strcmp(basis,'G')
        rbf = rbf_G; 
    elseif strcmp(basis,'MQ')
        rbf = rbf_MQ; 
    elseif strcmp(basis,'IQ')
        rbf = rbf_IQ;
    elseif strcmp(basis,'cubic')
        rbf = rbf_cubic;
    elseif strcmp(basis,'quintic')
        rbf = rbf_quintic;
    elseif strcmp(basis,'TPS')
        rbf = rbf_TPS;
    else 
        error('desired RBF not yet implemented!')
    end 

end