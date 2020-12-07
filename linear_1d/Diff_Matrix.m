%% Diff_Matrix
% Author: Jan Glaubitz 
% Date: 06.11.2018
%
% Construct the Vandermonde and differentiation matrices for the *regular/strong* RBF
% method.
%
%  INPUT:
%  rbf     : annonymous rbf function
%  DM      : distance matrix
%  ep      : rbf shape parameter
%  d       : difference matrix 
%
%  OUTPUT:
%  V      : Vandermonde matrix 
%  Ax     : RBF differentiation matrix
%  D      : Nodal differentiation matrix

%%
function [V, Ax, D] = Diff_Matrix(rbf, DM, ep, d )
    
    %% Vandermonde matrix of the rbf function
    V = rbf(ep,DM); % Vandermonde matrix 

    %% Calculating the derivative of the rbf function
    syms e R x
    dxrbf = matlabFunction(diff(rbf(e,R),R) * x/R, 'vars',[e R x]);
    
    %% Creating the support matrix Ax to calculate D
    Ax = dxrbf(ep,DM,d);
    
    %% Differentiation matrix D
    D = Ax*inv(V);

end