%% Mass_Vector_Matrix
% Author: Jan Glaubitz 
% Date: 07.02.2020
%
% Construct the mass matrix
%
%  INPUT:
%  basis    : basis functions 
%  dx_basis : derivatives of the basis functions 
%  x_L, x_R : boundary points 
%
%  OUTPUT:
%  int    : Vector containing the integrals of the basis functions
%  M      : Mass matrix

%%
function [int, M, S] = Mass_Vector_Matrix( basis, dx_basis, x_L, x_R )
    
    %% Vector containing the integrals of the basis functions 
    N = length(basis(0)); 
    int = zeros(N,1);
    syms x 
    int = integral( @(x) basis(x), x_L, x_R, 'ArrayValued', true );

    %% Mass matrix M 
    M = zeros(N,N); 
    prod = @(x) basis(x)*basis(x)'; % product b_j*b_k 
    M = integral( @(x) prod(x), x_L, x_R, 'ArrayValued', true );
    
    %% Stiffnes matrix 
    S = zeros(N,N); 
    prod = @(x) dx_basis(x)*basis(x)'; 
    S = integral( @(x) prod(x), -1 , 1, 'ArrayValued', true );

    
end