%% Mass_Matrix
% Author: Jan Glaubitz 
% Date: 24.06.2019
%
% Construct the mass matrix
%
%  INPUT:
%  rbf     : annonymous rbf function
%  ep      : rbf shape parameter
%  X       : centres 
%
%  OUTPUT:
%  int    : Vector containing the integrals of the basis functions
%  M      : Mass matrix

%%
function [int, M] = Mass_Matrix(rbf, ep, X )
    
    %% Vector containing the integrals of the basis functions 
    [s, N] = size(X); 
    int = zeros(N,1);
    syms x 
    for k=1:N 
        int(k) = integral( @(x) rbf(ep,abs(x-X(k))), -1 , 1); 
    end

    %% Calculating the mass matrix of the rbf function 
    M = zeros(N,N);
    for k=1:N 
        for n=1:N
            M(k,n) = integral( @(x) rbf( ep, abs(x-X(k)) ).*rbf( ep, abs(x-X(n)) ) , -1 , 1); 
        end
    end

end