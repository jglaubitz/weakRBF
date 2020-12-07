%% Diff_Matrix
% Author: Jan Glaubitz 
% Date: 24.06.2019
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
%  S       : stiffness matrix 

%%
function S = Stiffness_Matrix(rbf, ep, X )
    
    [s, N] = size(X);

    %% Calculating the derivative of the rbf function
    syms e R x
    dxrbf = matlabFunction(diff(rbf(e,R),R) * x/R, 'vars',[e R x]);
    
    %% Creating the support matrix Ax to calculate D
    %Ax = dxrbf(ep,DM,d);
    
    %% Calculating the mass matrix of the rbf function 
    S = zeros(N,N);
    for k=1:N 
        for n=1:N
            S(k,n) = integral( @(x) rbf( ep, abs(x-X(n)) ).*dxrbf( ep, abs( x-X(k) ), x-X(k) ) , -1 , 1); 
        end
    end

end