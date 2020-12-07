%% Const_Basis
% Author: Jan Glaubitz 
% Date: 06.11.2018
%
% Construct a basis for V_{N,d}
%
%  INPUT:
%  rbf     : annonymous rbf function
%  ep      : rbf shape parameter
%  X       : set of centers 
%  d       : polynomial degree
%
%  OUTPUT:
%  basis       : basis of V_{N,d}

%%
function [basis, dx_basis] = Const_Basis( rbf, ep, X, d, V_rbf )
    
    %% no polynomials included
    if d < 0 
        
        %% coefficients 
        N = length(X); % number of centers
        gamma = zeros(N,N);
        y = eye(N); 
        gamma = V_rbf\y; 
       
        %% New basis {b_j} 
        basis = @(x) (gamma')*rbf(ep,abs(x-X')); 
        
        %% Derivative 
        syms x 
        dx_basis = matlabFunction( diff( basis(x) , x ) ); 
        
    %% polynomials included     
    else
    
        % Vandermonde matrix 
        %dis = DifferenceMatrix(X); % matrix with differences between points
        %DM = DistanceMatrix(X,X); % matrix with distances between points
        %[V, Ax, D] = Diff_Matrix(rbf, DM, ep, dis ); % Vandermonde matrix
        
        % define the polynomial basis p_0,...,p_d 
        if d == 0
            p = @(x) x.^0; 
        elseif d == 1
            p = @(x) [x.^0 ; x]; 
        elseif d == 2
            p = @(x) [x.^0 ; x ; x.^2]; 
        elseif d == 3
            p = @(x) [x.^0 ; x ; x.^2 ; x.^3]; 
        else 
            'd to large';
        end
        
        % Polynomial matrix P 
        N = length(X); % number of centers
        P = p(X)'; % polynomial matrix
        
        % Matrix A = ( V P ; P^T 0 )
        A = [ V_rbf P; P' zeros(d+1,d+1) ];
        
        % Coefficients alpha_j, beta_j 
        alpha = zeros(N,N); 
        beta = zeros(d+1,N); 
        gamma = zeros(N+d+1,N);
        y = zeros(N+d+1,N); 
        y(1:N,1:N) = eye(N); 
        gamma = A\y; 
        alpha = gamma(1:N,:); 
        beta = gamma(N+1:N+d+1,:);
       
        %% New basis {b_j} 
        basis = @(x) (alpha')*rbf(ep,abs(x-X')) + (beta')*p(x); 
        
        %% Derivative 
        syms x 
        dx_basis = matlabFunction( diff( basis(x) , x ) ); 
        
    end
    
end