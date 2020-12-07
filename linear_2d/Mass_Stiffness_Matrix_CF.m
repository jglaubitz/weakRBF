%% Mass_Stiffness_Matrix_CF
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
%  S      : Stiffness matrix 

%%
function [M, Sx] = Mass_Stiffness_Matrix_CF( basis, dx_basis, dy_basis, x_L, x_R, integration )

    if strcmp(integration,'trapez') 
        J = 1000; 
        [ Y, w] = compute_trapezoidal_2d(J);
    else
        error('Desried integration method not implemented yet!')
    end

    %% Mass matrix M and stiffnes matrix Sx (w.r.t x)
    N = length(basis([0,0]));
    M = zeros(N,N); 
    Sx = zeros(N,N);
    F = basis(Y); 
    G = dx_basis(Y(:,1)',Y(:,2)');
    for n=1:N 
        for m=1:N 
            if n<=m
                M(n,m) = (F(n,:).*F(m,:))*w; 
            else 
                M(n,m) = M(m,n); 
            end 
            Sx(n,m) = (G(n,:).*F(m,:))*w;
            [m,n,N]
        end
    end
    
end