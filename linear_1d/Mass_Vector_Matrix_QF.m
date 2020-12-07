%% Mass_Vector_Matrix_QF
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
function [int, M, S] = Mass_Vector_Matrix_QF( basis, dx_basis, x_L, x_R, QR )

    if strcmp(QR,'trapez') 
        J = 1000; 
        [x, w] = globTrapez_1d(J);
    elseif strcmp(QR,'Gauss')
        J = 1000; 
        [x, w] = lgwt(J,x_L,x_R);
    else
        error('Desried integration method not implemented yet!')
    end

    %% Vector containing the integrals of the basis functions 
    N = length(basis(0)); 
    F = basis(x');
    int = zeros(N,1); 
    int = F*w;

    %% Mass matrix M 
    M = zeros(N,N); 
    F = basis(x');
    for n=1:N 
        for m=1:N 
            M(n,m) = (F(n,:).*F(m,:))*w; 
        end
    end
    
    %% Stiffnes matrix 
    S = zeros(N,N); 
    F = basis(x'); 
    G = dx_basis(x');
    for n=1:N 
        for m=1:N 
            S(n,m) = (G(n,:).*F(m,:))*w; 
        end
    end
    
end