%% Correction_Matrix_CF
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
function C = Correction_Matrix_CF( BC, basis, dx_basis, dy_basis, x_L, x_R, integration )

    if strcmp(integration,'trapez') 
        J = 1000; 
        [x, w] = globTrapez_1d(J);
    else
        error('Desried integration method not implemented yet!')
    end

    %% Mass matrix M and stiffnes matrix Sx (w.r.t x)
    N = length(basis([0,0]));
    C = zeros(N,N); 
    % periodic BC
    if strcmp(BC,'periodic')
        F =  basis([-ones(length(x),1),x]) - basis([ones(length(x),1),x]); 
        G = basis([ones(length(x),1),x]);
        for n=1:N 
            for m=1:N 
                C(n,m) = (F(n,:).*G(m,:))*w; 
                [m,n,N]
            end
        end
    % else     
    else
        error('Desired BC not yet implemented!')
    end
    
end