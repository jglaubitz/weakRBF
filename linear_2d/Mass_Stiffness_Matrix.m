%% Mass_Stiffness_Matrix
% Author: Jan Glaubitz 
% Date: Nov 25, 2020
%
% Construct the mass matrix
%
%  INPUT:
%  basis    : basis functions 
%  dx_basis, dy_basis : partial derivatives of the basis functions 
%  x_L, x_R : boundary points 
%
%  OUTPUT:
%  M : Mass matrix
%  S : Stiffness matrix

%%
function [M, S] = Mass_Stiffness_Matrix( basis, dx_basis, dy_basis, x_L, x_R )
    
    %% Auxilary stuff 
    N = length(basis([0,0]));

    %% Mass matrix M 
    M = zeros(N,N);  
    for j=1:N 
        for k=1:N 
            if k<=j 
                aux_vec1 = zeros(N,1); 
                aux_vec2 = zeros(N,1);
                aux_vec1(j) = 1; 
                aux_vec1(k) = 1; 
                syms x y
                prod = @(x,y) dot( basis([x,y]), aux_vec1 )*... 
                    dot( basis([x,y]), aux_vec2 ); % product b_j*b_k
                M(j,k) = integral2(prod, x_L, x_R, x_L, x_R); 
            else 
                M(j,k) = M(k,j); 
            end
        end
    end
    
    %% Stiffnes matrix 
    S = zeros(N,N); 

    
end