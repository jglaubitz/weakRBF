%% comp_num_flux
% Author: Jan Glaubitz 
% Date: 24.06.2019
%
% Compute the numerical flux vector
%
%  INPUT:
%  flux   : flux function 
%  BC     : boundary conditions 
%  IC     : initial condition 
%  u      : solution values
%
%  OUTPUT:
%  fnum   : numerical flux vector

%%
function [C] = comp_num_flux( flux1, flux2, flux3, g, BC, u1, u2, u3 )
    
    N = length(u1); 
    B = zeros(1,3);
    C = zeros(2,3);

    if strcmp(BC,'periodic') 
        % left boundary 
        u1_left = u1(N); 
        u2_left = u2(N); 
        u3_left = u3(N); 
        u1_right = u1(1); 
        u2_right = u2(1); 
        u3_right = u3(1); 
        B = LLF_flux_Euler( flux1, flux2, flux3, g, u1_left, u2_left, u3_left, u1_right, u2_right, u3_right ); 
        C(1,1) = B(1); 
        C(1,2) = B(2);
        C(1,3) = B(3);
        % right boundary 
        u1_left = u1(N); 
        u2_left = u2(N); 
        u3_left = u3(N); 
        u1_right = u1(1); 
        u2_right = u2(1); 
        u3_right = u3(1); 
        B = LLF_flux_Euler( flux1, flux2, flux3, g, u1_left, u2_left, u3_left, u1_right, u2_right, u3_right ); 
        C(2,1) = B(1); 
        C(2,2) = B(2);
        C(2,3) = B(3);
    end 

end