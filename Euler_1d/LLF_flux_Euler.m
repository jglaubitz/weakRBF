%% LLF_flux_Euler
% Author: Jan Glaubitz 
% Date: 25.06.2019
%
% Compute a numerical local Lax-Friedrichs flux for the Euler equations
%
%  INPUT:
%  flux   : flux function
%  a      : LHS values
%  b      : RHS values 
%
%  OUTPUT:
%  fnum   : numerical flux 

%%
function [B] = LLF_flux_euler( flux1, flux2, flux3, g, u1_left, u2_left, u3_left, u1_right, u2_right, u3_right )

    B = zeros(1,3);

    % velocities
    velocity_left = u2_left/u1_left; 
    velocity_right = u2_right/u1_right; 
    % pressure 
    pressure_left = (g-1)*( u3_left - 0.5*(velocity_left^2)*u1_left ); 
    pressure_right = (g-1)*( u3_right - 0.5*(velocity_right^2)*u1_right );
    % sound speed 
    a_left = sqrt( g*pressure_left/u1_left ); 
    a_right = sqrt( g*pressure_right/u1_right ); 
    
    % eigenvalues 
    EV1_left = velocity_left - a_left; 
    EV2_left = velocity_left; 
    EV3_left = velocity_left + a_left; 
    EV1_right = velocity_right - a_right; 
    EV2_right = velocity_right; 
    EV3_right = velocity_right + a_right; 

    % LLF flux 
    lambda = max( abs( [EV1_left; EV2_left; EV3_left; EV1_right; EV2_right; EV3_right] ) );
    B(1,1) = 0.5*( flux1(u1_left,u2_left,u3_left) + flux1(u1_right,u2_right,u3_right) ) - 0.5*lambda*( u1_right - u1_left ); 
    B(1,2) = 0.5*( flux2(u1_left,u2_left,u3_left) + flux2(u1_right,u2_right,u3_right) ) - 0.5*lambda*( u2_right - u2_left ); 
    B(1,3) = 0.5*( flux3(u1_left,u2_left,u3_left) + flux3(u1_right,u2_right,u3_right) ) - 0.5*lambda*( u3_right - u3_left ); 
    
end