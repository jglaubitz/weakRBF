%% initial_cond_2d
% Author: Jan Glaubitz 
% Date: Nov 17, 2020
%
% Set up the initial condition and reference solution 
%
%  INPUT:
%  Init_C : initial condition (string)
%  BC : boundary condition (string)
%
%  OUTPUT:
%  IC : initial condition (function) 
%  ref : reference solution (function)
%%
function [IC, ref] = initial_cond_2d( Init_C, BC )
    
    %% Initial condition
    IC_sin = @(x,y) sin(2*pi*x).*(0.5*sin(2*pi*y)-1); % function for the IC 1
    % select the right IC: sin, exp, disc
    if strcmp(Init_C,'sin')
        IC = IC_sin; 
    else 
        error('IC not yet implemented!')
    end

    %% Reference solution
    ref_zeroInflow = @(t,x,y) (x > t-1).*IC(x-t,y); % function for the IC 1 
    ref_periodic = @(t,x,y) IC(x-t,y);
    % select the right IC: sin, exp, disc
    if strcmp(BC,'inflow')
        ref = ref_zeroInflow; 
    elseif strcmp(BC,'periodic') 
        ref = ref_periodic;
    else 
        error('BC not yet implemented!')
    end
    
end