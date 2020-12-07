%% initial_cond
% Author: Jan Glaubitz 
% Date: Nov 17, 2020
%
% Set up the initial condition 
%
%  INPUT:
%  Init_C : initial condition (string)
%
%  OUTPUT:
%  IC : initial condition (function 
%%
function IC = initial_cond( Init_C )
    
    IC_sin = @(x) sin(2*pi*x)+1; % function for the IC 1
    IC_exp = @(x) exp(-20*x.^2); % function for the IC 2
    IC_cos2 = @(x) cos(4*pi*x).^2; % function for the IC 3 

    % select the right IC: sin, exp, disc
    if strcmp(Init_C,'sin')
        IC = IC_sin; 
    elseif strcmp(Init_C,'exp')
        IC = IC_exp; 
    elseif strcmp(Init_C,'cos^2')
        IC = IC_cos2;
    else 
        error('IC not yet implemented!')
    end

end