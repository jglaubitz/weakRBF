%% globTrapez_1d
% Jan Glaubitz, Oct 29, 2020 
%
% Computes the data points, quadrature weights, and result of the gloabal 
% trapezoidal rule in one dimension 
% INPUT
%  fun : underlying piecewise smooth function 
%  N : number of globael data points 
% OUTPUT 
%  x : data points
%  w : quadrature weights 
%  I : result of the corresponding QR 

function [x, w] = globTrapez_1d(N)

    x = linspace(-1,1,N)'; % data points 
    Delta_x = 2/(N-1); % distance between any two points 
    w = 2*ones(N,1); % initialize weights 
    w(1) = 1; % correct boundary weights 
    w(end) = 1; % correct boundary weights 
    w = 0.5*Delta_x*w; % weights 

end