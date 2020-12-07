%% compute_trapezoidal_2d
% Generates the vectors of Legender points and weights. 
% 
% INPUT: 
%  n:      number of points in every direction 
%
% OUTPUT: 
%  X: vector of data points
%  w: vector of cubature weights 

function [ X, w] = compute_trapezoidal_2d(n)

    % One dimensional points and weights 
    x = linspace(-1,1,n); % data points 
    Delta_x = 2/(n-1); % distance between any two points 
    w_1d = 2*ones(1,n); % initialize weights 
    w_1d(1) = 1; % correct boundary weights 
    w_1d(end) = 1; % correct boundary weights 
    w_1d = 0.5*Delta_x*w_1d;
    
    % Going over to two dimensions 
    [PointsX, PointsY] = meshgrid(x, x); 
    X = [ reshape(PointsX, numel(PointsX),1)'; 
          reshape(PointsY, numel(PointsY),1)']'; 
    % weights 
    [weightsX, weightsY] = meshgrid(w_1d, w_1d); 
     w_aux = [ reshape(weightsX, numel(weightsX),1)'; 
               reshape(weightsY, numel(weightsY),1)']'; 
     w = w_aux(:,1).*w_aux(:,2);    
           
end

