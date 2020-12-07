%% grid_points
% Author: Jan Glaubitz 
% Date: Nov 19, 2020 
%
% Construct the grid points
%
%  INPUT:
%  x_L, x_R : boundaries of the domain 
%  N : number of points in every direction 
%  points : type of points 
%
%  OUTPUT:
%  X : grid points 

%%
function [xx, yy, X] = grid_points_2d( x_L, x_R, N, points )
    
    x = linspace(x_L,x_R,N); % start with equidistant points in one dimension  
    [xx, yy] = meshgrid(x, x); 
    X = [ reshape(xx, numel(xx),1)'; 
          reshape(yy, numel(yy),1)']';
    
    if strcmp(points,'random') 
    	for n=1:N^2 
        	if X(n,1) ~= x_L && X(n,1) ~= x_R && X(n,2) ~= x_L && X(n,2) ~= x_R 
                X(n,:) = (x_R-x_L)*rand(1,2)+x_L;
            end
        end
    end
    
end