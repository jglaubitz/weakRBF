%% linear_strong_RBF
% Author: Jan Glaubitz 
% Date: Nov 17, 2020
%
% Solve the linear advection problem by the strong RBF method
%
%  INPUT:
%  BC : boundary condition 
%  T : end time 
%  CFL : CFL number 
%  x : grid points 
%  IC : initial condition 
%  rbf : basis function 
%  ep : shape parameter 
%
%  OUTPUT:
%  u : numerical solution 
%  momentum : values of the momentum over time 
%  energy : values of the energy over time 
%%
function [u, momentum, energy] = linear_strong_RBF( BC, T, CFL, x, IC, rbf, ep )
  
N = length(x); % number of collocation points 
u0 = IC(x); % initial condition 

%% Distance, Vandermone, and differentiation matrix 
D = DifferenceMatrix(x'); % matrix with differences between points
DM = DistanceMatrix(x',x'); % matrix with distances between points 
[V, Ax, D] = Diff_Matrix(rbf, DM, ep, D ); % Vandermonde and (nodal) Differentiation matrices 
[ int, M] = Mass_Matrix(rbf, ep, x'); % vector of integrals and mass matrix 

%% Time step sice 
dx = min( abs( x(2:end) - x(1:end-1) ) ); % spatial step length 
dt = CFL*dx; % time step

%% Time integration - SSPRK(3,3) 
t = 0; % start time 
u = u0; % initialize numerical solution 
u_hat = V\u; % modal values 
momentum(1,1) = t;
momentum(1,2) = dot(int,u_hat);
energy(1,1) = t; 
energy(1,2) = u_hat'*M*u_hat;

j = 1; % counter
while (t<T) 
    
    % time step 
    if T-t<dt 
        dt = T-t; 
    end 
    t = t+dt; 

    % Boundary conditions 
    if strcmp(BC,'inflow') 
        u(1) = IC(1-mod(t,2));
    elseif strcmp(BC,'periodic')
        u(1) = u(N); 
    end
    
    % Third order Runge-Kutta in time
    % First stage 
    u_s1 = u - dt*D*u; 
    % Second stage 
    u_s2 = (3/4)*u + (1/4)*u_s1 - (1/4)*dt*D*u_s1;  
    % Third (and final) stage 
    u_s3 = (1/3)*u + (2/3)*u_s2 - (2/3)*dt*D*u_s2; 
    % Update the solution vectors
    u = u_s3; 
    
    % Momentum and energy over time 
    u_hat = V\u; % modal values
    momentum(j,1) = t;
    momentum(j,2) = dot(int,u_hat);
    energy(j,1) = t; 
    energy(j,2) = u_hat'*M*u_hat;
    j = j+1;
    
end

end