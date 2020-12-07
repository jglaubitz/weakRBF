%% linear_strong_RBF_2d
% Author: Jan Glaubitz 
% Date: Nov 25, 2020
%
% Solve the linear advection problem by the strong RBF method
%
%  INPUT:
%  BC : boundary condition 
%  T : end time 
%  CFL : CFL number 
%  X : grid points 
%  IC : initial condition 
%  rbf : basis function 
%  ep : shape parameter 
%
%  OUTPUT:
%  u : numerical solution  
%%
function u = linear_strong_RBF_2d( BC, T, CFL, X, u0, rbf, ep )
  
N = length(u0); % number of collocation points 

%% Generater basis functions and the differentiation matrix 
DM = Tools_DistanceMatrix(X, X); % Matrix with distances between points
V  = rbf(ep,DM); % Vandermonde matrix of the RBF function
[basis, dx_basis, dy_basis] = Solve_EvaluateBasis(rbf, ep, X, -1, V); % basis 
A = basis(X)'; % new Vandermonde matrix (should be the identity matrix!)
Dx = dx_basis(X(:,1)',X(:,2)')'; % differentiation matrix 

%% Time step sice 
dx = min( DM+42*eye(N), [], 'all' ); % spatial step length 
dt = CFL*dx; % time step

%% Time integration - SSPRK(3,3) 
t = 0; % start time 
u = u0; % initialize numerical solution 

while (t<T) 
    
    % time step 
    if T-t<dt 
        dt = T-t; 
    end 
    t = t+dt; 

    % BC
    if strcmp(BC,'inflow') 
        x_L = min(X(:,1)); % left boundary
        u( X(:,1) == x_L ) = 0; % set all values at the west boundary to zero
    elseif strcmp(BC,'periodic') 
        x_L = min(X(:,1)); % left boundary 
        x_R = max(X(:,1)); % right boundary 
        u( X(:,1) == x_L ) = u( X(:,1) == x_R ); % periodic BC
    else
        error('Desired BC not implemented yet')
    end
    
    % Third order Runge-Kutta in time
    % First stage 
    u_s1 = u - dt*Dx*u; 
    % Second stage 
    u_s2 = (3/4)*u + (1/4)*u_s1 - (1/4)*dt*Dx*u_s1;  
    % Third (and final) stage 
    u_s3 = (1/3)*u + (2/3)*u_s2 - (2/3)*dt*Dx*u_s2; 
    % Update the solution vectors
    u = u_s3; 
    
end

end