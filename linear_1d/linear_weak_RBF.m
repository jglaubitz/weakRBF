%% linear_weak_RBF
% Author: Jan Glaubitz 
% Date: Nov 17, 2020
%
% Solve the linear advection problem by the weak RBF method
%
%  INPUT:
%  BC : boundary condition 
%  T : end time 
%  CFL : CFL number 
%  x : grid points 
%  IC : initial condition 
%  rbf : basis function 
%  ep : shape parameter 
%  d : polynomial degree
%  integration : way integration is performed (exact, trapez, Gauss, LS)
%
%  OUTPUT:
%  u : numerical solution 
%  momentum : values of the momentum over time 
%  energy : values of the energy over time 
%%
function [u, momentum, energy] = linear_weak_RBF( BC, T, CFL, x, IC, rbf, ep, d, integration )
  
N = length(x); % number of collocation points 
u0 = IC(x); % initial condition 

%% Construct RBF basis
DM = DistanceMatrix(x',x'); % matrix with distances between points 
V_rbf = rbf(ep,DM); % Vandermonde matrix of the RBF function
[basis, dx_basis] = Const_Basis( rbf, ep, x', d, V_rbf ); % basis including poynomails

%% Mass and stiffness matix 
if strcmp(integration,'exact') 
    [int, M, S] = Mass_Vector_Matrix( basis, dx_basis, -1, 1 ); % vector of integrals
else
    [int, M, S] = Mass_Vector_Matrix_QF( basis, dx_basis, -1, 1, integration ); % vector of integrals   
end

%% Other matrices 
% Restriction matrix
R = zeros(2,N); % restriction matrix 
R(1,:) = basis(-1)'; 
R(2,:) = basis(1)';
% Boundary matrix 
B = zeros(2,2);
B(1,1) = -1; 
B(2,2) = 1;
% Correction matrix 
C = inv(M)*(R')*B;

%% Time step sice 
dx = min( abs( x(2:end) - x(1:end-1) ) ); % spatial step length 
dt = CFL*dx; % time step

%% Time integration - SSPRK(3,3) 
t = 0; % start time 
u = u0; % initialize numerical solution 
momentum(1,1) = t;
momentum(1,2) = dot(int,u);
energy(1,1) = t; 
energy(1,2) = u'*M*u;

j = 1; % counter
while (t<T) 
    
    % time step 
    if T-t<dt 
        dt = T-t; 
    end 
    t = t+dt; 

    % Boundary conditions (numerical flux)
    if strcmp(BC,'inflow') 
        fnum = [IC(1-mod(t,2));u(N)];
    elseif strcmp(BC,'periodic')
        fnum = [u(N);u(N)]; 
    end
    
    % Third order Runge-Kutta in time 
    % Spacial discretisation: L(u) = inv(M)*S*u - C*fnum
    % First stage 
    u_s1 = u + dt*(inv(M)*S*u - C*fnum); 
    % Second stage 
    u_s2 = (3/4)*u + (1/4)*u_s1 + (1/4)*dt*(inv(M)*S*u_s1 - C*fnum);  
    % Third (and final) stage 
    u_s3 = (1/3)*u + (2/3)*u_s2 + (2/3)*dt*(inv(M)*S*u_s2 - C*fnum); 
    % Update the solution vectors
    u = u_s3; 
    
    % Momentum and energy over time 
    momentum(j,1) = t;
    momentum(j,2) = dot(int,u);
    energy(j,1) = t; 
    energy(j,2) = u'*M*u;
    j = j+1;
    
end

end