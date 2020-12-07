%% linear_weak_RBF_2d
% Author: Jan Glaubitz 
% Date: Nov 25, 2020
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
function u = linear_weak_RBF_2d( BC, T, CFL, X, u0, kernel, rbf, ep, points, d, integration )
  
N = length(u0); % number of collocation points 

%% Generater basis functions and the differentiation matrix 
DM = Tools_DistanceMatrix(X, X); % Matrix with distances between points
V  = rbf(ep,DM); % Vandermonde matrix of the RBF function
[basis, dx_basis, dy_basis] = Solve_EvaluateBasis(rbf, ep, X, -1, V); % basis 
A = basis(X)'; % new Vandermonde matrix (should be the identity matrix!)
Dx = dx_basis(X(:,1)',X(:,2)')'; % differentiation matrix 

%% Mass, stiffness, and correction matix 
if N==400 && d==0 && ( strcmp(kernel,'cubic') || strcmp(kernel,'quintic') )
    % Mass and stiffness matrix 
    load = matfile(['matrices/M_N=',num2str(N),'_',points,'_',kernel,'_d=',num2str(d),'.mat']);
    M = load.M;
    load = matfile(['matrices/Sx_N=',num2str(N),'_',points,'_',kernel,'_d=',num2str(d),'.mat']);
    Sx = load.Sx; 
    % Correction matrix 
    load = matfile(['matrices/C_N=',num2str(N),'_',points,'_',kernel,'_d=',num2str(d),'_',BC,'.mat']);
    C = load.C;
elseif strcmp(integration,'trapez') 
    % Mass and stiffness matrix 
    [M, Sx] = Mass_Stiffness_Matrix_CF( basis, dx_basis, dy_basis, -1, 1, integration );
    save( ['matrices/M_N=',num2str(N),'_',points,'_',kernel,'_d=',num2str(d),'.mat'], 'M' ); % safe mass matrix 
    save( ['matrices/Sx_N=',num2str(N),'_',points,'_',kernel,'_d=',num2str(d),'.mat'], 'Sx' ); % safe stiffness matrix
    % Correction matrix 
    C = Correction_Matrix_CF( BC, basis, dx_basis, dy_basis, -1, 1, integration ); 
    save( ['matrices/C_N=',num2str(N),'_',points,'_',kernel,'_d=',num2str(d),'_',BC,'.mat'], 'C' ); 
else
    error('Desired integration procedure not implemented yet!')  
end

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
    
    % Third order Runge-Kutta in time 
    % Spacial discretisation: ML(u) = Sx*u - C*u
    % First stage 
    u_t = M\( Sx*u + C*u ); % derivative of u w.r.t 
    u_s1 = u + dt*u_t; % update in time 
    % Second stage 
    u_t = M\( Sx*u_s1 + C*u_s1 ); % derivative of u w.r.t 
    u_s2 = (3/4)*u + (1/4)*u_s1 + (1/4)*dt*u_t;  
    % Third (and final) stage 
    u_t = M\( Sx*u_s2 + C*u_s2 ); % derivative of u w.r.t
    u_s3 = (1/3)*u + (2/3)*u_s2 + (2/3)*dt*u_t; 
    % Update the solution vectors
    u = u_s3; 
    
end

end