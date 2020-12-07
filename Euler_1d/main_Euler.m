%% Strong RBF collocation method for the Euler equations 
% Author: Jan Glaubitz 
% Date: 24.06.2019

% Domain is [0,1]

%%
clear, clc, close all

%% Setting up common variables 
Init_C = 'smooth_flow'; % smooth_flow
BC = 'periodic'; % inflow, periodic
T = 0.1; % final time 
basis = 'cubic'; % G, MQ, IQ, cubic, quintic, TPS 
ep = 20; % shape parameter
N = 20; % number of points 
P = 2; % degree up to which polynomials are included 

% set up Euler equations 
g = 3; % ratio of specific heats 
flux1 = @(u1,u2,u3) u2; % 1st flux component in the Euler equations 
flux2 = @(u1,u2,u3) 0.5*(3-g)*u2.^2./u1 + (g-1)*u3; % 2nd flux component in the Euler equations 
flux3 = @(u1,u2,u3) g*u2.*u3./u1 - 0.5*(g-1)*u2.^3./(u1.^2); % 3th flux component in the Euler equations 

%% radial basis functions 
rbf_G = @(ep,r) exp(-(ep*r).^2); % Gaussians
rbf_MQ = @(ep,r) sqrt(1 + (ep*r).^2); % multiquadrics
rbf_IQ = @(ep,r) 1./(1 + (ep*r).^2); % inverse quadrics 
rbf_cubic = @(ep,r) (ep*r).^3; % cubic 
rbf_quintic = @(ep,r) (ep*r).^5; % quintic
rbf_TPS = @(ep,r) (ep*r).^2 .* log(ep*r); % thin plate spline (TPS)
% Choose the RBF
if strcmp(basis,'G')
    rbf = rbf_G; 
elseif strcmp(basis,'MQ')
    rbf = rbf_MQ; 
elseif strcmp(basis,'IQ')
    rbf = rbf_IQ;
elseif strcmp(basis,'cubic')
    rbf = rbf_cubic;
elseif strcmp(basis,'quintic')
    rbf = rbf_quintic;
elseif strcmp(basis,'TPS')
    rbf = rbf_TPS;
else 
    fprintf('Wrong RBF!')
end 

%% Generating the collocation points 
X = linspace(-1,1,N)'; % equidistant collocation points 
dx = max( abs( X(2:N) - X(1:N-1) ) ); % Spatial step length 
dt = 0.1*dx; % Time step
XX = linspace(-1,1,10*N)'; % evaluation points 

%% Distance, Vandermone, and differentiation matrix 
d = DifferenceMatrix(X,X); % matrix with differences between points
DM = DistanceMatrix(X,X); % matrix with distances between points 
DM_XX = DistanceMatrix(X,XX); % distance matrix for the extrapolation matrix
[V, Ax, D] = Diff_Matrix(rbf, DM, ep, d ); % Vandermonde and (nodal) Differentiation matrices 
[ int, M] = Mass_Matrix(rbf, ep, X'); % vector of integrals and mass matrix 
Extrapol_matrix = rbf(ep,DM_XX'); % extrapolation matrix

%% Additional matrices for the weak RBF method (P=0)
S = Stiffness_Matrix(rbf, ep, X'); % stiffness matrix
R = zeros(2,N); % restriction matrix 
R(1,:) = V(1,:); 
R(2,:) = V(N,:);
B = zeros(2,2); % boundary matrix 
B(1,1) = -1; 
B(2,2) = 1; 
C = inv(M)*(R')*B; % correction matrix 

%% Matrices for the the weak RBF method with P>0
% Matrix for the polynomials and the RBF approx. including polynomials 
Q = zeros(N,P); 
if P == 1
    q = @(x) x.^0; 
elseif P == 2
    q = @(x) [ x.^0; x.^1]; 
elseif P == 3 
    q = @(x) [ x.^0; x.^1; x.^2 ]; 
elseif P == 4 
    q = @(x) [ x.^0; x.^1; x.^2; x.^3 ]; 
else 
    'P to large';
end
Q = q(X')';
A = [ V Q; Q' zeros(P,P) ];
% Coefficients alpha_n, beta_n 
alpha = zeros(N,N); 
beta = zeros(P,N); 
gamma = zeros(N+P,N);
y = zeros(N+P,N); 
y(1:N,1:N) = eye(N); 
gamma = A\y; 
alpha = gamma(1:N,:); 
beta = gamma(N+1:N+P,:);
% New basis {b_n} 
basis = @(x) (alpha')*rbf(ep,abs(x-X)) + (beta')*q(x); 
% Vandermonde matrix w.r.t. {b_n}_{n=1}^N
Vand_b = basis(X')';
% Compute the vector of integrals 
int_P1 = zeros(N,1); 
int_P1 = integral( @(x) basis(x), -1, 1, 'ArrayValued', true );
% Compute the mass matrix 
b_b = @(x) basis(x)*basis(x)'; 
M_P1 = zeros(N,N); 
M_P1 = integral( @(x) b_b(x), -1, 1, 'ArrayValued', true );
% Compute the stiffness matrix 
S_P1 = zeros(N,N); 
syms x
dxb = matlabFunction( diff( basis(x) , x ) );
b_dxb = @(x) dxb(x)*basis(x)'; 
S_P1 = integral( @(x) b_dxb(x), -1, 1, 'ArrayValued', true );
% Restriction matrix
R_P1 = zeros(2,N); % restriction matrix 
R_P1(1,:) = basis(-1)'; 
R_P1(2,:) = basis(1)';
% Boundary matrix 
B = zeros(2,2);
B(1,1) = -1; 
B(2,2) = 1;
% Correction matrix 
C_P1 = inv(M_P1)*(R_P1')*B; 
% Extrapolation matrix 
Extrapol_matrix_P1 = basis(XX')';

%% Generate initial condition 
IC_density = @(x) 1 + 0.5*sin(pi*x); % IC for the density 
IC_velocity = @(x) 0; % IC for the velocity 
IC_pressure = @(x) IC_density(x).^g; % IC for the pressure 
IC_energy = @(x) IC_pressure(x)/(g-1) + 0.5*(IC_velocity(x).^2).*IC_density(x); % energy is derived from the other ICs
% nodal values of the IC 
u1_IC = IC_density(X); % IC for u1 at the centres X 
u2_IC = IC_density(X).*IC_velocity(X); % IC for u2 at the centres X 
u3_IC = IC_pressure(X)/(g-1) + 0.5*(IC_velocity(X).^2).*IC_density(X); % IC for u2 at the centres X 
%% strong RBF method
u1_s = u1_IC; % nodal values 1st component
u2_s = u2_IC; % nodal values 2nd component 
u3_s = u3_IC; % nodal values 3th component 
u1_hat_s = V\u1_s; % modal coefficients 1st component
u2_hat_s = V\u2_s; % nodal values 2nd component  
u3_hat_s = V\u3_s; % nodal values 3th component 
%% weak RBF method P=0 
u1_P0 = u1_IC; % nodal values 1st component 
u2_P0 = u2_IC; % nodal values values 2nd component
u3_P0 = u3_IC; % nodal values 3th component
u1_hat_P0 = V\u1_P0; % modal coefficients 1st component
u2_hat_P0 = V\u2_P0; % modal coefficients values 2nd component 
u3_hat_P0 = V\u3_P0; % modal coefficients 3th component 
fnum1_P0 = zeros(2,1); % numerical flux 1st component 
fnum2_P0 = zeros(2,1); % numerical flux values 2nd component 
fnum3_P0 = zeros(2,1); % numerical flux 3th component 
%% weak RBF method P=1 
u1_P1 = u1_IC; % nodal values 1st component 
u2_P1 = u2_IC; % nodal values values 2nd component
u3_P1 = u3_IC; % nodal values 3th component
u1_hat_P1 = Vand_b\u1_P1; % modal coefficients 1st component 
u2_hat_P1 = Vand_b\u2_P1; % modal coefficients values 2nd component 
u3_hat_P1 = Vand_b\u3_P1; % modal coefficients 3th component 
fnum1_P1 = zeros(2,1); % numerical flux 1st component 
fnum2_P1 = zeros(2,1); % numerical flux values 2nd component 
fnum3_P1 = zeros(2,1); % numerical flux 3th component 

%% Time integration - SSPRK(3,3) 
t = 0; % start time
while (t<T) 
    if T-t<dt 
        dt = T-t; 
    end 
    t = t+dt; 
    
    %% Boundary conditions 
    % strong RBF method
    mean = 0.5*( u1_s(1) + u1_s(N) ); 
    u1_s(1) = mean; 
    u1_s(N) = mean; 
    mean = 0.5*( u2_s(1) + u2_s(N) ); 
    u2_s(1) = mean; 
    u2_s(N) = mean; 
    mean = 0.5*( u3_s(1) + u3_s(N) ); 
    u3_s(1) = mean; 
    u3_s(N) = mean; 
    % weak RBF method P=0 
    F = comp_num_flux( flux1, flux2, flux3, g, BC, u1_P0, u2_P0, u3_P0 );
    fnum1_P0 = F(:,1); 
    fnum2_P0 = F(:,2);
    fnum3_P0 = F(:,3);
    % weak RBF method P=1 
    F = comp_num_flux( flux1, flux2, flux3, g, BC, u1_P1, u2_P1, u3_P1 );
    fnum1_P1 = F(:,1); 
    fnum2_P1 = F(:,2);
    fnum3_P1 = F(:,3);
    
    %% Third order Runge-Kutta in time - strong RBF method 
    % 1st stage 
    f1_s = flux1(u1_s,u2_s,u3_s); % flux values 
    f2_s = flux2(u1_s,u2_s,u3_s); % flux values 
    f3_s = flux2(u1_s,u2_s,u3_s); % flux values 
    u1_s_s1 = u1_s - dt*D*f1_s; % update 
    u2_s_s1 = u2_s - dt*D*f2_s; % update 
    u3_s_s1 = u3_s - dt*D*f3_s; % update 
    % 2nd stage 
    f1_s = flux1(u1_s_s1,u2_s_s1,u3_s_s1); % flux values 
    f2_s = flux2(u1_s_s1,u2_s_s1,u3_s_s1); % flux values 
    f3_s = flux2(u1_s_s1,u2_s_s1,u3_s_s1); % flux values 
    u1_s_s2 = (3/4)*u1_s + (1/4)*u1_s_s1 - (1/4)*dt*D*f1_s; % update 
    u2_s_s2 = (3/4)*u2_s + (1/4)*u2_s_s1 - (1/4)*dt*D*f2_s; % update 
    u3_s_s2 = (3/4)*u3_s + (1/4)*u3_s_s1 - (1/4)*dt*D*f3_s; % update 
    % 3rd (and final) stage 
    f1_s = flux1(u1_s_s2,u2_s_s2,u3_s_s2); % flux values 
    f2_s = flux2(u1_s_s2,u2_s_s2,u3_s_s2); % flux values 
    f3_s = flux2(u1_s_s2,u2_s_s2,u3_s_s2); % flux values 
    u1_s_s3 = (1/3)*u1_s + (2/3)*u1_s_s2 - (2/3)*dt*D*f1_s; % update 
    u2_s_s3 = (1/3)*u2_s + (2/3)*u2_s_s2 - (2/3)*dt*D*f2_s; % update 
    u3_s_s3 = (1/3)*u3_s + (2/3)*u3_s_s2 - (2/3)*dt*D*f3_s; % update 
    % Update the solution vectors
    u1_s = u1_s_s3; % new solution values 
    u2_s = u2_s_s3; % new solution values 
    u3_s = u3_s_s3; % new solution values 
    
    %% Third order Runge-Kutta in time - weak RBF method P=0 
    % 1st stage 
    f1_P0 = flux1(u1_P0,u2_P0,u3_P0); % flux values
    f2_P0 = flux2(u1_P0,u2_P0,u3_P0); % flux values
    f3_P0 = flux3(u1_P0,u2_P0,u3_P0); % flux values
    f1_hat_P0 = V\f1_P0; % flux coefficients 
    f2_hat_P0 = V\f2_P0; % flux coefficients 
    f3_hat_P0 = V\f3_P0; % flux coefficients 
    u1_hat_P0_s1 = u1_hat_P0 + dt*( inv(M)*S*f1_hat_P0 - C*fnum1_P0 ); % update 
    u2_hat_P0_s1 = u2_hat_P0 + dt*( inv(M)*S*f2_hat_P0 - C*fnum2_P0 ); % update 
    u3_hat_P0_s1 = u3_hat_P0 + dt*( inv(M)*S*f3_hat_P0 - C*fnum3_P0 ); % update 
    % 2nd stage 
    u1_P0_s1 = V*u1_hat_P0_s1; % nodal values 
    u2_P0_s1 = V*u2_hat_P0_s1; % nodal values 
    u3_P0_s1 = V*u3_hat_P0_s1; % nodal values 
    f1_P0 = flux1(u1_P0_s1,u2_P0_s1,u3_P0_s1); % flux values
    f2_P0 = flux2(u1_P0_s1,u2_P0_s1,u3_P0_s1); % flux values
    f3_P0 = flux3(u1_P0_s1,u2_P0_s1,u3_P0_s1); % flux values
    f1_hat_P0 = V\f1_P0; % flux coefficients 
    f2_hat_P0 = V\f2_P0; % flux coefficients 
    f3_hat_P0 = V\f3_P0; % flux coefficients 
    u1_hat_P0_s2 = (3/4)*u1_hat_P0 + (1/4)*u1_hat_P0_s1 + (1/4)*dt*( inv(M)*S*f1_hat_P0 - C*fnum1_P0 ); % update
    u2_hat_P0_s2 = (3/4)*u2_hat_P0 + (1/4)*u2_hat_P0_s1 + (1/4)*dt*( inv(M)*S*f2_hat_P0 - C*fnum2_P0 ); % update
    u3_hat_P0_s2 = (3/4)*u3_hat_P0 + (1/4)*u3_hat_P0_s1 + (1/4)*dt*( inv(M)*S*f3_hat_P0 - C*fnum3_P0 ); % update
    % 3th (and final) stage 
    u1_P0_s2 = V*u1_hat_P0_s2; % nodal values 
    u2_P0_s2 = V*u2_hat_P0_s2; % nodal values 
    u3_P0_s2 = V*u3_hat_P0_s2; % nodal values 
    f1_P0 = flux1(u1_P0_s2,u2_P0_s2,u3_P0_s2); % flux values
    f2_P0 = flux2(u1_P0_s2,u2_P0_s2,u3_P0_s2); % flux values
    f3_P0 = flux3(u1_P0_s2,u2_P0_s2,u3_P0_s2); % flux values
    f1_hat_P0 = V\f1_P0; % flux coefficients 
    f2_hat_P0 = V\f2_P0; % flux coefficients 
    f3_hat_P0 = V\f3_P0; % flux coefficients 
    u1_hat_P0_s3 = (1/3)*u1_hat_P0 + (2/3)*u1_hat_P0_s2 + (2/3)*dt*( inv(M)*S*f1_hat_P0 - C*fnum1_P0 ); % update
    u2_hat_P0_s3 = (1/3)*u2_hat_P0 + (2/3)*u2_hat_P0_s2 + (2/3)*dt*( inv(M)*S*f2_hat_P0 - C*fnum2_P0 ); % update
    u3_hat_P0_s3 = (1/3)*u3_hat_P0 + (2/3)*u3_hat_P0_s2 + (2/3)*dt*( inv(M)*S*f3_hat_P0 - C*fnum3_P0 ); % update
    % Update the solution vectors
    u1_hat_P0 = u1_hat_P0_s3; % new modal coefficients 
    u2_hat_P0 = u2_hat_P0_s3; % new modal coefficients 
    u3_hat_P0 = u3_hat_P0_s3; % new modal coefficients 
    u1_P0 = V*u1_hat_P0; % new (nodal) solution values 
    u2_P0 = V*u2_hat_P0; % new (nodal) solution values 
    u3_P0 = V*u3_hat_P0; % new (nodal) solution values 
    
    %% Third order Runge-Kutta in time - weak RBF method P=1 
    % 1st stage 
    f1_P1 = flux1(u1_P1,u2_P1,u3_P1); % flux values
    f2_P1 = flux2(u1_P1,u2_P1,u3_P1); % flux values
    f3_P1 = flux3(u1_P1,u2_P1,u3_P1); % flux values
    f1_hat_P1 = Vand_b\f1_P1; % flux coefficients 
    f2_hat_P1 = Vand_b\f2_P1; % flux coefficients 
    f3_hat_P1 = Vand_b\f3_P1; % flux coefficients 
    u1_hat_P1_s1 = u1_hat_P1 + dt*( inv(M_P1)*S_P1*f1_hat_P1 - C_P1*fnum1_P1 ); % update 
    u2_hat_P1_s1 = u2_hat_P1 + dt*( inv(M_P1)*S_P1*f2_hat_P1 - C_P1*fnum2_P1 ); % update 
    u3_hat_P1_s1 = u3_hat_P1 + dt*( inv(M_P1)*S_P1*f3_hat_P1 - C_P1*fnum3_P1 ); % update 
    % 2nd stage 
    u1_P1_s1 = Vand_b*u1_hat_P1_s1; % nodal values 
    u2_P1_s1 = Vand_b*u2_hat_P1_s1; % nodal values 
    u3_P1_s1 = Vand_b*u3_hat_P1_s1; % nodal values 
    f1_P1 = flux1(u1_P1_s1,u2_P1_s1,u3_P1_s1); % flux values
    f2_P1 = flux2(u1_P1_s1,u2_P1_s1,u3_P1_s1); % flux values
    f3_P1 = flux3(u1_P1_s1,u2_P1_s1,u3_P1_s1); % flux values
    f1_hat_P1 = Vand_b\f1_P1; % flux coefficients 
    f2_hat_P1 = Vand_b\f2_P1; % flux coefficients 
    f3_hat_P1 = Vand_b\f3_P1; % flux coefficients 
    u1_hat_P1_s2 = (3/4)*u1_hat_P1 + (1/4)*u1_hat_P1_s1 + (1/4)*dt*( inv(M_P1)*S_P1*f1_hat_P1 - C_P1*fnum1_P1 ); % update
    u2_hat_P1_s2 = (3/4)*u2_hat_P1 + (1/4)*u2_hat_P1_s1 + (1/4)*dt*( inv(M_P1)*S_P1*f2_hat_P1 - C_P1*fnum2_P1 ); % update
    u3_hat_P1_s2 = (3/4)*u3_hat_P1 + (1/4)*u3_hat_P1_s1 + (1/4)*dt*( inv(M_P1)*S_P1*f3_hat_P1 - C_P1*fnum3_P1 ); % update
    % 3th (and final) stage 
    u1_P1_s2 = Vand_b*u1_hat_P1_s2; % nodal values 
    u2_P1_s2 = Vand_b*u2_hat_P1_s2; % nodal values 
    u3_P1_s2 = Vand_b*u3_hat_P1_s2; % nodal values 
    f1_P1 = flux1(u1_P1_s2,u2_P1_s2,u3_P1_s2); % flux values
    f2_P1 = flux2(u1_P1_s2,u2_P1_s2,u3_P1_s2); % flux values
    f3_P1 = flux3(u1_P1_s2,u2_P1_s2,u3_P1_s2); % flux values
    f1_hat_P1 = Vand_b\f1_P1; % flux coefficients 
    f2_hat_P1 = Vand_b\f2_P1; % flux coefficients 
    f3_hat_P1 = Vand_b\f3_P1; % flux coefficients 
    u1_hat_P1_s3 = (1/3)*u1_hat_P1 + (2/3)*u1_hat_P1_s2 + (2/3)*dt*( inv(M_P1)*S_P1*f1_hat_P1 - C_P1*fnum1_P1 ); % update
    u2_hat_P1_s3 = (1/3)*u2_hat_P1 + (2/3)*u2_hat_P1_s2 + (2/3)*dt*( inv(M_P1)*S_P1*f2_hat_P1 - C_P1*fnum2_P1 ); % update
    u3_hat_P1_s3 = (1/3)*u3_hat_P1 + (2/3)*u3_hat_P1_s2 + (2/3)*dt*( inv(M_P1)*S_P1*f3_hat_P1 - C_P1*fnum3_P1 ); % update
    % Update the solution vectors
    u1_hat_P1 = u1_hat_P1_s3; % new modal coefficients 
    u2_hat_P1 = u2_hat_P1_s3; % new modal coefficients 
    u3_hat_P1 = u3_hat_P1_s3; % new modal coefficients 
    u1_P1 = Vand_b*u1_hat_P1; % new (nodal) solution values 
    u2_P1 = Vand_b*u2_hat_P1; % new (nodal) solution values 
    u3_P1 = Vand_b*u3_hat_P1; % new (nodal) solution values 
    
end

% Numerical solution at the evaluation points 
% strong RBF method 
u1_hat_s = V\u1_s; % modal coefficients 
u2_hat_s = V\u2_s; % modal coefficients 
u3_hat_s = V\u3_s; % modal coefficients 
uu1_s = Extrapol_matrix*u1_hat_s; % nodal values at evaluation points 
uu2_s = Extrapol_matrix*u2_hat_s; % nodal values at evaluation points 
uu3_s = Extrapol_matrix*u3_hat_s; % nodal values at evaluation points 
% weak RBF method P=0 
uu1_P0 = Extrapol_matrix*u1_hat_P0; % nodal values at evaluation points 
uu2_P0 = Extrapol_matrix*u2_hat_P0; % nodal values at evaluation points 
uu3_P0 = Extrapol_matrix*u3_hat_P0; % nodal values at evaluation points 
% weak RBF method P=1 
uu1_P1 = Extrapol_matrix_P1*u1_hat_P1; % nodal values at evaluation points 
uu2_P1 = Extrapol_matrix_P1*u2_hat_P1; % nodal values at evaluation points 
uu3_P1 = Extrapol_matrix_P1*u3_hat_P1; % nodal values at evaluation points 

% Physical values 
density_s = u1_s; % density at centers
density_P0 = u1_P0; % density at centers
density_P1 = u1_P1; % density at centers 
ddensity_s = uu1_s; % density at evaluation points 
ddensity_P0 = uu1_P0; % density at evaluation points 
ddensity_P1 = uu1_P1; % density at evaluation points 
velocity_s = u2_s./u1_s; % velocity at centers
velocity_P0 = u2_P0./u1_P0; % velocity at centers
velocity_P1 = u2_P1./u1_P1; % velocity at centers
vvelocity_s = uu2_s./uu1_s; % velocity at evaluation points 
vvelocity_P0 = uu2_P0./uu1_P0; % velocity at evaluation points 
vvelocity_P1 = uu2_P1./uu1_P1; % velocity at evaluation points 
energy_s = u3_s; % energy at centers
energy_P0 = u3_P0; % energy at centers
energy_P1 = u3_P1; % energy at centers
eenergy_s = uu3_s; % energy at evaluation points 
eenergy_P0 = uu3_P0; % energy at evaluation points 
eenergy_P1 = uu3_P1; % energy at evaluation points 
pressure_s = (g-1)*( u3_s - 0.5*(u2_s.^2)./u1_s ); % pressure at centers
pressure_P0 = (g-1)*( u3_P0 - 0.5*(u2_P0.^2)./u1_P0 ); % pressure at centers
pressure_P1 = (g-1)*( u3_P1 - 0.5*(u2_P1.^2)./u1_P1 ); % pressure at centers
ppressure_s = (g-1)*( uu3_s - 0.5*(uu2_s.^2)./uu1_s ); % pressure at evaluation points 
ppressure_P0 = (g-1)*( uu3_P0 - 0.5*(uu2_P0.^2)./uu1_P0 ); % pressure at evaluation points 
ppressure_P1 = (g-1)*( uu3_P1 - 0.5*(uu2_P1.^2)./uu1_P1 ); % pressure at evaluation points 

% Compute physical reference values 
density_ref = zeros(N,1); % density at centers
velocity_ref = zeros(N,1); % velocity at centers
pressure_ref = zeros(N,1); % pressure at centers 
energy_ref = zeros(N,1); % internal energy at centers
% compute 
myfun1 = @(x,y) x + sqrt(3)*IC_density(y)*T - y; % parameterized function 
myfun2 = @(x,y) x - sqrt(3)*IC_density(y)*T - y; % parameterized function
for n=1:N 
    x = X(n); % parameter
    fun1 = @(y) myfun1(x,y); % function of y alone 
    fun2 = @(y) myfun2(x,y); % function of y alone 
    x1 = fzero( fun1, x ); % x1 
    x2 = fzero( fun2, x ); % x2
    density_ref(n) = 0.5*( IC_density(x1) + IC_density(x2) ); % density 
    velocity_ref(n) = sqrt(3)*( density_ref(n) - IC_density(x1) ); % velocity 
    pressure_ref(n) = density_ref(n).^g; % pressure (isentropic law)
    energy_ref(n) = pressure_ref(n)/(g-1) + 0.5*(velocity_ref(n).^2)*density_ref(n); % internal energy (equation of state)
end 
ddensity_ref = zeros(10*N,1); % density at evaluation points 
vvelocity_ref = zeros(10*N,1); % velocity at evaluation points
ppressure_ref = zeros(10*N,1); % pressure at evaluation points 
eenergy_ref = zeros(10*N,1); % internal energy at evaluation points
for n=1:10*N 
    x = XX(n); % parameter
    fun1 = @(y) myfun1(x,y); % function of y alone 
    fun2 = @(y) myfun2(x,y); % function of y alone 
    x1 = fzero( fun1, x ); % x1 
    x2 = fzero( fun2, x ); % x2
    ddensity_ref(n) = 0.5*( IC_density(x1) + IC_density(x2) ); % density 
    vvelocity_ref(n) = sqrt(3)*( ddensity_ref(n) - IC_density(x1) ); % velocity 
    ppressure_ref(n) = ddensity_ref(n).^g; % pressure 
    eenergy_ref(n) = ppressure_ref(n)/(g-1) + 0.5*(vvelocity_ref(n).^2)*ddensity_ref(n); % internal energy 
end 

% plot density
figure(1) 
p = plot(XX,ddensity_ref,'k:', 'LineWidth',2.5); 
hold on 
sz = 80;
scatter(X,density_s,sz,'rs', 'filled'); 
scatter(X,density_P1,sz,'bo', 'filled');
hold off
set(gca, 'FontSize', 24)  % Increasing ticks fontsize
%xlim([X(1),X(end)]) 
%ylim([-0.15,1.05]) 
xlabel('$x$','Interpreter','latex') 
ylabel('density $\rho$','Interpreter','latex')
id = legend('ref','usual RBF','weak RBF','Interpreter','latex','Location','best');
set(id, 'Interpreter','latex', 'FontSize',22)

% plot velocity
figure(2) 
p = plot(XX,vvelocity_ref,'k:', 'LineWidth',2.5); 
hold on 
sz = 80;
scatter(X,velocity_s,sz,'rs', 'filled'); 
scatter(X,velocity_P1,sz,'bo', 'filled');
hold off
set(gca, 'FontSize', 24)  % Increasing ticks fontsize
%xlim([x(1),x(end)]) 
%ylim([-0.15,1.05]) 
xlabel('$x$','Interpreter','latex') 
ylabel('velocity $u$','Interpreter','latex')
id = legend('ref','usual RBF','weak RBF','Interpreter','latex','Location','best');
set(id, 'Interpreter','latex', 'FontSize',22)

% plot pressure
figure(3) 
p = plot(XX,ppressure_ref,'k:', 'LineWidth',2.5); 
hold on 
sz = 80;
scatter(X,pressure_s,sz,'rs', 'filled'); 
scatter(X,pressure_P1,sz,'bo', 'filled');
hold off
set(gca, 'FontSize', 24)  % Increasing ticks fontsize
%xlim([x(1),x(end)]) 
%ylim([-0.15,1.05]) 
xlabel('$x$','Interpreter','latex') 
ylabel('pressure $p$','Interpreter','latex')
id = legend('ref','usual RBF','weak RBF','Interpreter','latex','Location','best');
set(id, 'Interpreter','latex', 'FontSize',22)