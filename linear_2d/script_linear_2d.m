% Script to compare the strong and weak RBF method for the linear advection
% equation in two dimensions 

% We use the domain [-1,1]^2 
clear, clc, close all 

%% Setting up common variables 
Init_C = 'sin'; % sin, exp, disc
BC = 'periodic'; % inflow, periodic
T = 2; % final time 
kernel = 'quintic'; % G, MQ, IQ, cubic, quintic
ep = 1; % shape parameter
N = 20; % number of points 
d = 0; % polynomial degree 
points = 'equid'; % equid, random
CFL = 0.1; % CFL number 
integration = 'trapez'; % way integration is performed (exact, trapez, Gauss, LS)

%% Generating the collocation points 
[xx, yy, X] = grid_points_2d(-1,1,N,points); % generate grid points 
if strcmp(points,'random')
    if N^2 == 400
        load = matfile(['matrices/X_N=',num2str(N^2),'_',points,'.mat']);
        X = load.X;
    else 
        save( ['matrices/X_N=',num2str(N^2),'_',points,'.mat'], 'X' );
    end
end 

%% set up RBF, IC, and reference solution  
rbf = basis_function( kernel );
[IC, ref] = initial_cond_2d( Init_C, BC ); 
u0 = IC(X(:,1),X(:,2));
u_ref = ref(T,X(:,1),X(:,2));

%% routine for strong and weak RBF method 
u_strong = linear_strong_RBF_2d( BC, T, CFL, X, u0, rbf, ep ); % strong RBF
u_weak = linear_weak_RBF_2d( BC, T, CFL, X, u0, kernel, rbf, ep, points, d, integration ); % weak RBF with constant 

%% Rearrange solutions 
if strcmp(points,'equid') 
    uu_ref = reshape( u_ref, N, N );; % for plotting
    uu_strong = reshape( u_strong, N, N ); 
    uu_weak = reshape( u_weak, N, N );
elseif strcmp(points,'random') 
    f = scatteredInterpolant(X(:,1),X(:,2),u_ref); 
    uu_ref = f(xx,yy); 
    f = scatteredInterpolant(X(:,1),X(:,2),u_strong); 
    uu_strong = f(xx,yy); 
    f = scatteredInterpolant(X(:,1),X(:,2),u_weak); 
    uu_weak = f(xx,yy); 
else
    error('Desried points not implemented yet!')
end

% plot refrence solution
figure(1)
s = mesh(xx,yy,uu_ref); 
s.EdgeColor = 'interp'; 
set(s, 'LineWidth',1.5); 
set(gca, 'FontSize', 20)  % Increasing ticks fontsize 
xlabel('$x$','Interpreter','latex') 
ylabel('$y$','Interpreter','latex') 
zlabel('$u$','Interpreter','latex') 
zlim([-2,2])
%str = sprintf( ['plots/testFun1_2d.fig'] );
%savefig(str); 

% plot strong RBF
figure(2)
s = mesh(xx,yy,uu_strong); 
s.EdgeColor = 'interp'; 
set(s, 'LineWidth',1.5); 
set(gca, 'FontSize', 20)  % Increasing ticks fontsize 
xlabel('$x$','Interpreter','latex') 
ylabel('$y$','Interpreter','latex') 
zlabel('$u$','Interpreter','latex') 
zlim([-2,2])
%str = sprintf( ['figures/linear_2d_usual_',kernel,'_',points,'_T=',num2str(T),'.fig'] );
%savefig(str); 

% plot weak RBF
figure(3)
s = mesh(xx,yy,uu_weak); 
s.EdgeColor = 'interp'; 
set(s, 'LineWidth',1.5); 
set(gca, 'FontSize', 20)  % Increasing ticks fontsize 
xlabel('$x$','Interpreter','latex') 
ylabel('$y$','Interpreter','latex') 
zlabel('$u$','Interpreter','latex') 
zlim([-2,2])
%str = sprintf( ['figures/linear_2d_weak_',kernel,'_',points,'_T=',num2str(T),'.fig'] );
%savefig(str);