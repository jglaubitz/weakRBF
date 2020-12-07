% We use the domain [-1,1] 
clear, clc, close all 

%% Setting up common variables 
Init_C = 'cos^2'; % sin, exp, cos^2
BC = 'inflow'; % inflow, periodic
T = 2; % final time 
basis = 'cubic'; % G, MQ, IQ, cubic, quintic
ep = 1; % shape parameter
N = 20; % number of points 
d = -1; % polynomial degree 
points = 'equid'; % equid, random
CFL = 0.1; % CFL number 
integration = 'exact'; % way integration is performed (exact, trapez, Gauss)

NN = []; ddx = [];
max_error_strong = []; max_error_weak_d0 = []; max_error_weak_d1 = []; 
L2_error_strong = []; L2_error_weak_d0 = []; L2_error_weak_d1 = [];
for N=20:20:100

    clear x u_strong u_weak_d0 u_weak_d1
    N
    
    %% Generating the collocation points 
    x = linspace(-1,1,N)'; % equidistant collocation points 
    if strcmp(points,'random')
        x(2:end-1) = 2*rand(1,N-2)-1; % random collocation points 
    end
    x = sort(x,'ascend'); 
    dx = max( abs( x(2:end) - x(1:end-1) ) ); % spatial step length

    %% set up RBF and IC 
    rbf = basis_function( basis );
    IC = initial_cond( Init_C ); 
    u0 = IC(x);

    %% routine for strong and weak RBF method 
    [u_strong, m_strong, e_strong] = linear_strong_RBF( BC, T, CFL, x, IC, rbf, ep ); % strong RBF
    [u_weak_d0, m_weak_d0, e_weak_d0] = linear_weak_RBF( BC, T, CFL, x, IC, rbf, ep, -1, integration ); % weak RBF without polynomials 
    [u_weak_d1, m_weak_d1, e_weak_d1] = linear_weak_RBF( BC, T, CFL, x, IC, rbf, ep, 0, integration ); % weak RBF with constant 

    %% Reference solution 
    u_ref = IC( mod(abs(x-T+1),2) - 1 ); % reference solution 
    %% Maximum error 
    NN = [NN;N]; % number of grid points 
    ddx = [ddx;dx]; % spatial resolution 
    error = max( abs( u_ref - u_strong ) ); % max error 
    max_error_strong = [max_error_strong;error]; % store max error for strong RBF method 
    error = max( abs( u_ref - u_weak_d0 ) ); % max error 
    max_error_weak_d0 = [max_error_weak_d0;error]; % store max error for weak RBF method 
    error = max( abs( u_ref - u_weak_d1 ) ); % max error 
    max_error_weak_d1 = [max_error_weak_d1;error]; % store max error for weak RBF method 
    %% Mean square error 
    error = norm( u_ref - u_strong )/sqrt(N); % mean square error 
    L2_error_strong = [L2_error_strong;error]; % store mean square error for strong RBF method 
    error = norm( u_ref - u_weak_d0 )/sqrt(N); % mean square error 
    L2_error_weak_d0 = [L2_error_weak_d0;error]; % store mean square error for weak RBF method 
    error = norm( u_ref - u_weak_d1 )/sqrt(N); % mean square error 
    L2_error_weak_d1 = [L2_error_weak_d1;error]; % store mean square error for weak RBF method 

end

%% Fit line 
c = polyfit(NN,max_error_weak_d0,1); % fit parameters c 
line_max_w = polyval(c,NN); % evaluate corresponding polynomial (linear)

%% plot maximum errors
figure(1) 
hold on 
sz = 80; 
if strcmp(points,'equid')
    scatter(NN,max_error_strong,sz,'rs', 'filled'); 
    scatter(NN,max_error_weak_d0,sz,'go', 'filled');
    scatter(NN,max_error_weak_d1,sz,'b^', 'filled');
else 
    scatter(ddx,max_error_strong,sz,'rs', 'filled'); 
    scatter(ddx,max_error_weak_d0,sz,'go', 'filled');
    scatter(ddx,max_error_weak_d1,sz,'b^', 'filled');
end
hold off
set(gca, 'FontSize', 24)  % Increasing ticks fontsize 
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
if strcmp(points,'equid')
    xlabel('$N$','Interpreter','latex') 
else 
    xlabel('$h$','Interpreter','latex')
end
ylabel('$\|u-u_N\|_\infty$','Interpreter','latex')
id = legend('usual RBF','weak RBF ($P=0$)','weak RBF ($P=1$)','Interpreter','latex','Location','best');
set(id, 'Interpreter','latex', 'FontSize',24)

%% plot L2 errors
figure(2) 
hold on 
sz = 80;
if strcmp(points,'equid')
    scatter(NN,L2_error_strong,sz,'rs', 'filled'); 
    scatter(NN,L2_error_weak_d0,sz,'go', 'filled');
    scatter(NN,L2_error_weak_d1,sz,'b^', 'filled');
else 
    scatter(ddx,L2_error_strong,sz,'rs', 'filled'); 
    scatter(ddx,L2_error_weak_d0,sz,'go', 'filled');
    scatter(ddx,L2_error_weak_d1,sz,'b^', 'filled');
end
hold off
set(gca, 'FontSize', 24)  % Increasing ticks fontsize
%xlim([NN(1),NN(end)]) 
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
if strcmp(points,'equid')
    xlabel('$N$','Interpreter','latex') 
else 
    xlabel('$h$','Interpreter','latex')
end
ylabel('$\|u-u_N\|_2$','Interpreter','latex')
id = legend('usual RBF','weak RBF ($P=0$)','weak RBF ($P=1$)','Interpreter','latex','Location','best');
set(id, 'Interpreter','latex', 'FontSize',24)