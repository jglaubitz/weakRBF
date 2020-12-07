% Script to compare the strong and weak RBF method for the linear advection
% equation 

% We use the domain [-1,1] 
clear, clc, close all 

%% Setting up common variables 
Init_C = 'exp'; % sin, exp, disc
BC = 'periodic'; % inflow, periodic
T = 2; % final time 
basis = 'quintic'; % G, MQ, IQ, cubic, quintic
ep = 5; % shape parameter
N = 20; % number of points 
d = -1; % polynomial degree 
points = 'equid'; % equid, random
CFL = 0.1; % CFL number 
integration = 'exact'; % way integration is performed (exact, trapez, Gauss)

%% Generating the collocation points 
x = linspace(-1,1,N)'; % equidistant collocation points 
if strcmp(points,'random')
    x(2:end-1) = 2*rand(1,N-2)-1; % random collocation points 
end
x = sort(x,'ascend');

%% set up RBF and IC 
rbf = basis_function( basis );
IC = initial_cond( Init_C ); 
u0 = IC(x);

%% routine for strong and weak RBF method 
[u_strong, m_strong, e_strong] = linear_strong_RBF( BC, T, CFL, x, IC, rbf, ep ); % strong RBF
[u_weak_d0, m_weak_d0, e_weak_d0] = linear_weak_RBF( BC, T, CFL, x, IC, rbf, ep, -1, integration ); % weak RBF without polynomials 
[u_weak_d1, m_weak_d1, e_weak_d1] = linear_weak_RBF( BC, T, CFL, x, IC, rbf, ep, 0, integration ); % weak RBF with constant 

%% Reference solution 
u_ref = u0; 
xx = linspace(-1,1,100*N)';
uu_ref = IC(xx);

% plot numerical solution 
figure(1) 
p = plot(xx,uu_ref,'k:', 'LineWidth',2.5); 
hold on 
sz = 80;
scatter(x,u_strong,sz,'rs', 'filled'); 
scatter(x,u_weak_d0,sz,'go', 'filled');
scatter(x,u_weak_d1,sz,'b^', 'filled');
hold off
set(gca, 'FontSize', 24)  % Increasing ticks fontsize
xlim([x(1),x(end)]) 
ylim([-0.15,1.05]) 
xlabel('$x$','Interpreter','latex') 
ylabel('$u$','Interpreter','latex')
id = legend('ref','usual','weak ($P=0$)','weak ($P=1$)','Interpreter','latex','Location','best');
set(id, 'Interpreter','latex', 'FontSize',22)

% plot momentum 
figure(2)
p = plot(m_strong(:,1),m_strong(:,2),'r--', m_weak_d0(:,1),m_weak_d0(:,2),'g-.', m_weak_d1(:,1),m_weak_d1(:,2),'b-'); 
set(p, 'LineWidth',2.5)
set(gca, 'FontSize', 24)  % Increasing ticks fontsize
xlabel('$t$','Interpreter','latex') 
ylabel('$\int u \, {\rm{d}}x$','Interpreter','latex')
id = legend('usual RBF','weak RBF ($P=0$)','weak RBF ($P=1$)','Interpreter','latex','Location','best');
set(id, 'Interpreter','latex', 'FontSize',22)

% plot energy 
figure(3)
p = plot(e_strong(:,1),e_strong(:,2),'r--', e_weak_d0(:,1),e_weak_d0(:,2),'g-.', e_weak_d1(:,1),e_weak_d1(:,2),'b-'); 
set(p, 'LineWidth',2.5)
set(gca, 'FontSize', 24)  % Increasing ticks fontsize
xlabel('$t$','Interpreter','latex') 
ylabel('$\|u\|_2^2$','Interpreter','latex')
id = legend('usual RBF','weak RBF ($P=0$)','weak RBF ($P=1$)','Interpreter','latex','Location','best');
set(id, 'Interpreter','latex', 'FontSize',22)