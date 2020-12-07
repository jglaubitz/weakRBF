%% Script for the error plots: Inflow BCs + quintic

clear, clc, close all 
points = 'equid'; 

%% Inflow, quintic - max error 

NN = [20  40  60  80 100]'; 
error_strong = [2.0275     0.31487    0.069455    0.076577    0.077225]'; 
error_weak_d0 = [0.77098    0.092093    0.036358    0.023928    0.018809]'; 
error_weak_d1 = [0.71554    0.099924    0.036512    0.023961    0.018833]';

x1 = NN(end);
x2 = NN(1); 
y1_s = error_strong(end);
y2_s = (x1/x2)^1 * y1_s;
y1_w = error_weak_d0(end);
y2_w = (x1/x2)^2 * y1_w; 

figure(1) 
hold on 
sz = 80; 
if strcmp(points,'equid')
    scatter(NN,error_strong,sz,'rs', 'filled'); 
    scatter(NN,error_weak_d0,sz,'go', 'filled');
    scatter(NN,error_weak_d1,sz,'b^', 'filled');
else 
    scatter(ddx,error_strong,sz,'rs', 'filled'); 
    scatter(ddx,error_weak_d0,sz,'go', 'filled');
    scatter(ddx,error_weak_d1,sz,'b^', 'filled');
end 
% 2nd order slope
plot([x1 x2],[y1_w y2_w],':k', 'LineWidth',2.5); 
text(22,10^(-0.8),'C(N-20)^{-2}','color','k','fontsize',20)
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
id = legend('usual RBF','weak RBF ($P=0$)','weak RBF ($P=1$)','Interpreter','latex','Location','northeast');
set(id, 'Interpreter','latex', 'FontSize',24)


%% Inflow, quintic - L2 error 

error_strong = [0.791     0.17436    0.039887    0.027733    0.027403]'; 
error_weak_d0 = [0.39869    0.037988    0.018682    0.012924    0.010013]'; 
error_weak_d1 = [0.39046    0.038858    0.018797    0.012951    0.010024]';

x1 = NN(end);
x2 = NN(1); 
y1_s = error_strong(end);
y2_s = (x1/x2)^1 * y1_s;
y1_w = error_weak_d0(end);
y2_w = (x1/x2)^2 * y1_w;  

figure(2) 
hold on 
sz = 80; 
if strcmp(points,'equid')
    scatter(NN,error_strong,sz,'rs', 'filled'); 
    scatter(NN,error_weak_d0,sz,'go', 'filled');
    scatter(NN,error_weak_d1,sz,'b^', 'filled');
else 
    scatter(ddx,error_strong,sz,'rs', 'filled'); 
    scatter(ddx,error_weak_d0,sz,'go', 'filled');
    scatter(ddx,error_weak_d1,sz,'b^', 'filled');
end 
% 2nd order slope
plot([x1 x2],[y1_w y2_w],':k', 'LineWidth',2.5); 
text(22,10^(-1.1),'C(N-20)^{-2}','color','k','fontsize',20)
hold off
set(gca, 'FontSize', 24)  % Increasing ticks fontsize 
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
if strcmp(points,'equid')
    xlabel('$N$','Interpreter','latex') 
else 
    xlabel('$h$','Interpreter','latex')
end
ylabel('$\|u-u_N\|_2$','Interpreter','latex')
id = legend('usual RBF','weak RBF ($P=0$)','weak RBF ($P=1$)','Interpreter','latex','Location','northeast');
set(id, 'Interpreter','latex', 'FontSize',24)