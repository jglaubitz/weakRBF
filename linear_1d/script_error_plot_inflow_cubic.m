%% Script for the error plots: Inflow BCs + cubic

clear, clc, close all 
points = 'equid'; 

%% Inflow, cubic - max error 

NN = [20  40  60  80 100]'; 
error_strong = [1.3542     0.43415     0.14914    0.057919    0.036106]'; 
error_weak_d0 = [0.69547    0.078926    0.041683    0.027366     0.02106]'; 
error_weak_d1 = [0.71225    0.076857    0.041017    0.027177    0.020979]';

x1 = NN(end);
x2 = NN(1); 
y1_s = error_strong(end);
y2_s = (x1/x2)^1 * y1_s;
y1_w = error_weak_d0(end);
y2_w = (x1/x2)^(2) * y1_w;  

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
text(25,10^(-0.35),'C(N-20)^{-2}','color','k','fontsize',20)
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
id = legend('usual RBF','weak RBF ($P=0$)','weak RBF ($P=1$)','Interpreter','latex','Location','southwest');
set(id, 'Interpreter','latex', 'FontSize',24)


%% Inflow, cubic - L2 error 

error_strong = [0.56687     0.20914    0.070211    0.029245    0.016112]'; 
error_weak_d0 = [0.37719    0.039827    0.020353    0.013836    0.010594]'; 
error_weak_d1 = [0.37935    0.038601    0.020077    0.013746    0.010555]';

x1 = NN(end);
x2 = NN(1); 
y1_s = error_strong(end);
y2_s = (x1/x2)^1 * y1_s;
y1_w = error_weak_d0(end);
y2_w = (x1/x2)^(2) * y1_w; 

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
text(25,10^(-0.65),'C(N-20)^{-2}','color','k','fontsize',20)
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