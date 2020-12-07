%% Script for the error plots: Periodic BCs + quintic

clear, clc, close all 
points = 'equid'; 

%% Periodic, quintic - max error 

NN = [20  40  60  80 100]'; 
error_strong = [2.0572     0.34604    0.096451     0.10333     0.30632]'; 
error_weak_d0 = [1.0796    0.079408      0.0331    0.019945    0.014225]'; 
error_weak_d1 = [1.064    0.082802    0.034367    0.020464    0.014492]';

x1 = NN(end);
x2 = NN(1); 
y1_s = error_strong(end);
y2_s = (x1/x2)^1 * y1_s;
y1_w = error_weak_d0(end);
y2_w = (x1/x2)^(2.5) * y1_w;  

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
plot([x1 x2],[y1_w y2_w],':k', 'LineWidth',2.5); 
text(25,10^(-0.2),'C(N-20)^{-2.5}','color','k','fontsize',20)
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


%% Periodic, quintic - L2 error 

error_strong = [0.82294     0.19464    0.058369    0.061536     0.14773]'; 
error_weak_d0 = [0.48332    0.027895    0.010638   0.0063155    0.004465]'; 
error_weak_d1 = [0.48014    0.028952    0.010772   0.0063455   0.0044775]';

x1 = NN(end);
x2 = NN(1); 
y1_s = error_strong(end);
y2_s = (x1/x2)^1 * y1_s;
y1_w = error_weak_d0(end);
y2_w = (x1/x2)^(2.5) * y1_w; 

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
plot([x1 x2],[y1_w y2_w],':k', 'LineWidth',2.5); 
text(25,10^(-0.7),'C(N-20)^{-2.5}','color','k','fontsize',20)
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
id = legend('usual RBF','weak RBF ($P=0$)','weak RBF ($P=1$)','Interpreter','latex','Location','best');
set(id, 'Interpreter','latex', 'FontSize',24)