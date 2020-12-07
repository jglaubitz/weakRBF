%% Script for the error plots: Periodic BCs + cubic

clear, clc, close all 
points = 'equid'; 

%% Periodic, cubic - max error 

NN = [20  40  60  80 100]'; 
error_strong = [1.3697     0.51185     0.22937     0.10244    0.059246]'; 
error_weak_d0 = [0.80844    0.058769    0.022821     0.01256   0.0088146]'; 
error_weak_d1 = [0.79529    0.055924     0.02206    0.012297    0.008678]';

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
text(25,10^(-0.4),'C(N-20)^{-2.5}','color','k','fontsize',20)
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


%% Periodic, cubic - L2 error 

error_strong = [0.58883     0.29125     0.11829    0.054275    0.030387]'; 
error_weak_d0 = [0.44627    0.030629    0.012724   0.0076686   0.0054573]'; 
error_weak_d1 = [0.44666    0.029167    0.012344    0.007524   0.0053875]';

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
text(25,10^(-0.6),'C(N-20)^{-2.5}','color','k','fontsize',20)
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
id = legend('usual RBF','weak RBF ($P=0$)','weak RBF ($P=1$)','Interpreter','latex','Location','southwest');
set(id, 'Interpreter','latex', 'FontSize',24)