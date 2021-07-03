%% Exercise 1

%-------------------------------------------------------------------------%
% Chaotic systems
%-------------------------------------------------------------------------%

% Date: 10/04/2021
% Author/s: Yi Qiang Ji Zhang
% Subject: Nonlinear Systems, Chaos and Control in Engineering
% Professor: Antonio Pons & Cristina Masoller

% Clear workspace, command window and close windows
clear;
close all;
clc;
% Set interpreter to latex
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

% Constants
sigma = 10;
beta = 8/3;
rho_vec = [21 24.15 30];
rho = rho_vec(3);

% Numerical data
h = 0.001; % Time steps of RK4
t_final = 50; % time units
N_steps = ceil(t_final/h); % Number of steps, rounds up with ceil

% Declaration of solution vectors
t = 0:h:t_final;
x = zeros(1, length(t));
y = zeros(1, length(t));
z = zeros(1, length(t));

% Function handle
f1 = @(t,x,y,z) sigma*(y - x);
f2 = @(t,x,y,z) rho*x - y -x*z;
f3 = @(t,x,y,z) -beta*z + x*y;


% Initial conditions
x0 = 0.1:0.1:0.3;
y0 = 0.1:0.1:0.3;
z0 = 0.1:0.1:0.3;


% For each x0
for j=1:length(x0)
    % Initial conditions
    t(1) = 0; % We begin at t=0 s
    x(j,1) = x0(j);
    y(j,1) = y0(j);
    z(j,1) = z0(j);
    
    % RK4 update loop
    for i=1:N_steps

        % Update function f1 and f2
        % Compute coefficients sub 1
        k1_f1 = f1(t(i)         ,x(j,i)                   ,y(j,i)               ,z(j,i));
        k1_f2 = f2(t(i)         ,x(j,i)                   ,y(j,i)               ,z(j,i));
        k1_f3 = f3(t(i)         ,x(j,i)                   ,y(j,i)               ,z(j,i));

        % Compute coefficients sub 2
        k2_f1 = f1(t(i) + h/2   ,x(j,i) + k1_f1*h/2       ,y(j,i) + k1_f2*h/2   ,z(j,i) + k1_f3*h/2);
        k2_f2 = f2(t(i) + h/2   ,x(j,i) + k1_f1*h/2       ,y(j,i) + k1_f2*h/2   ,z(j,i) + k1_f3*h/2);
        k2_f3 = f3(t(i) + h/2   ,x(j,i) + k1_f1*h/2       ,y(j,i) + k1_f2*h/2   ,z(j,i) + k1_f3*h/2);

        % Compute coefficients sub 3
        k3_f1 = f1(t(i) + h/2   ,x(j,i) + k2_f1*h/2       ,y(j,i) + k2_f2*h/2   ,z(j,i) + k2_f3*h/2);
        k3_f2 = f2(t(i) + h/2   ,x(j,i) + k2_f1*h/2       ,y(j,i) + k2_f2*h/2   ,z(j,i) + k2_f3*h/2);
        k3_f3 = f3(t(i) + h/2   ,x(j,i) + k2_f1*h/2       ,y(j,i) + k2_f2*h/2   ,z(j,i) + k2_f3*h/2);

        % Compute coefficients sub 4
        k4_f1 = f1(t(i) + h     ,x(j,i) + k3_f1*h         ,y(j,i) + k3_f2*h     ,z(j,i) + k3_f3*h);
        k4_f2 = f2(t(i) + h     ,x(j,i) + k3_f1*h         ,y(j,i) + k3_f2*h     ,z(j,i) + k3_f3*h);
        k4_f3 = f3(t(i) + h     ,x(j,i) + k3_f1*h         ,y(j,i) + k3_f2*h     ,z(j,i) + k3_f3*h);

        % Compute next step
        x(j,i+1) = x(j,i) + h/6*(k1_f1 + 2*k2_f1 + 2*k3_f1 + k4_f1);
        y(j,i+1) = y(j,i) + h/6*(k1_f2 + 2*k2_f2 + 2*k3_f2 + k4_f2);
        z(j,i+1) = z(j,i) + h/6*(k1_f3 + 2*k2_f3 + 2*k3_f3 + k4_f3);

    end
    figure(1)
    plot(t, x(j,:));
    hold on;
    figure(2)
    plot(t, y(j,:));
    hold on;
    figure(3)
    plot(t, z(j,:));
    hold on;
    figure(4)
    plot(x(j,:), z(j,:));
    hold on;
end


str0 = {strcat('$\rho = $' , num2str(rho))};
str0 = [str0 , strcat('$\rho = $' , num2str(rho))];

plot_pdf = figure(1);
xlabel('t')
ylabel('x')
title('\textbf{Lorenz System for}',str0{2})
legend('$x_0, y_0, z_0 = 0.1$','$x_0, y_0, z_0 = 0.2$','$x_0, y_0, z_0 = 0.3$',...
'location','best')
box on
grid minor
hold off;

% Save pdf
set(plot_pdf, 'Units', 'Centimeters');
pos = get(plot_pdf, 'Position');
set(plot_pdf, 'PaperPositionMode', 'Auto', 'PaperUnits', 'Centimeters', ...
    'PaperSize',[pos(3), pos(4)]);
print(plot_pdf, 'lorenz_system_rho_30_x.pdf', '-dpdf', '-r0');

% Save png
print(gcf,'lorenz_system_rho_30_x.png','-dpng','-r600');


plot_pdf2 = figure(2);
xlabel('t')
ylabel('y')
title('\textbf{Lorenz System for}',str0{2})
legend('$x_0, y_0, z_0 = 0.1$','$x_0, y_0, z_0 = 0.2$','$x_0, y_0, z_0 = 0.3$',...
'location','best')
box on
grid minor
hold off;

% Save pdf
set(plot_pdf2, 'Units', 'Centimeters');
pos = get(plot_pdf2, 'Position');
set(plot_pdf2, 'PaperPositionMode', 'Auto', 'PaperUnits', 'Centimeters', ...
    'PaperSize',[pos(3), pos(4)]);
print(plot_pdf2, 'lorenz_system_rho_30_y.pdf', '-dpdf', '-r0');

% Save png
print(gcf,'lorenz_system_rho_30_y.png','-dpng','-r600');

plot_pdf3 = figure(3);
xlabel('t')
ylabel('z')
title('\textbf{Lorenz System for}',str0{2})
legend('$x_0, y_0, z_0 = 0.1$','$x_0, y_0, z_0 = 0.2$','$x_0, y_0, z_0 = 0.3$',...
'location','best')
box on
grid minor
hold off;

% Save pdf
set(plot_pdf3, 'Units', 'Centimeters');
pos = get(plot_pdf3, 'Position');
set(plot_pdf3, 'PaperPositionMode', 'Auto', 'PaperUnits', 'Centimeters', ...
    'PaperSize',[pos(3), pos(4)]);
print(plot_pdf3, 'lorenz_system_rho_30_z.pdf', '-dpdf', '-r0');

% Save png
print(gcf,'lorenz_system_rho_30_z.png','-dpng','-r600');


plot_pdf4 = figure(4);
xlabel('x')
ylabel('z')
title('\textbf{Lorenz System Phase Portrait for}',str0{2})
box on
grid minor
hold off;

% Save pdf
set(plot_pdf4, 'Units', 'Centimeters');
pos = get(plot_pdf4, 'Position');
set(plot_pdf4, 'PaperPositionMode', 'Auto', 'PaperUnits', 'Centimeters', ...
    'PaperSize',[pos(3), pos(4)]);
print(plot_pdf4, 'lorenz_system_rho_30_xz.pdf', '-dpdf', '-r0');

% Save png
print(gcf,'lorenz_system_rho_30_xz.png','-dpng','-r600');




