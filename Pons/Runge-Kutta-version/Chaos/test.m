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
rho = [21 24.15 30];


% Function handle
f1 = @(t,x,y,z) sigma*(y - x);
f2 = @(t,x,y,z) rho(3)*x - y -x*z;
f3 = @(t,x,y,z) -beta*z + x*y;


% 1.1 Numerical data
dt = [0.001]; % Time steps (dt)
t_final = 50; % time units

for j=1:length(dt)
% N_steps(j) = ceil(t_final/h(j)); % Number of steps, rounds up with ceil
N_steps(j) = length(0:dt:t_final)-1;
end

% Initial conditions
x0 = 0.1:0.1:0.3;
y0 = 0.1:0.1:0.3;
z0 = 0.1:0.1:0.3;

% For each initial condition
for k=1:length(x0)
    % For each h (dt)
    for j=1:length(dt)
        t(j,1) = 0; % We begin at t=0 s
        x(j,1) = x0(k); % x(0)
        y(j,1) = y0(k); % y(0)
        z(j,1) = z0(k); % y(0)
            % Euler method's update loop
            for i=1:N_steps(j)
                t(j,i+1) = t(j,i) + dt(j);
                x(j,i+1) = x(j,i) + f1(t(j,i),x(j,i),y(j,i),z(j,i))*dt(j);
                y(j,i+1) = y(j,i) + f2(t(j,i),x(j,i),y(j,i),z(j,i))*dt(j);
                z(j,i+1) = z(j,i) + f3(t(j,i),x(j,i),y(j,i),z(j,i))*dt(j);
            end
            figure(1)
            plot(t(j,:), x(j,:));
            hold on;
            figure(2)
            plot(t(j,:), y(j,:));
            hold on;
            figure(3)
            plot(t(j,:), z(j,:));
            hold on;
            figure(4)
            plot(x(j,:), z(j,:));
            hold on;
    end
end

plot_pdf = figure(1);
xlabel('t')
ylabel('x')
title('\textbf{Lorenz System for $\rho = 30$}')
legend('$x_0, y_0, z_0 = 0.1$','$x_0, y_0, z_0 = 0.2$','$x_0, y_0, z_0 = 0.3$',...
'location','best')
box on
grid minor
hold off;


plot_pdf2 = figure(2);
xlabel('t')
ylabel('y')
title('\textbf{Lorenz System for $\rho = 30$}')
legend('$x_0, y_0, z_0 = 0.1$','$x_0, y_0, z_0 = 0.2$','$x_0, y_0, z_0 = 0.3$',...
'location','best')
box on
grid minor
hold off;


plot_pdf3 = figure(3);
xlabel('t')
ylabel('z')
title('\textbf{Lorenz System for $\rho = 30$}')
legend('$x_0, y_0, z_0 = 0.1$','$x_0, y_0, z_0 = 0.2$','$x_0, y_0, z_0 = 0.3$',...
'location','best')
box on
grid minor
hold off;

plot_pdf4 = figure(4);
xlabel('x')
ylabel('z')
title('\textbf{Lorenz System Phase Portrait for $\rho = 30$}')
box on
grid minor
hold off;







