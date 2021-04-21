%% Exercise 2

%-------------------------------------------------------------------------%
% Maximum Lyapunov Exponent
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
dt = 0.001; % Time steps (dt)
t_final = 50; % time units

for j=1:length(dt)
% N_steps(j) = ceil(t_final/h(j)); % Number of steps, rounds up with ceil
N_steps(j) = length(0:dt:t_final)-1;
end

% Initial conditions
x0 = [0.1 0.105];
y0 = [0.1 0.105];
z0 = [0.1 0.105];

% For each initial condition
for k=1:length(x0)
        t(1) = 0; % We begin at t=0 s
        x(k,1) = x0(k); % x(0)
        y(k,1) = y0(k); % y(0)
        z(k,1) = z0(k); % y(0)
            % Euler method's update loop
            for i=1:N_steps(j)
                t(i+1) = t(i) + dt;
                x(k,i+1) = x(k,i) + f1(t(i),x(k,i),y(k,i),z(k,i))*dt;
                y(k,i+1) = y(k,i) + f2(t(i),x(k,i),y(k,i),z(k,i))*dt;
                z(k,i+1) = z(k,i) + f3(t(i),x(k,i),y(k,i),z(k,i))*dt;
            end       
end

%% Plots

% X
plot_pdf = figure(1);
delta = log(abs(x(2,:)-x(1,:)));
plot(t,delta);
hold on;

init=2000;
final=length(t)/3-1000;


p = polyfit(t(init:final),delta(init:final),1);
f1 = polyval(p,t(init:final));
plot(t(init:final),f1,'LineWidth',2);

xlabel('Time units')
ylabel('$\ln{\left||\delta\right||}$')
title('\textbf{Maximum Lyapunov Exponent for $\rho = 30$}')
str = {strcat('$\lambda = $' , num2str(p(1)))};
str = [str , strcat('$\lambda = $' , num2str(p(1)))];
legend('x',str{2},'location','best');
box on
grid minor
hold off;


% Save pdf
set(plot_pdf, 'Units', 'Centimeters');
pos = get(plot_pdf, 'Position');
set(plot_pdf, 'PaperPositionMode', 'Auto', 'PaperUnits', 'Centimeters', ...
    'PaperSize',[pos(3), pos(4)]);
print(plot_pdf, 'MLE_rho_30_x.pdf', '-dpdf', '-r0');

% Save png
print(gcf,'MLE_rho_30_x.png','-dpng','-r600');



% Y
plot_pdf2 = figure(2);
delta = log(abs(y(2,:)-y(1,:)));
plot(t,delta);
hold on;

p = polyfit(t(init:final),delta(init:final),1);
f1 = polyval(p,t(init:final));
plot(t(init:final),f1,'LineWidth',2);

xlabel('Time units')
ylabel('$\ln{\left||\delta\right||}$')
title('\textbf{Maximum Lyapunov Exponent for $\rho = 30$}')
str = {strcat('$\lambda = $' , num2str(p(1)))};
str = [str , strcat('$\lambda = $' , num2str(p(1)))];
legend('y',str{2},'location','best');
box on
grid minor
hold off;

% Save pdf
set(plot_pdf2, 'Units', 'Centimeters');
pos = get(plot_pdf2, 'Position');
set(plot_pdf2, 'PaperPositionMode', 'Auto', 'PaperUnits', 'Centimeters', ...
    'PaperSize',[pos(3), pos(4)]);
print(plot_pdf2, 'MLE_rho_30_y.pdf', '-dpdf', '-r0');

% Save png
print(gcf,'MLE_rho_30_y.png','-dpng','-r600');


% Z
plot_pdf3 = figure(3);
delta = log(abs(z(2,:)-z(1,:)));
plot(t,delta);
hold on;

p = polyfit(t(init:final),delta(init:final),1);
f1 = polyval(p,t(init:final));
plot(t(init:final),f1,'LineWidth',2);

xlabel('Time units')
ylabel('$\ln{\left||\delta\right||}$')
title('\textbf{Maximum Lyapunov Exponent for $\rho = 30$}')
str = {strcat('$\lambda = $' , num2str(p(1)))};
str = [str , strcat('$\lambda = $' , num2str(p(1)))];
legend('z',str{2},'location','best');
box on
grid minor
hold off;


% Save pdf
set(plot_pdf3, 'Units', 'Centimeters');
pos = get(plot_pdf3, 'Position');
set(plot_pdf3, 'PaperPositionMode', 'Auto', 'PaperUnits', 'Centimeters', ...
    'PaperSize',[pos(3), pos(4)]);
print(plot_pdf3, 'MLE_rho_30_z.pdf', '-dpdf', '-r0');

% Save png
print(gcf,'MLE_rho_30_z.png','-dpng','-r600');


