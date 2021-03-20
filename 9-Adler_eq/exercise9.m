%% Adler Equation

%-------------------------------------------------------------------------%
% Exercise 9 : Solve Adler’s equation with w=1, a=0.99 and theta(0)=pi/2
%-------------------------------------------------------------------------%

% Date: 12/03/2021
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

%% Adler equation

% Constants and Parameters
w = 1;
a = 0.99;

% Function handle
f1 = @(t,x,a,w) w-a*sin(x);

% 1.1 Numerical data
dt = [0.01]; % Time steps (dt)
t_final = 150; % time units

for j=1:length(dt)
% N_steps(j) = ceil(t_final/h(j)); % Number of steps, rounds up with ceil
N_steps(j) = length(0:dt:t_final)-1;
end

% Initial conditions
x0 = pi/2;

% For each h (dt)
for j=1:length(dt)
    t(j,1) = 0; % We begin at t=0 s
    for k = 1:length(x0)
    x(j,1) = x0(k); % x(0) = 1

        % Euler method's update loop
        for i=1:N_steps(j)
            t(j,i+1) = t(j,i) + dt(j);
            x(j,i+1) = x(j,i) + f1(t(j,i),x(j,i),a,w)*dt(j);
        end
        plot_pdf = figure(1);
        plot(t(j,:),x(j,:));
        hold on;
    end
end
xlim([0 150])
xlabel('Time units')
ylabel('$\theta$')
title('\textbf{Adler equation}')
box on
grid minor
hold off;

% Save pdf
set(plot_pdf, 'Units', 'Centimeters');
pos = get(plot_pdf, 'Position');
set(plot_pdf, 'PaperPositionMode', 'Auto', 'PaperUnits', 'Centimeters', ...
    'PaperSize',[pos(3), pos(4)]);
print(plot_pdf, 'adler_eq_t_theta.pdf', '-dpdf', '-r0');

% Save png
print(gcf,'adler_eq_theta_t_theta.png','-dpng','-r600');
