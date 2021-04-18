%% Exercise 1

%-------------------------------------------------------------------------%
% Adler Equation
%-------------------------------------------------------------------------%

% Date: 09/03/2021
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

% Function handle
a = 1;
b = 2;
c = 3;
d = 1;


f1 = @(t,x,y) a*x + b*y;
f2 = @(t,x,y) c*x + d*y;


% 1.1 Numerical data
dt = [0.01]; % Time steps (dt)
t_final = 10; % time units

for j=1:length(dt)
% N_steps(j) = ceil(t_final/h(j)); % Number of steps, rounds up with ceil
N_steps(j) = length(0:dt:t_final)-1;
end

% Initial conditions
x0 = 0.2;
y0 = 1.5;

% For each h (dt)
for j=1:length(dt)
    t(j,1) = 0; % We begin at t=0 s
    x(j,1) = x0; % x(0)
    y(j,1) = y0; % y(0)
        % Euler method's update loop
        for i=1:N_steps(j)
            t(j,i+1) = t(j,i) + dt(j);
            x(j,i+1) = x(j,i) + f1(t(j,i),x(j,i),y(j,i))*dt(j);
            y(j,i+1) = y(j,i) + f2(t(j,i),x(j,i),y(j,i))*dt(j);
        end
        plot(t(j,:),x(j,:),t(j,:),y(j,:));
        hold on;
end

plot_pdf = figure(1);
xlabel('Time units')
ylabel('Value')
title('Two-dimensional linear systems')
legend('x','y','location','southwest')
box on
grid minor
hold off;

% Save pdf
set(plot_pdf, 'Units', 'Centimeters');
pos = get(plot_pdf, 'Position');
set(plot_pdf, 'PaperPositionMode', 'Auto', 'PaperUnits', 'Centimeters', ...
    'PaperSize',[pos(3), pos(4)]);
print(plot_pdf, 'exercise1_pre_2DLinearSystem.pdf', '-dpdf', '-r0');

% Save png
print(gcf,'exercise1_pre_2DLinearSystem.png','-dpng','-r600');


