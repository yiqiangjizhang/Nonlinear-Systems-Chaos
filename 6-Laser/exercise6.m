%% Laser turn on simulation

%-------------------------------------------------------------------------%
% Exercise 6_1: Simulate the 'turn on' when r is constant, r>r star=0.
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

%% Laser turn on

% Constants and Parameters
r = 1:1:3;

% Function handle
f1 = @(t,x,r) r*x - x^2;

% 1.1 Numerical data
dt = [0.001]; % Time steps (dt)
t_final = 20; % time units

V_steps = zeros(1,length(dt));

for j=1:length(dt)
% N_steps(j) = ceil(t_final/h(j)); % Number of steps, rounds up with ceil
N_steps(j) = length(0:dt:t_final)-1;
end

% Initial conditions
x0 = 0.01;

% For each h (dt)
for j=1:length(dt)
    
    t(j,1) = 0; % We begin at t=0 s
    for k = 1:length(x0)
    x(j,1) = x0(k); % x(0) = 1
    
    for k=1:length(r)
        % Euler method's update loop
        for i=1:N_steps(j)
            t(j,i+1) = t(j,i) + dt(j);
            x(j,i+1) = x(j,i) + f1(t(j,i),x(j,i),r(k))*dt(j);
        end
        plot_pdf = figure(1);
        plot(t(j,:),x(j,:));
        hold on;
    end
    end
end
xlim([0 20])
xlabel('Time units')
ylabel('x')
title('\textbf{Laser turn on}')
legend('$r=1$', '$r=2$', '$r=3$','location','east')
box on
grid minor
hold off;

% Save pdf
set(plot_pdf, 'Units', 'Centimeters');
pos = get(plot_pdf, 'Position');
set(plot_pdf, 'PaperPositionMode', 'Auto', 'PaperUnits', 'Centimeters', ...
    'PaperSize',[pos(3), pos(4)]);
print(plot_pdf, 'laser_turn_on_v.pdf', '-dpdf', '-r0');

% Save png
print(gcf,'laser_turn_on_v.png','-dpng','-r600');
