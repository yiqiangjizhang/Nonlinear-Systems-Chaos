%% Exercise 3

%-------------------------------------------------------------------------%
% Two-dimensional linear systems
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

% Function handle
f1 = @(t,x,y) -x + 4*x^3;
f2 = @(t,x,y) -2*x;


% 1.1 Numerical data
dt = [0.01]; % Time steps (dt)
t_final = 10; % time units

for j=1:length(dt)
% N_steps(j) = ceil(t_final/h(j)); % Number of steps, rounds up with ceil
N_steps(j) = length(0:dt:t_final)-1;
end

% Initial conditions
x0 = [-0.45 -0.40 -0.30 -0.25 -0.10 0.10 0.25 0.30 0.40 0.45];
y0 = [-0.45 -0.40 -0.30 -0.25 -0.10 0.10 0.25 0.30 0.40 0.45];

% For each initial condition
for k=1:length(x0)
    % For each h (dt)
    for j=1:length(dt)
        t(j,1) = 0; % We begin at t=0 s
        x(j,1) = x0(k); % x(0)
        y(j,1) = y0(k); % y(0)
            % Euler method's update loop
            for i=1:N_steps(j)
                t(j,i+1) = t(j,i) + dt(j);
                x(j,i+1) = x(j,i) + f1(t(j,i),x(j,i),y(j,i))*dt(j);
                y(j,i+1) = y(j,i) + f2(t(j,i),x(j,i),y(j,i))*dt(j);
            end
            plot(x(j,:),y(j,:));
            hold on;
    end
end

plot_pdf = figure(1);
xlabel('x')
ylabel('y')
title('\textbf{Phase portrait of the Nonlinear System}')
% legend('x','y','location','southwest')
box on
grid minor
hold off;


% Save pdf
set(plot_pdf, 'Units', 'Centimeters');
pos = get(plot_pdf, 'Position');
set(plot_pdf, 'PaperPositionMode', 'Auto', 'PaperUnits', 'Centimeters', ...
    'PaperSize',[pos(3), pos(4)]);
print(plot_pdf, 'exercise3_phase_portrait.pdf', '-dpdf', '-r0');

% Save png
print(gcf,'exercise3_phase_portrait.png','-dpng','-r600');

