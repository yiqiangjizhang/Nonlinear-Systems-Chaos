%% Adler Equation

%-------------------------------------------------------------------------%
% Exercise 9 : With w=1/sqrt(2), calculate the trajectory for an arbitrary
% initial condition.
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

%% Plot 1

% Constants and Parameters
w = 1/sqrt(2);
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
x0 = [0:0.1:6.5];

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

xlim([0 15])
ylim([0 8])
set(gca,'YTick',0:pi/2:3*pi) 
set(gca,'YTickLabel',{'$04','$\pi/2$','$\pi$','$3\pi/2$','$2\pi$','$5\pi/2$'})
xlabel('Time units')
ylabel('$\theta$')
title('\textbf{Adler equation}')
% legend()
box on
grid minor
hold off;

% Save pdf
set(plot_pdf, 'Units', 'Centimeters');
pos = get(plot_pdf, 'Position');
set(plot_pdf, 'PaperPositionMode', 'Auto', 'PaperUnits', 'Centimeters', ...
    'PaperSize',[pos(3), pos(4)]);
print(plot_pdf, 'adler_eq_theta_time_init_cond.pdf', '-dpdf', '-r0');

% Save png
print(gcf,'adler_eq_theta_time_init_cond.png','-dpng','-r600');
