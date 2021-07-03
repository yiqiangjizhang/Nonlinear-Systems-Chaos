%% Exercise 2: Euler Approximation

%-------------------------------------------------------------------------%
% First and Second order Euler method (not working)
%-------------------------------------------------------------------------%

% Date: 27/02/2021
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

%% Euler approximation

% 1.1 Numerical data
h = 0.1; % Time steps (dt)
t_final = 1; % time units
N_steps = ceil(t_final/h); % Number of steps, rounds up with ceil

% 1.3. Declaration of solution vectors
t = zeros(1, N_steps);
x = zeros(1, length(t));

% 1.4. Initial conditions
t(1) = 0; % We begin at t=0 s
x(1) = 1; % x(0) = 1

% Function handle
f = @(t,x) -x;

% Euler method's update loop
for i=1:N_steps
    t(i+1) = t(i) + h;
    x(i+1) = x(i) + f(t(i),x(i))*h;
end

% Exponential error
a = 1;
time = 1;
expo = a*(exp(-time));

error = abs(expo-x(end));
disp(error);


plot_pdf = figure(1);
plot(t,x);
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
xlabel('Time units')
ylabel('x')
title('Plot')
%legend('u','v')
box on
grid minor
hold off;


% Save pdf
set(plot_pdf, 'Units', 'Centimeters');
pos = get(plot_pdf, 'Position');
set(plot_pdf, 'PaperPositionMode', 'Auto', 'PaperUnits', 'Centimeters', ...
    'PaperSize',[pos(3), pos(4)]);
print(plot_pdf, 'euler_approx_dt.pdf', '-dpdf', '-r0');

% Save png
print(gcf,'euler_approx_dt.png','-dpng','-r600');
