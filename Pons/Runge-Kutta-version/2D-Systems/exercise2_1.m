%% Exercise 2

%-------------------------------------------------------------------------%
% Two-dimensional linear systems
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
b = -2;
c = 3;
d = 1;

% Numerical data
h = 0.001; % Time steps of RK4
t_final = 10; % time units
N_steps = ceil(t_final/h); % Number of steps, rounds up with ceil

% Declaration of solution vectors
t = 0:h:t_final;
x = zeros(1, length(t));
y = zeros(1, length(t));

% Initial conditions
t(1) = 0; % We begin at t=0 s
x(1) = 0.2;
y(1) = +1.5;

% Function handle
f1 = @(t,x,y) a*x + b*y;
f2 = @(t,x,y) c*x + d*y;


% RK4 update loop
for i=1:N_steps
    
    % Update function f1 and f2
    % Compute coefficients sub 1
    k1_f1 = f1(t(i)         ,x(i)                   ,y(i)               );
    k1_f2 = f2(t(i)         ,x(i)                   ,y(i)               );
    
    % Compute coefficients sub 2
    k2_f1 = f1(t(i) + h/2   ,x(i) + k1_f1*h/2       ,y(i) + k1_f2*h/2   );
    k2_f2 = f2(t(i) + h/2   ,x(i) + k1_f1*h/2       ,y(i) + k1_f2*h/2   );
    
    % Compute coefficients sub 3
    k3_f1 = f1(t(i) + h/2   ,x(i) + k2_f1*h/2       ,y(i) + k2_f2*h/2   );
    k3_f2 = f2(t(i) + h/2   ,x(i) + k2_f1*h/2       ,y(i) + k2_f2*h/2   );
    
    % Compute coefficients sub 3
    k4_f1 = f1(t(i) + h     ,x(i) + k3_f1*h         ,y(i) + k3_f2*h     );
    k4_f2 = f2(t(i) + h     ,x(i) + k3_f1*h         ,y(i) + k3_f2*h     );
    
    % Compute next step
    x(i+1) = x(i) + h/6*(k1_f1 + 2*k2_f1 + 2*k3_f1 + k4_f1);
    y(i+1) = y(i) + h/6*(k1_f2 + 2*k2_f2 + 2*k3_f2 + k4_f2);
    
end

% Perturbation
plot_pdf = figure(1);
plot(t,x,t,y);
xlabel('Time units')
ylabel('Value')
title('\textbf{Two-dimensional linear systems}')
legend('x','y','location','southwest')
box on
grid minor
hold off;


% Save pdf
set(plot_pdf, 'Units', 'Centimeters');
pos = get(plot_pdf, 'Position');
set(plot_pdf, 'PaperPositionMode', 'Auto', 'PaperUnits', 'Centimeters', ...
    'PaperSize',[pos(3), pos(4)]);
print(plot_pdf, 'exercise2_1_2DLinearSystem.pdf', '-dpdf', '-r0');

% Save png
print(gcf,'exercise2_1_2DLinearSystem.png','-dpng','-r600');



plot_pdf2 = figure(2);
plot(x,y);
xlabel('x')
ylabel('y')
title('\textbf{Two-dimensional linear systems}')
box on
grid minor

% Save pdf
set(plot_pdf2, 'Units', 'Centimeters');
pos = get(plot_pdf2, 'Position');
set(plot_pdf2, 'PaperPositionMode', 'Auto', 'PaperUnits', 'Centimeters', ...
    'PaperSize',[pos(3), pos(4)]);
print(plot_pdf2, 'exercise2_1_2DLinearSystem_xy.pdf', '-dpdf', '-r0');

% Save png
print(gcf,'exercise2_1_2DLinearSystem_xy.png','-dpng','-r600');

