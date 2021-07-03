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



% Numerical data
h = 0.001; % Time steps of RK4
t_final = 10; % time units
N_steps = ceil(t_final/h); % Number of steps, rounds up with ceil

% Declaration of solution vectors
t = 0:h:t_final;
x = zeros(1, length(t));
y = zeros(1, length(t));

% Function handle
f1 = @(t,x,y) -x + 4*x^3;
f2 = @(t,x,y) -2*x;


% Initial conditions
x0 = [-0.45 -0.40 -0.30 -0.25 -0.10 0.10 0.25 0.30 0.40 0.45];
y0 = [-0.45 -0.40 -0.30 -0.25 -0.10 0.10 0.25 0.30 0.40 0.45];


% For each h (dt)
for j=1:length(x0)
    % Initial conditions
    t(1) = 0; % We begin at t=0 s
    x(j,1) = x0(j);
    y(j,1) = y0(j);
    
    % RK4 update loop
    for i=1:N_steps

        % Update function f1 and f2
        % Compute coefficients sub 1
        k1_f1 = f1(t(i)         ,x(j,i)                   ,y(j,i)               );
        k1_f2 = f2(t(i)         ,x(j,i)                   ,y(j,i)               );

        % Compute coefficients sub 2
        k2_f1 = f1(t(i) + h/2   ,x(j,i) + k1_f1*h/2       ,y(j,i) + k1_f2*h/2   );
        k2_f2 = f2(t(i) + h/2   ,x(j,i) + k1_f1*h/2       ,y(j,i) + k1_f2*h/2   );

        % Compute coefficients sub 3
        k3_f1 = f1(t(i) + h/2   ,x(j,i) + k2_f1*h/2       ,y(j,i) + k2_f2*h/2   );
        k3_f2 = f2(t(i) + h/2   ,x(j,i) + k2_f1*h/2       ,y(j,i) + k2_f2*h/2   );

        % Compute coefficients sub 3
        k4_f1 = f1(t(i) + h     ,x(j,i) + k3_f1*h         ,y(j,i) + k3_f2*h     );
        k4_f2 = f2(t(i) + h     ,x(j,i) + k3_f1*h         ,y(j,i) + k3_f2*h     );

        % Compute next step
        x(j,i+1) = x(j,i) + h/6*(k1_f1 + 2*k2_f1 + 2*k3_f1 + k4_f1);
        y(j,i+1) = y(j,i) + h/6*(k1_f2 + 2*k2_f2 + 2*k3_f2 + k4_f2);

    end
        plot_pdf = figure(1);
        plot(x(j,:),y(j,:));
        hold on;
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

