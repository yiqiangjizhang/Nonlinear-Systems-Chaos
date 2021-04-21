%% Exercise 4

%-------------------------------------------------------------------------%
% Direction Field
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

%% Map

% Field
x = -2.5:0.01:2.5;
y = -2.5:0.01:2.5;

dxdt = zeros(1,length(x));
dydt = zeros(1,length(y));

% Lines
i = 1;
while i <= length(x)
    dxdt(i) = x(i)+exp(0);
    dydt(i) = 0;
    i = i+1;
end

% Position of the quiver vectors
pos_x = size(length(x),length(y));
pos_y = size(length(x),length(y));


for i = 1:20:length(x)
   for j = 1:20:length(y) 
       pos_x(i,j) = x(i);
       pos_y(i,j) = y(j);
       value_dx_dt(i,j) = x(i) + exp(-y(j));
       value_dy_dt(i,j) = -y(j);        
   end
end

% Figure plot
plot_pdf = figure(1);
plot(x,dxdt,'b',y,dydt,'r');
hold on
qv = quiver(pos_x,pos_y,value_dx_dt,value_dy_dt,'c');
set(qv,'AutoScale','on', 'AutoScaleFactor', 100)
xlim([-2.5 1.5])
ylim([-1.5 1.5])
xlabel('x')
ylabel('y')
legend("dx/dt = 0","dy/dt = 0", 'location', 'northeast');
title('\textbf{Direction Field}');
grid on;
grid minor;
box on;

% Save pdf
set(plot_pdf, 'Units', 'Centimeters');
pos = get(plot_pdf, 'Position');
set(plot_pdf, 'PaperPositionMode', 'Auto', 'PaperUnits', 'Centimeters', ...
    'PaperSize',[pos(3), pos(4)]);
print(plot_pdf, 'direction_field.pdf', '-dpdf', '-r0');

% Save png
print(gcf,'direction_field.png','-dpng','-r600');

%% Trajectories


% Numerical data
h = 0.001; % Time steps of RK4
t_final = 2; % time units
N_steps = ceil(t_final/h); % Number of steps, rounds up with ceil

% Declaration of solution vectors
t = 0:h:t_final;
x = zeros(1, length(t));
y = zeros(1, length(t));


% Function handle
f1 = @(t,x,y) x + exp(-y);
f2 = @(t,x,y) -y;


% Initial conditions
x0 = [0 25 50 100];
y0 = [0 25 50 100];


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
        plot_pdf2 = figure(2);
        plot(t,x(j,:), 'r', t,y(j,:), 'b');
        hold on;
end


plot_pdf2 = figure(2);
xlabel('Time units')
ylabel('Value')
title('\textbf{Trajectories with various initial conditions}')
legend('x','y','location','best')
box on
grid minor
hold off;


% Save pdf
set(plot_pdf2, 'Units', 'Centimeters');
pos = get(plot_pdf2, 'Position');
set(plot_pdf2, 'PaperPositionMode', 'Auto', 'PaperUnits', 'Centimeters', ...
    'PaperSize',[pos(3), pos(4)]);
print(plot_pdf2, 'trajectories.pdf', '-dpdf', '-r0');

% Save png
print(gcf,'trajectories.png','-dpng','-r600');






