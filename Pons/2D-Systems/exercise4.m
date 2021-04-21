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

% Function handle
f1 = @(t,x,y) x + exp(-y);
f2 = @(t,x,y) -y;

% 1.1 Numerical data
dt = 0.001; % Time steps (dt)
t_final = 2; % time units

for j=1:length(dt)
% N_steps(j) = ceil(t_final/h(j)); % Number of steps, rounds up with ceil
N_steps(j) = length(0:dt:t_final)-1;
end


% Initial conditions
x0 = [0 25 50 100];
y0 = [0 25 50 100];

% For each initial condition
for k=1:length(x0)
        t(1) = 0; % We begin at t=0 s
        x(k,1) = x0(k); % x(0)
        y(k,1) = y0(k); % y(0)
            % Euler method's update loop
            for i=1:N_steps(j)
                t(i+1) = t(i) + dt;
                x(k,i+1) = x(k,i) + f1(t(i),x(k,i),y(k,i))*dt;
                y(k,i+1) = y(k,i) + f2(t(i),x(k,i),y(k,i))*dt;
            end
            plot_pdf2 = figure(2);
            plot(t,x(k,:), 'r', t,y(k,:), 'b');
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






