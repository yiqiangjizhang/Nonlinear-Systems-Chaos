%-------------------------------------------------------------------------%
% Exercise 4: Plot the direction field
%-------------------------------------------------------------------------%

% Date: 04/03/2021
% Author/s: Iván Sermanoukian Molina
% Subject: Nonlinear systems, chaos and control in engineering
% Professor: Antonio Pons & Cristina Masoller

% Clear workspace, command window and close windows
clc
clear all
close all

% LaTeX configuration
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

%% Map

x = -2.5:0.01:1.5;
y = -2.5:0.01:1.5;
dxdt = zeros(1,length(x));
dydt = zeros(1,length(y));
i = 1;
while i <= length(x)
    dxdt(i) = x(i)+exp(0);
    dydt(i) = 0;
    i = i+1;
end

% Position of the quiver vectors
position_x = size(length(x),length(y));
position_y = size(length(x),length(y));

for i = 1:20:length(x)
   for j = 1:20:length(y) 
       position_x(i,j) = x(i);
       position_y(i,j) = y(j);
       val_dx_dt(i,j) = x(i)+exp(-y(j));
       val_dy_dt(i,j) = -y(j);        
   end
end

% Figure plot
h = figure(1);
plot(x,dxdt,'b',y,dydt,'r');
hold on
qv = quiver(position_x,position_y,val_dx_dt,val_dy_dt,'c');
xlim([-2.5 1.5])
ylim([-1.5 1.5])
grid on;
grid minor;
box on;

tlt = title("Direction Field");
ylab = ylabel("y");
xlab = xlabel("x");
leg  = legend("dx/dt = 0","dy/dt = 0", 'location', 'southeast');
set(qv,'AutoScale','on', 'AutoScaleFactor', 100)
set(tlt,'fontsize',16);
set(xlab,'fontsize',16);
set(ylab,'fontsize',16);
set(leg,'fontsize',16);

%% Save data

% % Save as pdf
% print2pdf(h,"Fig_4");
