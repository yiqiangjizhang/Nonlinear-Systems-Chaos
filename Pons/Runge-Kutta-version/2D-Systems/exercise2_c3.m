%% Exercise 2

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


%% Comparison

alpha = 3; 
beta  = 0.5; 

load('Ex2_S1_IC1.mat');
x_S1_IC1 = alpha*x;
y_S1_IC1 = alpha*y;

load('Ex2_S1_IC2.mat');
x_S1_IC2 = beta*x;
y_S1_IC2 = beta*y;

load('Ex2_S2_IC1.mat');
x_S2_IC1 = alpha*x;
y_S2_IC1 = alpha*y;

load('Ex2_S2_IC2.mat');
x_S2_IC2 = beta*x;
y_S2_IC2 = beta*y;

x_S1_compare = x_S1_IC1 + x_S1_IC2;
y_S1_compare = y_S1_IC1 + y_S1_IC2;
x_S2_compare = x_S2_IC1 + x_S2_IC2;
y_S2_compare = y_S2_IC1 + y_S2_IC2;

load('Ex2_S1_IC3.mat');
x_S1_IC3 = x;
y_S1_IC3 = y;

load('Ex2_S2_IC3.mat');
x_S2_IC3 = x;
y_S2_IC3 = y;

% Time domain

plot_pdf = figure(1);
plot(t,x_S1_IC3,t,y_S1_IC3,t,x_S1_compare,t,y_S1_compare);
grid on;
grid minor;
box on;
title("\textbf{Nonlinear System 1}");
ylabel("Value");
xlabel("Time units");
legend("$x_3$","$y_3$","$\alpha x_1 + \beta x_2$","$\alpha y_1 + \beta y_2$", 'location', 'best');

% Save pdf
set(plot_pdf, 'Units', 'Centimeters');
pos = get(plot_pdf, 'Position');
set(plot_pdf, 'PaperPositionMode', 'Auto', 'PaperUnits', 'Centimeters', ...
    'PaperSize',[pos(3), pos(4)]);
print(plot_pdf, 'exercise2_system1.pdf', '-dpdf', '-r0');

% Save png
print(gcf,'exercise2_system1.png','-dpng','-r600');


plot_pdf2 = figure(2);
plot(t,x_S2_IC3,t,y_S2_IC3,t,x_S2_compare,t,y_S2_compare);
grid on;
grid minor;
box on;
title("\textbf{Nonlinear System 2}");
ylabel("Value");
xlabel("Time units");
legend("$x_3$","$y_3$","$\alpha x_1 + \beta x_2$","$\alpha y_1 + \beta y_2$", 'location', 'best');

% Save pdf
set(plot_pdf2, 'Units', 'Centimeters');
pos = get(plot_pdf2, 'Position');
set(plot_pdf2, 'PaperPositionMode', 'Auto', 'PaperUnits', 'Centimeters', ...
    'PaperSize',[pos(3), pos(4)]);
print(plot_pdf2, 'exercise2_system2.pdf', '-dpdf', '-r0');

% Save png
print(gcf,'exercise2_system2.png','-dpng','-r600');



% Phase portrait

plot_pdf3 = figure(3);
plot(x_S1_IC3,y_S1_IC3,x_S1_compare,y_S1_compare);
grid on;
grid minor;
box on;
title("\textbf{Nonlinear System 1}");
ylabel("y");
xlabel("x");
legend("Sum","Combination");


% Save pdf
set(plot_pdf3, 'Units', 'Centimeters');
pos = get(plot_pdf3, 'Position');
set(plot_pdf3, 'PaperPositionMode', 'Auto', 'PaperUnits', 'Centimeters', ...
    'PaperSize',[pos(3), pos(4)]);
print(plot_pdf3, 'exercise2_system1_combination.pdf', '-dpdf', '-r0');

% Save png
print(gcf,'exercise2_system1_combination.png','-dpng','-r600');


plot_pdf4 = figure(4);
plot(x_S2_IC3,y_S2_IC3,x_S2_compare,y_S2_compare);
grid on;
grid minor;
box on;
title("\textbf{Nonlinear System 2}");
ylabel("y");
xlabel("x");
legend("Sum","Combination");

% Save pdf
set(plot_pdf4, 'Units', 'Centimeters');
pos = get(plot_pdf4, 'Position');
set(plot_pdf4, 'PaperPositionMode', 'Auto', 'PaperUnits', 'Centimeters', ...
    'PaperSize',[pos(3), pos(4)]);
print(plot_pdf4, 'exercise2_system2_combination.pdf', '-dpdf', '-r0');

% Save png
print(gcf,'exercise2_system2_combination.png','-dpng','-r600');

