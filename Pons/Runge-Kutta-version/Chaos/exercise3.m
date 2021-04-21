%% Exercise 3

%-------------------------------------------------------------------------%
% Henon Map
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

% Constants
a = 1.4;
b = 0.3;

% Initial conditions
x(1) = 0;
y(1) = 0;

% For 10000 iterations
N = 10000;

for i=1:N
    x(i+1) = y(i) + 1 - a*(x(i)^2);
    y(i+1) = b*x(i);
end

% Plot
plot_pdf = figure(1);
plot(x,y,'.');
xlabel('x')
ylabel('y')
title('\textbf{Henon Map}')
box on
grid minor
hold off;

% Save pdf
set(plot_pdf, 'Units', 'Centimeters');
pos = get(plot_pdf, 'Position');
set(plot_pdf, 'PaperPositionMode', 'Auto', 'PaperUnits', 'Centimeters', ...
    'PaperSize',[pos(3), pos(4)]);
print(plot_pdf, 'henon_map.pdf', '-dpdf', '-r0');

% Save png
print(gcf,'henon_map.png','-dpng','-r600');

