%% Logistic map
%
%-------------------------------------------------------------------------%
% Exercise 1: Logistic map for r=3.5
%-------------------------------------------------------------------------%

% Date: 25/01/2021
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

%% Logistic map

r = 3.5; % Parameter r
it = 10; % Number of iterations
x = zeros(1,it); % Initialize vector x

% Initial condition
x(1) = 0.2;

% Iterate function
for i = 2:it
    x(i) = r*x(i-1)*(1-x(i-1));
end

plot_pdf = figure(1);
plot(x);
title("\textbf{Logistic map for} $r=3.5$");
xlabel("Iterations $i$");
ylabel("$x(i)$");
grid on;
grid minor;

% Save pdf
set(plot_pdf, 'Units', 'Centimeters');
pos = get(plot_pdf, 'Position');
set(plot_pdf, 'PaperPositionMode', 'Auto', 'PaperUnits', 'Centimeters', ...
    'PaperSize',[pos(3), pos(4)]);
print(plot_pdf, 'logistic_map_r35.pdf', '-dpdf', '-r0');

% Save png
print(gcf,'logistic_map_r35.png','-dpng','-r600');
