%% Exercise 1: Logistic map

%-------------------------------------------------------------------------%
% Logistic map
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

%% Logistic map plot

% Declare function
logisticMap = @(x,r) r*x*(1-x);

% Number of iterations for each r
numIter = 400;
% Range of r
rMin = 3.5;
rMax = 4;
% Number of points of R
numR = 501;
% Initial condition
x0 = 0.5;

% Last 30% constant values (when it is in steady state)
k = 0.3*numIter;

% Create vector of R's
R = linspace(rMin, rMax, numR);
% Length of R
m = size(R, 2);

% Define points
points = zeros(numIter-k, m);

% Loop through each R
for i = 1:m
    % Calculate logistic map points 
    x = iterateFunction(x0, @(x) logisticMap(x,R(i)), numIter);
    % Matrix of points (for each colum R, get the points of logistic map)
    points(:,i) = x(k+1:end);
end

plot_pdf = figure(1);
hold on;
title("\textbf{Logistic map}");
plot(R, points, '.k');
box on
grid on
grid minor
xlabel("$r$");
ylabel("x(i)");
hold off;


% Save pdf
set(plot_pdf, 'Units', 'Centimeters');
pos = get(plot_pdf, 'Position');
set(plot_pdf, 'PaperPositionMode', 'Auto', 'PaperUnits', 'Centimeters', ...
    'PaperSize',[pos(3), pos(4)]);
print(plot_pdf, 'logistic_map_r_35_4.pdf', '-dpdf', '-r1000');

% Save png
print(gcf,'logistic_map_r_35_4.png','-dpng','-r600');
