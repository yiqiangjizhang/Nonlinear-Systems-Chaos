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
f1 = @(x,r) r*x + x^3 - x^5;

r = -100:10:100;
x = -100:10:100;

x_dot = zeros(length(r),length(x));

for j=1:length(r)
    for i=1:length(x)
        x_dot(j,i) = f1(x(i),r(j));
    end
    plot(x,x_dot);
    hold on
end
title("\textbf{ODE plot}");
legend()
box on
grid on
grid minor
xlabel("$x$");
ylabel("$\dot{x}$");


figure(2)
plot(r,x)
% 
% % Number of iterations for each r
% numIter = 400;
% % Range of r
% rMin = 0;
% rMax = 4.5;
% % Number of points of R
% numR = 501;
% % Initial condition
% x0 = 0;
% 
% % Last 30% constant values (when it is in steady state)
% k = 0.3*numIter;
% 
% % Create vector of R's
% R = linspace(rMin, rMax, numR);
% % Length of R
% m = size(R, 2);
% 
% % Define points
% points = zeros(numIter-k, m);
% 
% % Loop through each R
% for i = 1:m
%     % Calculate logistic map points 
%     x = iterateFunction(x0, @(x) f1(x,R(i)), numIter);
%     % Matrix of points (for each colum R, get the points of logistic map)
%     points(:,i) = x(k+1:end);
% end
% 
% figure(1);
% hold on;
% title("\textbf{Logistic map}");
% plot(R, x, '.k');
% box on
% grid on
% grid minor
% xlabel("$r$");
% ylabel("Iterations");
% hold off;