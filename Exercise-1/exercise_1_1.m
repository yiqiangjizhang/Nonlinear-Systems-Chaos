%% Exercise: Runge Kutta 4th Order

%-------------------------------------------------------------------------%
% Problem 1: Program that solves a system of 2 ODE using RK4
%-------------------------------------------------------------------------%

% Date: 25/01/2021
% Author/s: Yi Qiang Ji Zhang
% Subject: Nonlinear Systems, Chaos and Control in Engineering
% Professor: Antonio Pons & Cristina Masoller

% Clear workspace, command window and close windows
clear;
close all;
clc;

r = 0;

it = 100; % Number of iterations for each r
r_N = 100; % Number of points tried in 'r'

x = zeros(r_N,it); % Initialize vector x


vector_r = linspace(1,4,r_N);

% Initial condition
x(1:length(vector_r),1) = 0.2;


for j = 1:length(vector_r)
    r = vector_r(j); % Get the current r
    
    for i = 1:it
        x(j,i+1) = r*x(j, i)*(1-x(j,i));
    end
        
end

% figure(1)
% plot(x,'.');
figure(2)
for i=1:r_N
    hold on
    plot(vector_r(i),x(i,:),'.k')
end
