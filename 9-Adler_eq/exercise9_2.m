%% Exercise 9

%-------------------------------------------------------------------------%
% Adler Equation
%-------------------------------------------------------------------------%

% Date: 09/03/2021
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

%% Plot 1

% Constants and Parameters
w = 1;
a = 0:0.01:0.99;

% Function handle
f1 = @(t,x,a,w) w-a*sin(x);

% 1.1 Numerical data
dt = [0.01]; % Time steps (dt)
t_final = 150; % time units

for j=1:length(dt)
% N_steps(j) = ceil(t_final/h(j)); % Number of steps, rounds up with ceil
N_steps(j) = length(0:dt:t_final)-1;
end

% Initial conditions
x0 = pi/2;

for j = 1:length(a)
    x(j,1) = x0; % x(0) = 1
    t(j,1) = 0; % We begin at t=0 s
% Euler method's update loop
    for i=1:N_steps
        t(j,i+1) = t(j,i) + dt;
        x(j,i+1) = x(j,i) + f1(t(j,i),x(j,i),a(j),w)*dt;
    end
T(j) = ((t(j,end) -t(j,1))/(x(j,end) - x(j,1)))*(2*pi);

figure(1)
plot(a(j),T(j),'o');
hold on;
end


xlabel('a')
ylabel('$T$')
title('Plot')
box on
grid minor

% Analytical
T_analytical = 2*pi./sqrt(w^2 - a.^2);
plot(a,T_analytical)


