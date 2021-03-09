%% Exercise 6_3

%-------------------------------------------------------------------------%
% Simulate the equation with r increasing linearly in time. 
% Consider different variation rate (v) and/or different 
% initial value of the parameter (r0).
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

%% Plot 2

% Constants and Parameters
r0 = -1;
v = [0.1 0.01]; % Variation rate


% Time domain
t0 = 0; % Initial time
dt = [0.01]; % Time steps (dt)
t_final = 1000; % time units
N_steps = length(t0:dt:t_final)-1; % Number of time steps

% Function handle
f1 = @(t,x,r) r*x - x^2;

% Initial conditions
t(1:length(v),1) = 0;
x(1:length(v),1) = 0.01;
r_x(1:length(v),1) = r0;

% Loop through each h
for j = 1:length(v)
    % Loop through each r

    for i=1:N_steps
        r_x(j,i+1) = r_x(j,i) + v(j)*dt;
        t(j,i+1) = t(j,i) + dt;
        x(j,i+1) = x(j,i) + f1(t(j,i),x(j,i),r_x(j,i))*dt;
    end

plot(r_x(j,:),x(j,:));
hold on
end

xlim([-1 4])
xlabel('r')
ylabel('x')
title('Plot')
legend()
box on
grid minor
hold off;
