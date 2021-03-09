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
r = -1:0.01:1;
h = [0 0.01];

t0 = 0;
dt = [0.1]; % Time steps (dt)
t_final = 50; % time units

N_steps = length(t0:dt:t_final)-1;

f1 = @(t,x,r,h) r*x - x^2 + h;

t(1:length(r),1) = 0;
x(1:length(r),1) = 0.01;

% Loop through each h
for k = 1:length(h)
    % Loop through each r
    for j=1:length(r)
        % Loop through time steps
        for i=1:N_steps
            t(j,i+1) = t(j,i) + dt;
            x(j,i+1) = x(j,i) + f1(t(j,i),x(j,i),r(j),h(k))*dt;
        end
    x_r(k,j)=x(j,end);
    end

plot(r,x_r(k,:));
hold on
end


xlabel('r')
ylabel('x')
title('Plot')
legend()
box on
grid minor
hold off;
