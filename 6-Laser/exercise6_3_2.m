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
r0 = [-1 -2 -3];
v = 0.1; % Variation rate


% Time domain
t0 = 0; % Initial time
dt = [0.01]; % Time steps (dt)
t_final = 1000; % time units
N_steps = length(t0:dt:t_final)-1; % Number of time steps

% Function handle
f1 = @(t,x,r) r*x - x^2;

% Initial conditions
t(1:length(r0),1) = 0;
x(1:length(r0),1) = 0.01;

% Loop through each h
for j = 1:length(r0)
    % Loop through each r
    r_x(1:length(r0),1) = r0(j);
    for i=1:N_steps
        t(j,i+1) = t(j,i) + dt;
        r_x(j,i+1) = r_x(j,i) + v*dt;
        x(j,i+1) = x(j,i) + f1(t(j,i),x(j,i),r_x(j,i))*dt;
    end
figure(1)
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

% Duplicate plot
a1 = gca;
f2 = figure;
a2 = copyobj(a1,f2);
xlabel('r')
ylabel('x')
title('Plot')
legend()
box on
xlim([-3 1])
ylim([-5e-3 15e-3])

% https://es.mathworks.com/matlabcentral/answers/110325-how-to-copy-a-plot
