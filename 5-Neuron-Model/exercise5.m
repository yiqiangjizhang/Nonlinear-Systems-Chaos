%% Exercise 3: Fixed Points and Stability
%
%-------------------------------------------------------------------------%
% First and Second order Euler method
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

% Constants and Parameters
C = 10;
g_Na = 74;
I = 0;
V_1_2 = 1.5;
g_L = 19;
k = 16;
E_L = -67;
E_Na = 60;


f2 = @(V) (1 + exp((V_1_2 - V)/k))^(-1);
f1 = @(t,V) I - g_L*(V - E_L) - g_Na*f2(V)*(V - E_Na);


% 1.1 Numerical data
h = [0.001]; % Time steps (dt)
t_final = 1; % time units

V_steps = zeros(1,length(h));
for j=1:length(h)
% N_steps(j) = ceil(t_final/h(j)); % Number of steps, rounds up with ceil
N_steps(j) = length(0:h(j):t_final)-1;
end

% Initial conditions
initial_cond = -60:10:40;

% For each h (dt)
for j=1:length(h)
    
    t(j,1) = 0; % We begin at t=0 s
    for k = 1:length(initial_cond)
    V(j,1) = initial_cond(k); % x(0) = 1
    
        % Euler method's update loop
        for i=1:N_steps(j)
            t(j,i+1) = t(j,i) + h(j);
            V(j,i+1) = V(j,i) + f1(t(j,i),V(j,i))*h(j);
        end

    plot(t(j,:),V(j,:));
    hold on;
    end
end
xlabel('Time')
ylabel('Numerical Solution')
title('Plot')
legend()
box on
grid minor
hold off;

