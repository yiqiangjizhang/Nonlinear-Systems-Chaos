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


for tau=0:0.5:1
    solve_delay(tau);
end




xlabel('Time')
ylabel('$\theta$')
title('Plot')
legend()
box on
grid minor
hold off;

function solve_delay1(tau)
tau = 1;
ic = [0.5];
tspan = [0 100];
lambda = 1.8;
sol = dde23(@f,tau,ic,tspan);
plot(sol.x,sol.y(1,:),'r-')
hold on;
function v(t,y,Z);
        v = [lambda*y(1).*(1-Z(1))];
end
end



