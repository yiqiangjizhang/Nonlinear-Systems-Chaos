%% Matlab Introducxtion
%
%-------------------------------------------------------------------------%
% Problem 1: Program that solves a system of 2 ODE
%-------------------------------------------------------------------------%

% Date: 25/01/2021
% Author/s: Yi Qiang Ji Zhang
% Subject: Nonlinear Systems, Chaos and Control in Engineering
% Professor: Antonio Pons & Cristina Masoller

% Clear workspace, command window and close windows
clear;
close all;
clc;


r = 3.5;

it = 10; % Number of iterations

x = zeros(1,it); % Initialize vector x

% Initial condition
x(1) = 0.2;

for i = 2:it
    x(i) = r*x(i-1)*(1-x(i-1));
end

plot(x);



