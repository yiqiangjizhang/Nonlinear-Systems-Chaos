%-------------------------------------------------------------------------%
% Exercise 2: Comparison
%-------------------------------------------------------------------------%

% Date: 14/04/2021
% Author/s: Iván Sermanoukian Molina
% Subject: Nonlinear systems, chaos and control in engineering
% Professor: Antonio Pons & Cristina Masoller

% Clear workspace, command window and close windows
clc
clear all
close all

% LaTeX configuration
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

%% Comparison

alpha = 1.9; 
beta  = 14; 

load('TDS_Ex2_S1_IC1.mat');
x_sol_S1_IC1 = alpha*x_sol;
y_sol_S1_IC1 = alpha*y_sol;

load('TDS_Ex2_S1_IC2.mat');
x_sol_S1_IC2 = beta*x_sol;
y_sol_S1_IC2 = beta*y_sol;

load('TDS_Ex2_S2_IC1.mat');
x_sol_S2_IC1 = alpha*x_sol;
y_sol_S2_IC1 = alpha*y_sol;

load('TDS_Ex2_S2_IC2.mat');
x_sol_S2_IC2 = beta*x_sol;
y_sol_S2_IC2 = beta*y_sol;

x_S1_compare = x_sol_S1_IC1 + x_sol_S1_IC2;
y_S1_compare = y_sol_S1_IC1 + y_sol_S1_IC2;
x_S2_compare = x_sol_S2_IC1 + x_sol_S2_IC2;
y_S2_compare = y_sol_S2_IC1 + y_sol_S2_IC2;

load('TDS_Ex2_S1_IC3.mat');
x_sol_S1_IC3 = x_sol;
y_sol_S1_IC3 = y_sol;

load('TDS_Ex2_S2_IC3.mat');
x_sol_S2_IC3 = x_sol;
y_sol_S2_IC3 = y_sol;

% Time domain

h = figure(1);
plot(t_sol,x_sol_S1_IC3,t_sol,y_sol_S1_IC3,t_sol,x_S1_compare,t_sol,y_S1_compare);
grid on;
grid minor;
box on;
tlt = title("System 1");
ylab = ylabel("variable");
xlab = xlabel("time");
leg  = legend("$x_3$","$y_3$","$\alpha x_1 + \beta x_2$","$\alpha y_1 + \beta y_2$", 'location', 'best');
set(tlt,'fontsize',16);
set(xlab,'fontsize',16);
set(ylab,'fontsize',16);
set(leg,'fontsize',16);

h2 = figure(2);
plot(t_sol,x_sol_S2_IC3,t_sol,y_sol_S2_IC3,t_sol,x_S2_compare,t_sol,y_S2_compare);
grid on;
grid minor;
box on;
tlt = title("System 2");
ylab = ylabel("variable");
xlab = xlabel("time");
leg  = legend("$x_3$","$y_3$","$\alpha x_1 + \beta x_2$","$\alpha y_1 + \beta y_2$", 'location', 'best');
set(tlt,'fontsize',16);
set(xlab,'fontsize',16);
set(ylab,'fontsize',16);
set(leg,'fontsize',16);

% Phase portrait

h3 = figure(3);
plot(x_sol_S1_IC3,y_sol_S1_IC3,x_S1_compare,y_S1_compare);
grid on;
grid minor;
box on;
tlt = title("System 1");
ylab = ylabel("y");
xlab = xlabel("x");
set(tlt,'fontsize',16);
set(xlab,'fontsize',16);
set(ylab,'fontsize',16);

h4 = figure(4);
plot(x_sol_S2_IC3,y_sol_S2_IC3,x_S2_compare,y_S2_compare);
grid on;
grid minor;
box on;
tlt = title("System 2");
ylab = ylabel("y");
xlab = xlabel("x");
set(tlt,'fontsize',16);
set(xlab,'fontsize',16);
set(ylab,'fontsize',16);

%% Save data

% Save as pdf
% print2pdf(h,"Fig_2_1_txty_compare");
% print2pdf(h2,"Fig_2_2_txty_compare");
% print2pdf(h3,"Fig_2_1_xy_compare");
% print2pdf(h4,"Fig_2_2_xy_compare");
