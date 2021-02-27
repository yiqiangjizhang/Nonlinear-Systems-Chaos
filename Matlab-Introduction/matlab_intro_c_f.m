%% Matlab Introduction
% Plots c) and f)
%-------------------------------------------------------------------------%
% Problem 1: Program that solves a system of 2 ODE
%-------------------------------------------------------------------------%

% Date: 25/02/2021
% Author/s: Yi Qiang Ji Zhang
% Subject: Nonlinear Systems, Chaos and Control in Engineering
% Professor: Antonio Pons & Cristina Masoller

% Clear workspace, command window and close windows
clear;
close all;
clc;

% c) Oscillatory
% Constants
% The following constants were obtained using 2 points and the eq. of a line:
epsilon = 0.01;
alpha = 0.2;
I = 0.1;

%% Plot system of equations

% X axis
u = -1.5:0.1:1.5;
v = -1.5:0.1:1.5;
% Define derivatives
v_dudt = zeros(1,length(u));
v_dvdt = zeros(1,length(u));

% Solve equation when du/dt = 0 and dv/dt =0
i = 1;
while i <=length(u)
    v_dudt(i) = u(i)-u(i).^3;
    v_dvdt(i) = 1/alpha *(u(i)+I);
    i = i+1;
end

% Plot Figure
figure(1)
plot(u,v_dudt,'b',u,v_dvdt,'r');
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
ylim([-1 1]);
xlabel('u')
ylabel('v')
legend('du/dt=0','dv/dt=0')
title('Plot c) Oscillatory')
box on
grid minor
hold on;

% Quiver
[X,Y] = meshgrid(-1.5:0.1:1.5,-1.5:0.1:1.5); % Creates a mesh of points

% For each point (u,v) calculate the (du/dt,dv/dt)
for i=1:size(X,2)
    for j=1:size(Y,1)
        value_dudt(i,j) = u(i)-u(i).^3-v(j);
        value_dvdt(i,j) = epsilon*(u(i)-alpha*v(j)+I);
    end
end

% Plot quiver
quiver(Y,X,value_dudt,value_dvdt);
hold on;


%% ODE45 solver

% Time
tspan = linspace(0,1000,10e2);

% Function handle
funct = @(t,x) [x(1)-x(1).^3-x(2)
                epsilon*(x(1) - alpha*(x(2)) + I)];
            
% ODE45 solver
[t,x] = ode45(funct,tspan,[-1 0]);

% Figure
figure(2)
plot(t,x(:,1), t,x(:,2));
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
xlabel('Time')
ylabel('Concentration')
title('Plot e) (without perturbation)')
legend('u','v')
box on
grid minor

%% Runge Kutta 4th order solver

% System of equations
%{
du/dt = u-u^3-v;
dv/dt = epsilon*(u-alpha*v+I);
%}

% 1.1 Numerical data
h = 0.1; % Time steps of RK4
t_final = 1000; % time units
N_steps = ceil(t_final/h); % Number of steps, rounds up with ceil

% 1.3. Declaration of solution vectors
t = 0:h:t_final;
u = zeros(1, length(t));
v = zeros(1, length(t));

% 1.4. Initial conditions
t(1) = 0; % We begin at t=0 s
u(1) = -1;
v(1) = 0;

% 1.5. Definition of Function Handles
f1 = @(t,u,v) u-u^3-v; % First function
f2 = @(t,u,v) epsilon*(u-alpha*v+I); % Second function


% RK4 update loop
for i=1:N_steps
    
    % Update function f1 and f2
    % Compute coefficients sub 1
    k1_f1 = f1(t(i)         ,u(i)                   ,v(i)               );
    k1_f2 = f2(t(i)         ,u(i)                   ,v(i)               );
    
    % Compute coefficients sub 2
    k2_f1 = f1(t(i) + h/2   ,u(i) + k1_f1*h/2       ,v(i) + k1_f2*h/2   );
    k2_f2 = f2(t(i) + h/2   ,u(i) + k1_f1*h/2       ,v(i) + k1_f2*h/2   );
    
    % Compute coefficients sub 3
    k3_f1 = f1(t(i) + h/2   ,u(i) + k2_f1*h/2       ,v(i) + k2_f2*h/2   );
    k3_f2 = f2(t(i) + h/2   ,u(i) + k2_f1*h/2       ,v(i) + k2_f2*h/2   );
    
    % Compute coefficients sub 3
    k4_f1 = f1(t(i) + h     ,u(i) + k3_f1*h         ,v(i) + k3_f2*h     );
    k4_f2 = f2(t(i) + h     ,u(i) + k3_f1*h         ,v(i) + k3_f2*h     );
    
    % Compute next step
    u(i+1) = u(i) + h/6*(k1_f1 + 2*k2_f1 + 2*k3_f1 + k4_f1);
    v(i+1) = v(i) + h/6*(k1_f2 + 2*k2_f2 + 2*k3_f2 + k4_f2);
    
    % Perturbation
    if (i==t_final/h*0.4)
        u(i+1) = 1;
    end
end

figure(3);
plot(t,u, t,v);
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
xlabel('Time')
ylabel('Concentration')
title('Plot e) (with perturbation)')
legend('u','v')
box on
grid minor
hold off;













