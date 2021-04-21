%% Mandelbrot set 

%-------------------------------------------------------------------------%
% Mandelbrot set
%-------------------------------------------------------------------------%

% Date: 09/03/2021
% Author/s: Group 1
% Subject: High Performance Computing
% Professor: Manel Soria & Arnau Miro

% Clear workspace, command window and close windows
clear;
close all;
clc;
% Set interpreter to latex
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');


%%
% 
% % X plane (real(c))
% x0 = -1.5; % x inferior limit
% xf = 1.5; % x superior limit
% dx = 0.01; % Step size
% 
% % X plane (real(c))
% y0 = -1.5;
% yf = 1.5;
% dy = 0.01; % Step size
% 
% % Vectors
% x = x0:dx:xf;
% y = y0:dy:yf;
% 
% % Create a 2D mesh
% [X,Y] = meshgrid(x,y);
% 
% % 1.5. Definition of Function Handles
% f1 = @(t,x,y) x-e^y; % First function
% f2 = @(t,x,y) -y; % Second function
% 
% 
% 
% 
% % For each h (dt)
% for j=1:length(dt)
%     t(j,1) = 0; % We begin at t=0 s
%     x(j,1) = x0; % x(0)
%     y(j,1) = y0; % y(0)
%         % Euler method's update loop
%         for i=1:N_steps(j)
%             t(j,i+1) = t(j,i) + dt(j);
%             x(j,i+1) = x(j,i) + f1(t(j,i),x(j,i),y(i,j))*dt(j);
%             y(j,i+1) = y(j,i) + f2(t(j,i),x(j,i),y(i,j))*dt(j);
%         end
%         plot(t(j,:),x(j,:));
%         hold on;
% end
% 
% 
% 
% 

% Quiver
% [X,Y] = meshgrid(-1.5:0.1:1.5,-1.5:0.1:1.5); % Creates a mesh of points

x = -2.5:0.01:1.5;
y = -2.5:0.01:1.5;

i=1;
while i <= length(x)
    dxdt(i) = x(i) + exp(0);
    dydt(i) = 0;
    i = i+1;
end

pos_x = size(length(x), length(y));
pos_y = size(length(x), length(y));

for i=1:20:length(x)
    for j=1:20:length(y)
        pos_x(i,j) = x(i);
        pos_y(i,j) = y(j);
        value_dxdt(i,j) = x(i) + exp(-y(j));
        value_dydt(i,j) = -y(j);
    end
end

% Plot quiver
quiver(value_dxdt, value_dydt,10);
hold on
plot(x,dxdt,y,dydt);
hold off;

% For each point (u,v) calculate the (du/dt,dv/dt)
% for i=1:size(X,1)
%     for j=1:size(Y,1)
%         value_dudt(i,j) = i+exp(-j);
%         value_dvdt(i,j) = -j;
% 
%     end
% end


% x = -1.5:0.1:1.5;
% y = -1.5:0.1:1.5;
% 
% for i=1:length(x)
%     for j=1:length(y)
%         value_dudt(i,j) = i+exp(-j);
%         value_dvdt(i,j) = -j;
%     end
% end
% 
% 
% % Plot quiver
% quiver(x_var,y_var, x,y);
% hold on;
% 
% 





