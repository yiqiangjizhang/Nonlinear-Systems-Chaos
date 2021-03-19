%% Exercise 5: Neuron model
%
%-------------------------------------------------------------------------%
% Simulation of a Neuron model with different I
%-------------------------------------------------------------------------%

% Date: 02/03/2021
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

%% Neuron model

% Constants and Parameters
C = 10;
g_Na = 74;
I = -40:10:40;
V_1_2 = 1.5;
g_L = 19;
k = 16;
E_L = -67;
E_Na = 60;

% Functions
f2 = @(V) (1 + exp((V_1_2 - V)/k))^(-1);
f1 = @(t,V,I) I - g_L*(V - E_L) - g_Na*f2(V)*(V - E_Na);

% Numerical data
h = [0.001]; % Time steps (dt)
t_final = 10; % time units

for j=1:length(h)
% N_steps(j) = ceil(t_final/h(j)); % Number of steps, rounds up with ceil
N_steps(j) = length(0:h(j):t_final)-1;
end

% Initial conditions
initial_cond = -80;

% For each h (dt)
for j=1:length(h)
  
    t(j,1) = 0; % We begin at t=0 s
    for I_counter = 1:length(I)
    for k = 1:length(initial_cond)
    V(j,1) = initial_cond(k); % x(0) = 1
    
        % Euler method's update loop
        for i=1:N_steps(j)
            t(j,i+1) = t(j,i) + h(j);
            V(j,i+1) = V(j,i) + f1(t(j,i),V(j,i),I(I_counter))*h(j);
        end
    end
    plot_pdf = figure(1);
    plot(t(j,:),V(j,:));
    hold on;
    end
end
xlim([0 1])
xlabel('Time units')
ylabel('Membrane Potential $V(t)$ [$\mathrm{mV}$]')
title('\textbf{Neuron model for $V_0=-80$ and different values of I}')
legend('$I=-40$', '$I=-30$', '$I=-20$', '$I=-10$', '$I=0$', '$I=+10$', '$I=+20$', ...
    '$I=+30$', '$I=+40$', 'location','east');
box on
grid minor
hold off;

% Save pdf
set(plot_pdf, 'Units', 'Centimeters');
pos = get(plot_pdf, 'Position');
set(plot_pdf, 'PaperPositionMode', 'Auto', 'PaperUnits', 'Centimeters', ...
    'PaperSize',[pos(3), pos(4)]);
print(plot_pdf, 'neuron_model_V0_m80.pdf', '-dpdf', '-r0');

% Save png
print(gcf,'neuron_model_V0_m80.png','-dpng','-r600');


