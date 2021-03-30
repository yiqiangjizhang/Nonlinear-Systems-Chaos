%% Exercise 10

%-------------------------------------------------------------------------%
% Delayed logistic equation
%-------------------------------------------------------------------------%

% Date: 12/03/2021
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

%% Delayed Adler equation

% Paramter tau
tau = [1];

% For different tau values
for i=1:1:length(tau)
    solve_delay(tau(i));
end


% Plot
plot_pdf = figure(1);
xlabel('Time units')
ylabel('$y$')
title('\textbf{Delayed logistic equation for $\tau = 1$}')
% legend()
box on
grid minor

% Save pdf
set(plot_pdf, 'Units', 'Centimeters');
pos = get(plot_pdf, 'Position');
set(plot_pdf, 'PaperPositionMode', 'Auto', 'PaperUnits', 'Centimeters', ...
    'PaperSize',[pos(3), pos(4)]);
print(plot_pdf, 'delayed_logistic_eq.pdf', '-dpdf', '-r0');

% Save png
print(gcf,'delayed_logistic_eq.png','-dpng','-r600');


% Function definition
function solve_delay(tau)
    ic = [0.5];
    tspan = [0 100];
    lambda = 1.8;
    sol = dde23(@f,tau,ic,tspan);
    plot(sol.x,sol.y(1,:),'r-')
    hold on;
    function v=f(t,y,Z)
        v = [lambda*y(1).*(1-Z(1))];
    end
end
