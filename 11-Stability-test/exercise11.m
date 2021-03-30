%% Exercise 11

%-------------------------------------------------------------------------%
% Stability Diagram:
% Test the stability diagram by simulating the equation with D=0 
% and different values of the parameters c and tau.
%-------------------------------------------------------------------------%

% Date: 16/03/2021
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
c = [0.4, -0.4, -1.1, -1.1];
tau = [1, 25, 30, 0.5];

% For different tau values
for i=1:1:length(tau)
    solve_delay(tau(i),c(i));
end


% Plot
plot_pdf = figure(1);
xlabel('Time units')
ylabel('$x$')
title('\textbf{Stability Diagram}')
legend('$\mathrm{c} = +0.4$, $\tau = 1$', '$\mathrm{c} = -0.4$, $\tau = 25$', ...
    '$\mathrm{c} = -1.1$, $\tau = 30$', '$\mathrm{c} = -1.1$, $\tau = 0.5$', ...
    'location', 'east')
box on
grid minor

% Save pdf
set(plot_pdf, 'Units', 'Centimeters');
pos = get(plot_pdf, 'Position');
set(plot_pdf, 'PaperPositionMode', 'Auto', 'PaperUnits', 'Centimeters', ...
    'PaperSize',[pos(3), pos(4)]);
print(plot_pdf, 'stability_diagram.pdf', '-dpdf', '-r0');

% Save png
print(gcf,'stability_diagram.png','-dpng','-r600');


% Function definition
function solve_delay(tau,c)
    ic = [0.5];
    tspan = [0 500];
    lambda = 1.8;
    sol = dde23(@f,tau,ic,tspan);
    plot(sol.x,sol.y(1,:))
    hold on;
    function v=f(t,y,Z)
        v = [y(1)-y(1).^3 + c.*Z(1)];
    end
end
