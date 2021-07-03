%% Exercise 4: Fixed points and stability analysis

%-------------------------------------------------------------------------%
% Find the fixed points and compute their stability
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

%% Fixed points and stability

% Limits
initial = -0.3;
final= 0.2;
r1 = initial:0.001:0;
r2 = -1/4:0.001:final;
r3 = -1/4:0.001:0;
r4 = 0:0.001:final;

% Fixed points
x1_1 = 0*r1;
x1_2 = 0*r4;
x2 = ((1+(1+4*r2).^(1/2))/2).^(1/2);
x3 = ((1+(1+4*r2).^(1/2))/2).^(1/2);
x4 = ((1-(1+4*r3).^(1/2))/2).^(1/2);
x5 = ((1-(1+4*r3).^(1/2))/2).^(1/2);


% Plot of subcritial Pitchfork Bifurcation
plot_pdf = figure(1);
plot(r1,x1_1);
hold on;
box on
grid on
grid minor
title("\textbf{Subcritical Pitchfork Bifurcation}");
xlabel("$r$");
ylabel("$x$");
plot(r2, x2, r2, x3);
plot(r3, x4, '--', r3, x5, '--');
plot(r4, x1_2, '--');

% Plot symmetry
plot(r3, -x4, '--', r3, -x5, '--');
plot(r2, -x2, r2, -x3);
hold off


% Save pdf
set(plot_pdf, 'Units', 'Centimeters');
pos = get(plot_pdf, 'Position');
set(plot_pdf, 'PaperPositionMode', 'Auto', 'PaperUnits', 'Centimeters', ...
    'PaperSize',[pos(3), pos(4)]);
print(plot_pdf, 'subcritical_bif_hysteresis.pdf', '-dpdf', '-r0');

% Save png
print(gcf,'subcritical_bif_hysteresis.png','-dpng','-r600');
