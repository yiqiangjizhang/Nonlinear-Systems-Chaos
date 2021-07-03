%% Logistic map
%
%-------------------------------------------------------------------------%
% Exercise 1: Logistic map (NOT WORKING)
%-------------------------------------------------------------------------%

% Date: 25/01/2021
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

%% Logistic map

r = 0; % Parameter r

it = 100; % Number of iterations for each r
r_N = 100; % Number of points tried in 'r'

x = zeros(r_N,it); % Initialize vector x

vector_r = linspace(1,4,r_N); 

% Initial condition
x(1:length(vector_r),1) = 0.2;


for j = 1:length(vector_r)
    r = vector_r(j); % Get the current r
    
    for i = 1:it
        x(j,i+1) = r*x(j, i)*(1-x(j,i));
    end
        
end

% figure(1)
% plot(x,'.');
plot_pdf = figure(2);
for i=1:r_N
    hold on
    plot(vector_r(i),x(i,:),'.k')
end

title("\textbf{Logistic map}");
xlabel("$r$");
ylabel("Iterations");
grid on;
grid minor;

% % Save pdf
% set(plot_pdf, 'Units', 'Centimeters');
% pos = get(plot_pdf, 'Position');
% set(plot_pdf, 'PaperPositionMode', 'Auto', 'PaperUnits', 'Centimeters', ...
%     'PaperSize',[pos(3), pos(4)]);
% print(plot_pdf, 'logistic_map_r35.pdf', '-dpdf', '-r0');
% 
% % Save png
% print(gcf,'logistic_map_r35.png','-dpng','-r600');

