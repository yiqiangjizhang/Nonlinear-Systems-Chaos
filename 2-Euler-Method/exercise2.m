%% Exercise 2: Euler Approximation

%-------------------------------------------------------------------------%
% First and Second order Euler method
%-------------------------------------------------------------------------%

% Date: 27/02/2021
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

%% Euler approximation

% 1.1 Numerical data
h = [1, 0.1, 0.01, 0.001, 0.0001]; % Time steps (dt)
t_final = 1; % time units

N_steps = zeros(1,length(h));
for j=1:length(h)
% N_steps(j) = ceil(t_final/h(j)); % Number of steps, rounds up with ceil
N_steps(j) = length(0:h(j):t_final)-1;
end

% Function handle
f = @(t,x) -x;

% Analytical expresssion
C = 1;
f_exact = @(dt) C*(exp(-dt));


% For each h (dt)
for j=1:length(h)
    
    t(j,1) = 0; % We begin at t=0 s
    x(j,1) = 1; % x(0) = 1
    
    % Euler method's update loop
    for i=1:N_steps(j)
        t(j,i+1) = t(j,i) + h(j);
        x(j,i+1) = x(j,i) + f(t(j,i),x(j,i))*h(j);
    end
    
    % Error calculation
    dt = t(j,end);
    sol_exact = f_exact(dt);
    disp(sol_exact);
    error(j) = abs(sol_exact-x(j,end));
    disp(error);
    plot_pdf = figure(1);
    plot(t(j,:),x(j,:));
    hold on;

end
xlabel('Time units')
ylabel('x')
title('Plot')
legend('dt = 1', 'dt = 0.1', 'dt = 0.01', 'dt = 0.001', 'dt = 0.0001')
box on
grid minor
hold off;

% Save pdf
set(plot_pdf, 'Units', 'Centimeters');
pos = get(plot_pdf, 'Position');
set(plot_pdf, 'PaperPositionMode', 'Auto', 'PaperUnits', 'Centimeters', ...
    'PaperSize',[pos(3), pos(4)]);
print(plot_pdf, 'euler_approx_dt.pdf', '-dpdf', '-r0');

% Save png
print(gcf,'euler_approx_dt.png','-dpng','-r600');


plot_pdf2 = figure(2);
loglog(h,error);
xlabel('Error')
ylabel('dt')
title('\textbf{Plot}')
box on
grid minor


% Save pdf
set(plot_pdf2, 'Units', 'Centimeters');
pos = get(plot_pdf2, 'Position');
set(plot_pdf2, 'PaperPositionMode', 'Auto', 'PaperUnits', 'Centimeters', ...
    'PaperSize',[pos(3), pos(4)]);
print(plot_pdf2, 'dt_vs_error.pdf', '-dpdf', '-r0');

% Save png
print(gcf,'dt_vs_error.png','-dpng','-r600');



