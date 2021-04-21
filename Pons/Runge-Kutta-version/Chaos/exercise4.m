%% Exercise 4

%-------------------------------------------------------------------------%
% Chaotic systems
%-------------------------------------------------------------------------%

% Date: 10/04/2021
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

% Constants
a = 0:0.1:0.4;
b = 2;
c = 4;

% Numerical data
h = 0.001; % Time steps of RK4
t_final = 50; % time units
N_steps = ceil(t_final/h); % Number of steps, rounds up with ceil

% Declaration of solution vectors
t = 0:h:t_final;
x = zeros(1, length(t));
y = zeros(1, length(t));
z = zeros(1, length(t));

% Function handle
f1 = @(t,x,y,z,a) -y-z;
f2 = @(t,x,y,z,a) x + a*y;
f3 = @(t,x,y,z,a) b + z*(x - c);


% Initial conditions
x0 = [1];
y0 = [1];
z0 = [1];

% For each a
for l=1:length(a)
    % For each x0
    for j=1:length(x0)
        % Initial conditions
        t(1) = 0; % We begin at t=0 s
        x(j,1) = x0(j);
        y(j,1) = y0(j);
        z(j,1) = z0(j);

        % RK4 update loop
        for i=1:N_steps

            % Update function f1 and f2
            % Compute coefficients sub 1
            k1_f1 = f1(t(i)         ,x(j,i)                   ,y(j,i)               ,z(j,i), a(l));
            k1_f2 = f2(t(i)         ,x(j,i)                   ,y(j,i)               ,z(j,i), a(l));
            k1_f3 = f3(t(i)         ,x(j,i)                   ,y(j,i)               ,z(j,i), a(l));

            % Compute coefficients sub 2
            k2_f1 = f1(t(i) + h/2   ,x(j,i) + k1_f1*h/2       ,y(j,i) + k1_f2*h/2   ,z(j,i) + k1_f3*h/2, a(l));
            k2_f2 = f2(t(i) + h/2   ,x(j,i) + k1_f1*h/2       ,y(j,i) + k1_f2*h/2   ,z(j,i) + k1_f3*h/2, a(l));
            k2_f3 = f3(t(i) + h/2   ,x(j,i) + k1_f1*h/2       ,y(j,i) + k1_f2*h/2   ,z(j,i) + k1_f3*h/2, a(l));

            % Compute coefficients sub 3
            k3_f1 = f1(t(i) + h/2   ,x(j,i) + k2_f1*h/2       ,y(j,i) + k2_f2*h/2   ,z(j,i) + k2_f3*h/2, a(l));
            k3_f2 = f2(t(i) + h/2   ,x(j,i) + k2_f1*h/2       ,y(j,i) + k2_f2*h/2   ,z(j,i) + k2_f3*h/2, a(l));
            k3_f3 = f3(t(i) + h/2   ,x(j,i) + k2_f1*h/2       ,y(j,i) + k2_f2*h/2   ,z(j,i) + k2_f3*h/2, a(l));
            % Compute coefficients sub 4
            k4_f1 = f1(t(i) + h     ,x(j,i) + k3_f1*h         ,y(j,i) + k3_f2*h     ,z(j,i) + k3_f3*h, a(l));
            k4_f2 = f2(t(i) + h     ,x(j,i) + k3_f1*h         ,y(j,i) + k3_f2*h     ,z(j,i) + k3_f3*h, a(l));
            k4_f3 = f3(t(i) + h     ,x(j,i) + k3_f1*h         ,y(j,i) + k3_f2*h     ,z(j,i) + k3_f3*h, a(l));

            % Compute next step
            x(j,i+1) = x(j,i) + h/6*(k1_f1 + 2*k2_f1 + 2*k3_f1 + k4_f1);
            y(j,i+1) = y(j,i) + h/6*(k1_f2 + 2*k2_f2 + 2*k3_f2 + k4_f2);
            z(j,i+1) = z(j,i) + h/6*(k1_f3 + 2*k2_f3 + 2*k3_f3 + k4_f3);

        end
        figure(1)
        plot(t, x(j,:));
        hold on;
        figure(2)
        plot(t, y(j,:));
        hold on;
        figure(3)
        plot(t, z(j,:));
        hold on;
        figure(4)
        plot(t, x(j,:), t, y(j,:),t, y(j,:));
        hold on;
        plot_pdf5 = figure(5);
        plot3(x,y,z);
        hold on;
        plot_pdf6 = figure(6);
        pspectrum(x(j,:),t)
        hold on;
        plot_pdf6 = figure(7);
        pspectrum(y(j,:),t)
        hold on;
        plot_pdf6 = figure(8);
        pspectrum(z(j,:),t)
        hold on;
    end
end


plot_pdf = figure(1);
xlabel('Time units')
ylabel('x')
title('\textbf{R\"ossler System $x$ vs $t$}')
legend('$a=0$','$a=0.1$','$a=0.2$','$a=0.3$','$a=0.4$','location','best')
box on
grid minor
hold off;

% Save pdf
set(plot_pdf, 'Units', 'Centimeters');
pos = get(plot_pdf, 'Position');
set(plot_pdf, 'PaperPositionMode', 'Auto', 'PaperUnits', 'Centimeters', ...
    'PaperSize',[pos(3), pos(4)]);
print(plot_pdf, 'rossler_system_x.pdf', '-dpdf', '-r0');

% Save png
print(gcf,'rossler_system_x.png','-dpng','-r600');


plot_pdf2 = figure(2);
xlabel('Time units')
ylabel('y')
title('\textbf{R\"ossler System $y$ vs $t$}')
legend('$a=0$','$a=0.1$','$a=0.2$','$a=0.3$','$a=0.4$','location','best')
box on
grid minor
hold off;

% Save pdf
set(plot_pdf2, 'Units', 'Centimeters');
pos = get(plot_pdf2, 'Position');
set(plot_pdf2, 'PaperPositionMode', 'Auto', 'PaperUnits', 'Centimeters', ...
    'PaperSize',[pos(3), pos(4)]);
print(plot_pdf2, 'rossler_system_y.pdf', '-dpdf', '-r0');

% Save png
print(gcf,'rossler_system_y.png','-dpng','-r600');

plot_pdf3 = figure(3);
xlabel('Time units')
ylabel('z')
title('\textbf{R\"ossler System $z$ vs $t$}')
legend('$a=0$','$a=0.1$','$a=0.2$','$a=0.3$','$a=0.4$','location','best')
box on
grid minor
hold off;

% Save pdf
set(plot_pdf3, 'Units', 'Centimeters');
pos = get(plot_pdf3, 'Position');
set(plot_pdf3, 'PaperPositionMode', 'Auto', 'PaperUnits', 'Centimeters', ...
    'PaperSize',[pos(3), pos(4)]);
print(plot_pdf3, 'rossler_system_z.pdf', '-dpdf', '-r0');

% Save png
print(gcf,'rossler_system_z.png','-dpng','-r600');

plot_pdf4 = figure(4);
xlabel('Time units')
ylabel('x, y, z')
title('\textbf{R\"ossler System $x, y, z$ vs $t$}')
legend('$a=0$','$a=0.1$','$a=0.2$','$a=0.3$','$a=0.4$','location','best')
box on
grid minor
hold off;

% Save pdf
set(plot_pdf4, 'Units', 'Centimeters');
pos = get(plot_pdf4, 'Position');
set(plot_pdf4, 'PaperPositionMode', 'Auto', 'PaperUnits', 'Centimeters', ...
    'PaperSize',[pos(3), pos(4)]);
print(plot_pdf4, 'rossler_system_xyz.pdf', '-dpdf', '-r0');

% Save png
print(gcf,'rossler_system_xyz.png','-dpng','-r600');


plot_pdf5 = figure(5);
xlabel('x')
ylabel('y')
zlabel('z')
legend('$a=0$','$a=0.1$','$a=0.2$','$a=0.3$','$a=0.4$','location','best')
title('\textbf{R\"ossler System}')
grid on;
grid minor

% Save pdf
set(plot_pdf5, 'Units', 'Centimeters');
pos = get(plot_pdf5, 'Position');
set(plot_pdf5, 'PaperPositionMode', 'Auto', 'PaperUnits', 'Centimeters', ...
    'PaperSize',[pos(3), pos(4)]);
print(plot_pdf5, 'rossler_system_3D_all.pdf', '-dpdf', '-r0');

% Save png
print(gcf,'rossler_system_3D_all.png','-dpng','-r600');


plot_pdf6 = figure(6);
grid minor
legend('$a=0$','$a=0.1$','$a=0.2$','$a=0.3$','$a=0.4$','location','best')

% Save pdf
set(plot_pdf6, 'Units', 'Centimeters');
pos = get(plot_pdf6, 'Position');
set(plot_pdf6, 'PaperPositionMode', 'Auto', 'PaperUnits', 'Centimeters', ...
    'PaperSize',[pos(3), pos(4)]);
print(plot_pdf6, 'rossler_system_power_spec_x_all.pdf', '-dpdf', '-r0');

% Save png
print(gcf,'rossler_system_power_spec_x_all.png','-dpng','-r600');

plot_pdf7 = figure(7);
grid minor
legend('$a=0$','$a=0.1$','$a=0.2$','$a=0.3$','$a=0.4$','location','best')

% Save pdf
set(plot_pdf7, 'Units', 'Centimeters');
pos = get(plot_pdf7, 'Position');

set(plot_pdf7, 'PaperPositionMode', 'Auto', 'PaperUnits', 'Centimeters', ...
    'PaperSize',[pos(3), pos(4)]);
print(plot_pdf7, 'rossler_system_power_spec_y_all.pdf', '-dpdf', '-r0');

% Save png
print(gcf,'rossler_system_power_spec_y_all.png','-dpng','-r600');

plot_pdf8 = figure(8);
grid minor
legend('$a=0$','$a=0.1$','$a=0.2$','$a=0.3$','$a=0.4$','location','best')

% Save pdf
set(plot_pdf8, 'Units', 'Centimeters');
pos = get(plot_pdf8, 'Position');

set(plot_pdf8, 'PaperPositionMode', 'Auto', 'PaperUnits', 'Centimeters', ...
    'PaperSize',[pos(3), pos(4)]);
print(plot_pdf8, 'rossler_system_power_spec_z_all.pdf', '-dpdf', '-r0');

% Save png
print(gcf,'rossler_system_power_spec_z_all.png','-dpng','-r600');






