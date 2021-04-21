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

% Function handle
f1 = @(t,x,y,z,a) -y-z;
f2 = @(t,x,y,z,a) x + a*y;
f3 = @(t,x,y,z,a) b + z*(x - c);


% Numerical data
dt = [0.001]; % Time steps (dt)
t_final = 50; % time units

for j=1:length(dt)
% N_steps(j) = ceil(t_final/h(j)); % Number of steps, rounds up with ceil
N_steps(j) = length(0:dt:t_final)-1;
end

% Initial conditions
x0 = [1];
y0 = [1];
z0 = [1];

% For each a
for l=1:length(a)
    % For each initial condition
    for k=1:length(x0)
        % For each h (dt)
        for j=1:length(dt)
            t(j,1) = 0; % We begin at t=0 s
            x(j,1) = x0(k); % x(0)
            y(j,1) = y0(k); % y(0)
            z(j,1) = z0(k);
                % Euler method's update loop
                for i=1:N_steps(j)
                    t(j,i+1) = t(j,i) + dt(j);
                    x(j,i+1) = x(j,i) + f1(t(j,i),x(j,i),y(j,i),z(j,i),a(l))*dt(j);
                    y(j,i+1) = y(j,i) + f2(t(j,i),x(j,i),y(j,i),z(j,i),a(l))*dt(j);
                    z(j,i+1) = z(j,i) + f3(t(j,i),x(j,i),y(j,i),z(j,i),a(l))*dt(j);
                end
                plot_pdf = figure(1);
                plot(t(j,:), x(j,:));
                hold on;
                plot_pdf2 = figure(2);
                plot(t(j,:), y(j,:));
                hold on;
                plot_pdf3 = figure(3);
                plot(t(j,:), z(j,:));
                hold on;
                plot_pdf4 = figure(4);
                plot(t(j,:), x(j,:), t(j,:), y(j,:),t(j,:), y(j,:));
                hold on;
                plot_pdf5 = figure(5);
                plot3(x,y,z);
                hold on;
                plot_pdf6 = figure(6);
                pspectrum(x(k,:),t)
                hold on;
                plot_pdf6 = figure(7);
                pspectrum(y(k,:),t)
                hold on;
                plot_pdf6 = figure(8);
                pspectrum(z(k,:),t)
                hold on;
        end
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

% % Save pdf
% set(plot_pdf, 'Units', 'Centimeters');
% pos = get(plot_pdf, 'Position');
% set(plot_pdf, 'PaperPositionMode', 'Auto', 'PaperUnits', 'Centimeters', ...
%     'PaperSize',[pos(3), pos(4)]);
% print(plot_pdf, 'rossler_system_x.pdf', '-dpdf', '-r0');
% 
% % Save png
% print(gcf,'rossler_system_x.png','-dpng','-r600');


plot_pdf2 = figure(2);
xlabel('Time units')
ylabel('y')
title('\textbf{R\"ossler System $y$ vs $t$}')
legend('$a=0$','$a=0.1$','$a=0.2$','$a=0.3$','$a=0.4$','location','best')
box on
grid minor
hold off;

% % Save pdf
% set(plot_pdf2, 'Units', 'Centimeters');
% pos = get(plot_pdf2, 'Position');
% set(plot_pdf2, 'PaperPositionMode', 'Auto', 'PaperUnits', 'Centimeters', ...
%     'PaperSize',[pos(3), pos(4)]);
% print(plot_pdf2, 'rossler_system_y.pdf', '-dpdf', '-r0');
% 
% % Save png
% print(gcf,'rossler_system_y.png','-dpng','-r600');

plot_pdf3 = figure(3);
xlabel('Time units')
ylabel('z')
title('\textbf{R\"ossler System $z$ vs $t$}')
legend('$a=0$','$a=0.1$','$a=0.2$','$a=0.3$','$a=0.4$','location','best')
box on
grid minor
hold off;

% % Save pdf
% set(plot_pdf3, 'Units', 'Centimeters');
% pos = get(plot_pdf3, 'Position');
% set(plot_pdf3, 'PaperPositionMode', 'Auto', 'PaperUnits', 'Centimeters', ...
%     'PaperSize',[pos(3), pos(4)]);
% print(plot_pdf3, 'rossler_system_z.pdf', '-dpdf', '-r0');
% 
% % Save png
% print(gcf,'rossler_system_z.png','-dpng','-r600');

plot_pdf4 = figure(4);
xlabel('Time units')
ylabel('x, y, z')
title('\textbf{R\"ossler System $x, y, z$ vs $t$}')
legend('$a=0$','$a=0.1$','$a=0.2$','$a=0.3$','$a=0.4$','location','best')
box on
grid minor
hold off;

% % Save pdf
% set(plot_pdf4, 'Units', 'Centimeters');
% pos = get(plot_pdf4, 'Position');
% set(plot_pdf4, 'PaperPositionMode', 'Auto', 'PaperUnits', 'Centimeters', ...
%     'PaperSize',[pos(3), pos(4)]);
% print(plot_pdf4, 'rossler_system_xyz.pdf', '-dpdf', '-r0');
% 
% % Save png
% print(gcf,'rossler_system_xyz.png','-dpng','-r600');


plot_pdf5 = figure(5);
xlabel('x')
ylabel('y')
zlabel('z')
legend('$a=0$','$a=0.1$','$a=0.2$','$a=0.3$','$a=0.4$','location','best')
title('\textbf{R\"ossler System}')
grid on;
grid minor
% 
% % Save pdf
% set(plot_pdf5, 'Units', 'Centimeters');
% pos = get(plot_pdf5, 'Position');
% set(plot_pdf5, 'PaperPositionMode', 'Auto', 'PaperUnits', 'Centimeters', ...
%     'PaperSize',[pos(3), pos(4)]);
% print(plot_pdf5, 'rossler_system_3D_all.pdf', '-dpdf', '-r0');
% 
% % Save png
% print(gcf,'rossler_system_3D_all.png','-dpng','-r600');


plot_pdf6 = figure(6);
grid minor
legend('$a=0$','$a=0.1$','$a=0.2$','$a=0.3$','$a=0.4$','location','best')
% 
% % Save pdf
% set(plot_pdf6, 'Units', 'Centimeters');
% pos = get(plot_pdf6, 'Position');
% set(plot_pdf6, 'PaperPositionMode', 'Auto', 'PaperUnits', 'Centimeters', ...
%     'PaperSize',[pos(3), pos(4)]);
% print(plot_pdf6, 'rossler_system_power_spec_x_all.pdf', '-dpdf', '-r0');
% 
% % Save png
% print(gcf,'rossler_system_power_spec_x_all.png','-dpng','-r600');

plot_pdf7 = figure(7);
grid minor
legend('$a=0$','$a=0.1$','$a=0.2$','$a=0.3$','$a=0.4$','location','best')

% % Save pdf
% set(plot_pdf7, 'Units', 'Centimeters');
% pos = get(plot_pdf7, 'Position');
% 
% set(plot_pdf7, 'PaperPositionMode', 'Auto', 'PaperUnits', 'Centimeters', ...
%     'PaperSize',[pos(3), pos(4)]);
% print(plot_pdf7, 'rossler_system_power_spec_y_all.pdf', '-dpdf', '-r0');

% % Save png
% print(gcf,'rossler_system_power_spec_y_all.png','-dpng','-r600');

plot_pdf8 = figure(8);
grid minor
legend('$a=0$','$a=0.1$','$a=0.2$','$a=0.3$','$a=0.4$','location','best')

% % Save pdf
% set(plot_pdf8, 'Units', 'Centimeters');
% pos = get(plot_pdf8, 'Position');
% 
% set(plot_pdf8, 'PaperPositionMode', 'Auto', 'PaperUnits', 'Centimeters', ...
%     'PaperSize',[pos(3), pos(4)]);
% print(plot_pdf8, 'rossler_system_power_spec_z_all.pdf', '-dpdf', '-r0');
% 
% % Save png
% print(gcf,'rossler_system_power_spec_z_all.png','-dpng','-r600');






