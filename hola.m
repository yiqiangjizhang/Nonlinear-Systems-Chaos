clear;
close all;
clc;

set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

logisticMap = @(x,r) r*x*(1-x);

numIter = 400;
rMin = 1.5;
rMax = 4.5;
numC = 501;
x0 = 0.5;
k = 100;


R = linspace(rMin, rMax, numC);
m = size(R, 2);
points = zeros(numIter-k, m);

for i = 1:m
    x = iterateFunction(x0, @(x) logisticMap(x,R(i)), numIter);
    points(:,i) = x(k+1:end);
end

figure(1);
hold on;
title("\textbf{Logistic map}");
plot(R, points, '.k');
xlabel("$R$");
ylabel("Iterations");
hold off;