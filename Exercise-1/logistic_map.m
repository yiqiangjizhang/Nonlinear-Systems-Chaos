clear;
close all;
clc;

set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

logisticMap = @(x,c) c*x*(1-x);

numIter = 400;
cMin = 1.5;
cMax = 4.5;
numC = 501;
x0 = 0.5;
k = 100;


C = linspace(cMin, cMax, numC);
m = size(C, 2);
points = zeros(numIter-k, m);

for i = 1:m
    x = iterateFunction(x0, @(x) logisticMap(x,C(i)), numIter);
    points(:,i) = x(k+1:end);
end

figure(1);
hold on;
title("\textbf{Logistic map}");
plot(C, points, '.k');
xlabel("$C$");
ylabel("Iterations");
hold off;