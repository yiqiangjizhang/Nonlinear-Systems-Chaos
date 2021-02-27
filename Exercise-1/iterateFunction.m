function x = iterateFunction(x0, logisticMap, n)

x = zeros(n,1);
x(1) = x0;
for i = 2:n
    x(i) = logisticMap(x(i-1));
end

end