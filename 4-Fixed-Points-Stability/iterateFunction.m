function x = iterateFunction(x0, logisticMap, n)

% Define solution vector
x = zeros(n,1);
% Set initial condition (each R)
x(1) = x0;

% Iterate logistic map function
for i = 2:n
    x(i) = f1(x(i-1));
end

end