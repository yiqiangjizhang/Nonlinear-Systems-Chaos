function solve_delay1(tau)
ic = [0.5];
tspan = [0 100];
lambda = 1.8;
sol = dde23(@f,tau,ic,tspan);
plot(sol.x,sol.y(1,:),'r-')
hold on;
function v=f(t,y,Z)
    v = lambda*y(1).*(1-Z(1));
end
end