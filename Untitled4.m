% Logistics Map 
%       Classic chaos example. Plots semi-stable values of
%       x(n+1) = r*x(n)*(1-x(n)) as r increases to 4.
%
% Michael Hanchak, Dayton OH, USA, 2011
clear
scale = 10000; % determines the level of rounding
maxpoints = 200; % determines maximum values to plot
N = 3000; % number of "r" values to simulate
a = 2.0; % starting value of "r"
b = 4; % final value of "r"... anything higher diverges.
rs = linspace(a,b,N); % vector of "r" values
M = 500; % number of iterations of logistics equation
% Loop through the "r" values
for j = 1:length(rs)
    
    r=rs(j); % get current "r"
    x=zeros(M,1); % allocate memory
    x(1) = 0.5; % initial condition (can be anything from 0 to 1)
    
    for i = 2:M, % iterate
        x(i) = r*x(i-1)*(1-x(i-1));
    end
    % only save those unique, semi-stable values
    %out{j} = unique(round(scale*x(end-maxpoints:end)));
end
% Rearrange cell array into a large n-by-2 vector for plotting
data = [];
for k = 1:length(rs)
    n = length(out{k});
    data = [data;  rs(k)*ones(n,1),out{k}];
end
% Plot the data
figure(97);clf
h=plot(data(:,1),data(:,2)/scale,'k.');
set(h,'markersize',1)
axis tight
set(gca,'units','normalized','position',[0 0 1 1])
set(gcf,'color','white')
axis off