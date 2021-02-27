close all 
clear all
avalues=2.8:0.01:4;
N=100; 
a=avalues; 
x=0.1;
X=zeros(N,length(a));
for n=1:0.3*N 
    x=a.*x.*(1-x); 
    X(n,:)=x;
end
figure (9), hold on
for n=0.3*N:N 
    x=a.*x.*(1-x); 
    X(n,:)=x;
    plot(a,x,'.','MarkerSize',0.01) 
    axis ([2.8 4 0 1]) 
end
hold off