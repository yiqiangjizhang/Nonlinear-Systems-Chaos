b=0.5
t=1:1:10
N_an(1)=5
c = -log(b*N_an(1))
for i=1:length(t)
N_an(i) = exp(-c*exp(-b*i))/b
end 
plot(t,N_an)