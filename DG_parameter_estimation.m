% gradient decent for finding optimal parmater alpha

u2=mean_exit_time(1.5); %observation  

alpha=1;          % initial guess
delta_alp=0.0001;

u1=mean_exit_time(alpha); 
b=alpha+delta_alp;
U1=mean_exit_time(b);


while norm(u1-u2)>10^(-5) && alpha<2 && alpha>0
    parameter = alpha
    alpha=alpha - 0.01* (norm(U1-u2)-norm(u1-u2))/delta_alp;
    u1=mean_exit_time(alpha);
    b=alpha+delta_alp;
    U1=mean_exit_time(b);
    
end




%-----------------------------------------------

u2=mean_exit_time(1.5); %observation  

for n=1:19
alpha=0.1*n;          % initial guess
u1=mean_exit_time(alpha); 
y(n)=norm(u1-u2);
end
plot(0.1:0.1:1.9,y)