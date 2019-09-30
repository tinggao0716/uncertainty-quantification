% meant exit time from a bounded domain (-r,rr)
% levy motion
% here we didn't rescale the domain 



function U=escapenonsym(Alpha)


alpha=Alpha;
eps= 1;
d= 0;
f=@(x)  x-x.^3;
r= 0.1;              %-left
rr= 0.1;           %right



%for n=1:3
h=0.01;
J=rr/h;          %right
L=r/h;           %left
x=(-r-rr):h:(r+rr);

C=-zeta(alpha-1)*h^(2-alpha);          % correction term u''(x)
cons=alpha*gamma((1+alpha)/2)/(2^(1-alpha)*sqrt(pi)*gamma(1-alpha/2)); 

b=zeros(L+J-1,1);
a=zeros(L+J-1,1);
c=zeros(L+J-1,1); 

%nonintegral part
for i=J+2:L+2*J
    b(i-J-1)= -cons*2*C*eps/(h^2) - d/(h^2) - ...
        cons*eps*(1/(x(i)+r)^alpha+1/(rr-x(i))^alpha)/(alpha); % coefficient of u_j
    
    a(i-J-1)= cons*C*eps/(h^2) + d/(2*h^2) - f(x(i))/(2*h);  % coefficient of u_(j-1)
    c(i-J-1)= cons*C*eps/(h^2) + d/(2*h^2) + f(x(i))/(2*h);  % coefficient of u_(j+1)   
end
A=diag([b;1])+diag([a(2:end);0],-1)+diag(c(1:end),1);

if alpha<1 || alpha==1
 A(1,1)= cons*C*eps/(h^2) + d/(2*h^2)- cons*eps*(1/(x(J+2)+r)^alpha+1/(rr-x(J+2))^alpha)/(alpha) - f(x(J+2))/(h);
 A(1,2)= -cons*2*C*eps/(h^2) - d/(h^2) + f(x(J+2))/(h);
 A(1,3)= cons*C*eps/(h^2) + d/(2*h^2);
 A(end-1,end-1)= cons*C*eps/(h^2) + d/(2*h^2)- cons*eps*(1/(r+x(L+2*J))^alpha+1/(rr-x(L+2*J))^alpha)/(alpha) + f(x(L+2*J))/(h);
 A(end-1,end-2)= -cons*2*C*eps/(h^2) - d/(h^2) - f(x(L+2*J))/(h);
 A(end-1,end-3)= cons*C*eps/(h^2) + d/(2*h^2);
 A(end-1,end)=0;
end




%integral part when x>0
for j= -L+1:J-1
% coefficient of u_j
 b(j+L)= (- sum(1./abs(x(J+1-j:L+J)).^(1+alpha)) -sum(1./abs(x(L+J+2:2*J+L+1-j)).^(1+alpha))  +...
         .5/abs(x(J+1-j))^(1+alpha) +  .5/abs(x(2*J+L+1-j))^(1+alpha)  )*h;  
end
A=A+cons*eps.*(diag([b;0]));


% coefficient of u_(j+k) 
B=zeros(size(A));
B(1,:)=[0 1./abs(x(L+J+2:2*L+2*J-1)).^(1+alpha)   0.5/abs(x(2*L+2*J))^(1+alpha)].*h;
B(J+L-1,:)=[1./abs(x(3:L+J)).^(1+alpha)  0  0.5/abs(x(L+J+2))^(1+alpha)].*h;
for j=-L+2:J-2
B(L+j,:)=[1./abs(x(J+2-j:L+J)).^(1+alpha)  0   1./abs(x(L+J+2:2*J+L-j)).^(1+alpha)  .5/abs(x(2*J+L+1-j))^(1+alpha)].*h;
end
A = cons*eps.*B+A;


rhs=[-eps*cons/alpha./(rr-x(J+2:L+2*J)').^alpha;1];
U=A\rhs;

%  rescale:
X=-r:h:rr;


%plot(X,[0;U],'b--')


