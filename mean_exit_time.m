% meant exit time from a bounded domain 
% levy motion


function U=mean_exit_time(Alp)


alpha=Alp;
eps= 1;
d=0;
f=@(x) -x;
r=2; 


J=100;
h=1/J;
x=-2:h:2;
C= -zeta(alpha-1)*h^(2-alpha);          % correction term u''(x)
cons=alpha*gamma((1+alpha)/2)/(2^(1-alpha)*sqrt(pi)*gamma(1-alpha/2)); 

b=zeros(2*J-1,1);
a=zeros(2*J-1,1);
c=zeros(2*J-1,1); 

%nonintegral part
for i=J+2:3*J
    b(i-J-1)= -cons*2*C*eps/(h^2*r^(alpha)) - d/(h^2*r^2) - cons*eps*(1/(1+x(i))^alpha+1/(1-x(i))^alpha)/(alpha*r^(alpha)); % coefficient of u_j
    a(i-J-1)= cons*C*eps/(h^2*r^(alpha)) + d/(2*h^2*r^2) - f(x(i))/(2*h*r);  % coefficient of u_(j-1)
    c(i-J-1)= cons*C*eps/(h^2*r^(alpha)) + d/(2*h^2*r^2) + f(x(i))/(2*h*r);  % coefficient of u_(j+1)   
end
A=diag(b)+diag(a(2:end),-1)+diag(c(1:end-1),1);
%Jt=2*J-1;
%A=spdiags([[a(2:end); 0]  b  [0; c(1:end-1)] ],-1:1,Jt,Jt);

       
       
if alpha<1
A(1,1)= + d/(2*h^2*r^2)- cons*eps*(1/(1+x(J+2))^alpha+1/(1-x(J+2))^alpha)/(alpha*r^(alpha)) - f(x(J+2))/(h*r);
A(1,2)= - d/(h^2*r^2) + f(x(J+2))/(h*r);
A(1,3)= + d/(2*h^2*r^2);
A(end,end)= + d/(2*h^2*r^2)- cons*eps*(1/(1+x(3*J))^alpha+1/(1-x(3*J))^alpha)/(alpha*r^(alpha)) + f(x(3*J))/(h*r);
A(end,end-1)= - d/(h^2*r^2) - f(x(3*J))/(h*r);
A(end,end-2)= + d/(2*h^2*r^2);
end


%integral part when x>0
for j= 0:J-1
% coefficient of u_j
 b(j+J)= (- sum(1./abs(x(J+1-j:2*J)).^(1+alpha)) -sum(1./abs(x(2*J+2:3*J+1-j)).^(1+alpha))  +...
         .5/abs(x(J+1-j))^(1+alpha) +  .5/abs(x(3*J+1-j))^(1+alpha)  )*h/r^alpha;  
     
% coefficient of u_(j-1)
a(j+J)= (sum(x(J+1+j:2*J).*0.5./abs(x(J+1+j:2*J)).^(1+alpha))+...                
         sum(x(2*J+2:3*J+1-j).*0.5./abs(x(2*J+2:3*J+1-j)).^(1+alpha)) -...
         .5 * (x(J+1+j)/2/abs(x(J+1+j))^(1+alpha)+x(3*J+1-j)/2/abs(x(3*J+1-j))^(1+alpha) ) )/r^alpha;
% coefficient of u_(j+1)
c(j+J)= -a(j+J);
end

%integral part when x<0     
for j=-J+1:-1
% coefficient of u_j
b(j+J)= (- sum(1./abs(x(J+1-j:2*J)).^(1+alpha)) -sum(1./abs(x(2*J+2:3*J+1-j)).^(1+alpha)) +...
        .5/abs(x(J+1-j))^(1+alpha) +  .5/abs(x(3*J+1-j))^(1+alpha)     )*h/r^alpha;   

% coefficient of u_(j-1)
a(j+J)= (sum(x(J+1-j:2*J)*0.5./abs(x(J+1-j:2*J)).^(1+alpha))+...
          sum(x(2*J+2:3*J+1+j)*0.5./abs(x(2*J+2:3*J+1+j)).^(1+alpha)) -...
        .5 * (x(J+1-j)/2/abs(x(J+1-j))^(1+alpha)+x(3*J+1+j)/2/abs(x(3*J+1+j))^(1+alpha))  )/r^alpha;

% coefficient of u_(j+1)
c(j+J)= -a(j+J);
end
 A=A+cons*eps.*(diag(b)+diag(a(2:end),-1)+diag(c(1:end-1),1));

% %coefficient of u_(j+k) 
% B=zeros(size(A));
% for j=-J+1:J-1
% B(J+j,:)=[1./abs(x(J+2-j:2*J)).^(1+alpha)  0  1./abs(x(2*J+2:3*J-j)).^(1+alpha)].*h./r^alpha;
% end


B=zeros(size(A));
B(1,:)=[0 1./abs(x(2*J+2:4*J-1)).^(1+alpha)].*h./r^alpha;
B(2*J-1,:)=[1./abs(x(3:2*J)).^(1+alpha) 0].*h./r^alpha;
for j=-J+2:J-2
B(J+j,:)=[1./abs(x(J+2-j:2*J)).^(1+alpha)  0  1./abs(x(2*J+2:3*J-j)).^(1+alpha)].*h./r^alpha;
end



A = cons*eps.*B+A;

U=A\(-1*ones(length(A),1));
U=[0; U; 0];
%  rescale:
%X=-r:(h*r):r;







