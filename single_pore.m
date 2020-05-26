clear all
clc
close all
mu1=3.3e-3; mu2=0; L=1e-6; sig=1.6; small=1e-2; r=5.5e-6;
x1=small; t1=small; t2=10; p1=0; p2=0; A=pi*r^2; pc=2*pi*r*(sig/A);
dp=p1-p2+pc;
%dp=p1-p2;
dt=0.05; n=1000; x=zeros(1,n+1); xe=zeros(1,n+1); t=zeros(1,n+1);
x(1)=x1;
t(1)=t1;
xe(1)=x1;
q=0;
v=0;
for i=1:n
a=8*(mu1*x(i) + mu2*(L-x(i))) ;  b=pi*r^4;
q=dp*b/a; v=q/A;
x(i+1)=x(i)+v*dt;
t(i+1)=t(i)+dt;
xe(i+1)=sqrt((r*sig*t(i+1))/(2*mu1));
end
plot(t,x,t,xe)
fprintf('Infiltration at t=%f should be %f',t(n+1),x(n+1))
ma=[t;x;xe];
mb=ma';
csvwrite('validation.csv',mb);