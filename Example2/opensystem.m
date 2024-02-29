
clc;
clear all;
T=10;
N=10000;
stepsize=T/N;
a=ones(1,N/50);e=ones(1,N/100);c=ones(1,N/25);
b=2*ones(1,N/50);h=2*ones(1,N/50);f=2*ones(1,N/50); 
a1=ones(1,N/50);e1=ones(1,N/25);c1=ones(1,N/25);
b1=2*ones(1,N/25);h1=2*ones(1,N/25);f1=2*ones(1,N/50);

aa=0.1*ones(1,N/50);ee=0.2*ones(1,N/100);cc=0.3*ones(1,N/25);
bb=2.1*ones(1,N/50);hh=2.2*ones(1,N/50);ff=2.3*ones(1,N/50); 

aa1=0.1*ones(1,N/50);ee1=0.2*ones(1,N/25);cc1=0.3*ones(1,N/25);
bb1=2.1*ones(1,N/25);hh1=2.2*ones(1,N/25);ff1=2.3*ones(1,N/50); 
signal11=[aa ee cc bb hh ff aa1 ee1 cc1 bb hh ff aa ee cc bb1 bb1 hh1 ff1];
signal1=[signal11,signal11,2];
signal=[a e c b h f a1 e1 c1 b h f a e c b1 b1 h1 f1];
ssignal=[signal,signal,2];
ccontrolsignal=[a e c b h f a1 e1 c1 b h f a e c b1 b1 h1 f1 ];
controlsignal=[ccontrolsignal,ccontrolsignal,2];

st=0:stepsize:T;
X = linspace(0,T);
figure 
stairs(st,ssignal,'r','LineWidth',1),
xx=xlabel('Time(s)'),
set(xx,'Interpreter','Latex','FontName','Times New Roman','FontSize',12),
yy=ylabel('Switching signal $\sigma(t)$'),
set(yy,'Interpreter','Latex','FontName','Times New Roman','FontSize',12),
axis([0 10 0.5 2.5]),
set(gca,'FontSize',12)
%------------------------------------------------------------------
A1=[-1,3;2,2.5];
A2=[-2,2;3,-1.5];
B1=[0.1,0.2;0.2,0.3];
B2=[0.3,0.2;0.2,0.5];
C1=[0.1,0.2];
C2=[0.05,0.15];
D1=0.1;
D2=0.2;
F1=[0.2;0.3];
F2=[0.3;0.5];
Q1=[1.2,0;0,0.9];
Q2=[2,0;0,1.5];

k10 =[0 0;0 0];k11=[0 0;0 0];
k12=[0 0;0 0];k20=[0 0;0 0];
k21=[0 0;0 0];k22=[0 0;0 0];

x_0=[3;1]; % the initial state
x=[x_0,zeros(2,N)];
x0=[x_0,zeros(2,N)];
x1=[x_0,zeros(2,N)];
t=0:stepsize:T;
E=eye(2);
%本文扰动F1*exp(-0.2*k)
for k=1:1:N
        if signal1(k)==0.1
            u(:,k)=k10*x(:,k);
            x(:,1+k)= x(:,k)+(A1)*x(:,k)*stepsize+B1* u(:,k)*stepsize+F1*exp(-0.2*k)*stepsize;
             y(:,k)=C1*x(:,k)+D1*exp(-0.2*k);
        elseif signal1(k)==0.2
           u(:,k)=k11*x(:,k);
            x(:,1+k)= x(:,k)+(A1)*x(:,k)*stepsize+B1* u(:,k)*stepsize+F1*exp(-0.2*k)*stepsize;
         y(:,k)=C1*x(:,k)+D1*exp(-0.2*k);
        elseif signal1(k)==0.3
            u(:,k)=k12*x(:,k);
           x(:,1+k)= x(:,k)+(A1)*x(:,k)*stepsize+B1* u(:,k)*stepsize+F1*exp(-0.2*k)*stepsize;
             y(:,k)=C1*x(:,k)+D1*exp(-0.2*k);
           if signal1(k+1)==2.1
                 x(:,1+k)=(E+Q2)*x(:,1+k);
            end
       elseif signal1(k)==2.1
            u(:,k)=k20*x(:,k);
            x(:,1+k)= x(:,k)+(A2)*x(:,k)*stepsize+B2* u(:,k)*stepsize+F2*exp(-0.2*k)*stepsize;
        y(:,k)=C2*x(:,k)+D2*exp(-0.2*k);
        elseif signal1(k)==2.2
           u(:,k)=k21*x(:,k);
            x(:,1+k)= x(:,k)+(A2)*x(:,k)*stepsize+B2* u(:,k)*stepsize+F2*exp(-0.2*k)*stepsize;
        y(:,k)=C2*x(:,k)+D2*exp(-0.2*k);
        elseif signal1(k)==2.3
             u(:,k)=k22*x(:,k);
            x(:,1+k)= x(:,k)+(A2)*x(:,k)*stepsize+B2* u(:,k)*stepsize+F2*exp(-0.2*k)*stepsize;
              y(:,k)=C2*x(:,k)+D2*exp(-0.2*k);
            if signal1(k+1)==0.1
                 x(:,1+k)=(E+Q1)*x(:,1+k);
             end
        end
end

u(:,1+k)=k10*x(:,k);
y(:,1+k)=C2*x(:,1+k)+D2*exp(-0.2*(k+1));
%-----------------------------系统状态-------------------------------------------
figure;
plot(t,x(1,1:end),'r',t,x(2,1:end),'--r','LineWidth',1.5);
grid off
% axis([0 10 -0 3]),
LE=legend('$x_1$','$x_2$')
set(LE,'Interpreter','Latex','FontName','Times New Roman','FontSize',11),
 xx=xlabel('Time(s)'),
 set(xx,'Interpreter','Latex','FontName','Times New Roman','FontSize',12),
yy=ylabel('States $x(t)$'),
 set(yy,'Interpreter','Latex','FontName','Times New Roman','FontSize',12),
grid off
