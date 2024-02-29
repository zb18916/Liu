clc;
clear all;
T=10;
N=10000;
stepsize=T/N;
% %------------------------------------考虑无扰切换（不同的驻留时间）(颠簸变大）-------------------------------
a=ones(1,N/50);e=ones(1,N/100);c=ones(1,N/25);%0.2  0.1  0.4  
b=2*ones(1,N/50);h=2*ones(1,N/50);f=2*ones(1,N/50); %0.2 0.2 0.2 

a1=ones(1,N/50);e1=ones(1,N/25);c1=ones(1,N/25);%0.2 0.4  0.4  子系统一大于1.6
b1=2*ones(1,N/25);h1=2*ones(1,N/25);f1=2*ones(1,N/50); %0.4 0.4 0.2  子系统二大于2.5

aa=0.1*ones(1,N/50);ee=0.2*ones(1,N/100);cc=0.3*ones(1,N/25);%0.2  0.1  0.4  子系统一大于0.1066
bb=2.1*ones(1,N/50);hh=2.2*ones(1,N/50);ff=2.3*ones(1,N/50); %0.1 0.1 0.1   子系统二大于0.1580

aa1=0.1*ones(1,N/50);ee1=0.2*ones(1,N/25);cc1=0.3*ones(1,N/25);%0.4  0.2  0.4  子系统一大于0.1066
bb1=2.1*ones(1,N/25);hh1=2.2*ones(1,N/25);ff1=2.3*ones(1,N/50); %0.4 0.4 0.2   子系统二大于0.1580
signal11=[aa ee cc bb hh ff aa1 ee1 cc1 bb hh ff aa ee cc bb1 bb1 hh1 ff1];%总共N个
signal1=[signal11,signal11,2];
signal=[a e c b h f a1 e1 c1 b h f a e c b1 b1 h1 f1];%总共N个
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
%----------------------本文状态--------------------------------------------
A1=[-3.7401,1.4001;0.4752,-4.5586];
A2=[-4.1804,1.5321;0.4244,-4.9290];
B1=[231.5508;657.384];
B2=[231.8138;674.5305];
C1=[1,0];
C2=[1,0];
D1=0.5;
D2=1.5;
F1=[3;2];
F2=[1;1.8];
Q1=[0.012,0;0,0.09];
Q2=[0.008,0;0,0.011];
G=2;
k1=1.4;k2=3;
xi=0.99;
g1=3.2;g2=3.95;
L1=4991;L2=9000;
b=0.97;

%-----------------Controller1--------------------
k10 =[ -0.00058806, -0.005926];
k11=[-0.00059399, -0.0059861];
k12=[-0.0006000, -0.0060466];
k20=[-0.00058806, -0.0063187];
k21=[-0.00059399, -0.0063827];
k22=[-0.00059999, -0.006447];
%-------------------------Controller2-----------------------------
k100 =[ -0.00070848, -0.0059262];
k111=[-0.000715636, -0.0059861];
k122=[-0.00072286, -0.0060466];
k200=[-0.00061665, -0.0060303];
k211=[-0.00062288, -0.0060912];
k222=[-0.00062917, -0.0061528];
%------------------------Controller3---------------
 k01=[-0.0005999   -0.005659];
 k02=[-0.0006000  -0.005952];

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
%--Controller2
for k=1:1:N
        if signal1(k)==0.1
            u0(:,k)=k100*x0(:,k);
            x0(:,1+k)= x0(:,k)+(A1)*x0(:,k)*stepsize+B1* u0(:,k)*stepsize+F1*exp(-0.2*k)*stepsize;
             y0(:,k)=C1*x0(:,k)+D1*exp(-0.2*k);
        elseif signal1(k)==0.2
           u0(:,k)=k111*x0(:,k);
            x0(:,1+k)= x0(:,k)+(A1)*x0(:,k)*stepsize+B1* u0(:,k)*stepsize+F1*exp(-0.2*k)*stepsize;
         y0(:,k)=C1*x0(:,k)+D1*exp(-0.2*k);
        elseif signal1(k)==0.3
            u0(:,k)=k122*x0(:,k);
           x0(:,1+k)= x0(:,k)+(A1)*x0(:,k)*stepsize+B1* u0(:,k)*stepsize+F1*exp(-0.2*k)*stepsize;
             y0(:,k)=C1*x0(:,k)+D1*exp(-0.2*k);
           if signal1(k+1)==2.1
                 x0(:,1+k)=(E+Q2)*x0(:,1+k);
            end
       elseif signal1(k)==2.1
            u0(:,k)=k200*x0(:,k);
           
            x0(:,1+k)= x0(:,k)+(A2)*x0(:,k)*stepsize+B2* u0(:,k)*stepsize+F2*exp(-0.2*k)*stepsize;
        y0(:,k)=C2*x0(:,k)+D2*exp(-0.2*k);
        elseif signal1(k)==2.2
           u0(:,k)=k211*x0(:,k);
            x0(:,1+k)= x0(:,k)+(A2)*x0(:,k)*stepsize+B2* u0(:,k)*stepsize+F2*exp(-0.2*k)*stepsize;
        y0(:,k)=C2*x0(:,k)+D2*exp(-0.2*k);
        elseif signal1(k)==2.3
             u0(:,k)=k222*x0(:,k);
            x0(:,1+k)= x0(:,k)+(A2)*x0(:,k)*stepsize+B2* u0(:,k)*stepsize+F2*exp(-0.2*k)*stepsize;
              y0(:,k)=C2*x0(:,k)+D2*exp(-0.2*k);
            if signal1(k+1)==0.1
                 x0(:,1+k)=(E+Q1)*x0(:,1+k);
             end
        end
end
%--Controller3
for k=1:1:N
        if signal1(k)==0.1
            u1(:,k)=k01*x1(:,k);
            x1(:,1+k)= x1(:,k)+(A1)*x1(:,k)*stepsize+B1* u1(:,k)*stepsize+F1*exp(-0.2*k)*stepsize;
             y1(:,k)=C1*x1(:,k)+D1*exp(-0.2*k);
        elseif signal1(k)==0.2
           u1(:,k)=k01*x1(:,k);
            x1(:,1+k)=x1(:,k)+(A1)*x1(:,k)*stepsize+B1*u1(:,k)*stepsize+F1*exp(-0.2*k)*stepsize;
         y1(:,k)=C1*x1(:,k)+D1*exp(-0.2*k);
        elseif signal1(k)==0.3
            u1(:,k)=k01*x1(:,k);
           x1(:,1+k)= x1(:,k)+(A1)*x1(:,k)*stepsize+B1* u1(:,k)*stepsize+F1*exp(-0.2*k)*stepsize;
             y1(:,k)=C1*x1(:,k)+D1*exp(-0.2*k);
           if signal1(k+1)==2.1
                 x1(:,1+k)=(E+Q2)*x1(:,1+k);
            end
       elseif signal1(k)==2.1
            u1(:,k)=k02*x1(:,k);
           
            x1(:,1+k)= x1(:,k)+(A2)*x1(:,k)*stepsize+B2* u1(:,k)*stepsize+F2*exp(-0.2*k)*stepsize;
        y1(:,k)=C2*x1(:,k)+D2*exp(-0.2*k);
        elseif signal1(k)==2.2
           u1(:,k)=k02*x1(:,k);
            x1(:,1+k)= x1(:,k)+(A2)*x1(:,k)*stepsize+B2* u1(:,k)*stepsize+F2*exp(-0.2*k)*stepsize;
        y1(:,k)=C2*x1(:,k)+D2*exp(-0.2*k);
        elseif signal1(k)==2.3
             u1(:,k)=k02*x1(:,k);
            x1(:,1+k)= x1(:,k)+(A2)*x1(:,k)*stepsize+B2* u1(:,k)*stepsize+F2*exp(-0.2*k)*stepsize;
              y1(:,k)=C2*x1(:,k)+D2*exp(-0.2*k);
            if signal1(k+1)==0.1
                 x1(:,1+k)=(E+Q1)*x1(:,1+k);
             end
        end
end
u(:,1+k)=k10*x(:,k);
y(:,1+k)=C2*x(:,1+k)+D2*exp(-0.2*(k+1));
u0(:,1+k)=k100*x0(:,k);
y0(:,1+k)=C2*x0(:,1+k)+D2*exp(-0.2*(k+1));
u1(:,1+k)=k01*x1(:,k);
y1(:,1+k)=C2*x1(:,1+k)+D2*exp(-0.2*(k+1));
figure;
plot(t,u(1,1:end),'r','LineWidth',1.5);
hold on;
plot(t,u0(1,1:end),'--k','LineWidth',1.5);
hold on;
plot(t,u1(1,1:end),'-.b','LineWidth',1.5);
% title('控制输入u')
% legend('本文u','无无扰切换u','恒定增益u')
grid off
axis([0 10 -0.008 0]),
LE=legend('$W_f^*(t)$ under Controller 1','$W_f^*(t)$ under Controller 2','$W_f^*(t)$ under Controller 3')
set(LE,'Interpreter','Latex','FontName','Times New Roman','FontSize',11),
xx=xlabel('Time(s)'),
set(xx,'Interpreter','Latex','FontName','Times New Roman','FontSize',12),
yy=ylabel('The fuel flow $W_f^*(t)(pps/s)$'),
set(yy,'Interpreter','Latex','FontName','Times New Roman','FontSize',12),
% %-------------------局部放大图-1--------------------------------
axes('Position',[0.3,0.3,0.25,0.2]); % 生成子图  % 生成子图，第一个位置控制左右移动，第二个位置控制上下移动，第三个位置控制 宽窄
plot(t,u(1,1:end),'r','LineWidth',1.5);
hold on;grid on;
plot(t,u0(1,1:end),'--k','LineWidth',1.5);
hold on;
plot(t,u1(1,1:end),'-.b','LineWidth',1.5);
grid off
xlim([0.69,0.71]);
ylim([-0.00017,-0.00012]);
X11=[0.3,0.2];
Y11=[0.54,0.90];
annotation('arrow',X11,Y11)
% % % %-------------------局部放大图-2--------------------------------
axes('Position',[0.38,0.55,0.25,0.2]); % 生成子图  % 生成子图，第一个位置控制左右移动，第二个位置控制上下移动，第三个位置控制 宽窄
plot(t,u(1,1:end),'r','LineWidth',1.5);
hold on;grid on;
plot(t,u0(1,1:end),'--k','LineWidth',1.5);
hold on;
plot(t,u1(1,1:end),'-.b','LineWidth',1.5);
%plot(t,u1(1,1:end),);
grid off
xlim([1.295,1.305]);
ylim([-0.000011,-0.0000085]);
X11=[0.35,0.24];
Y11=[0.75,0.92];
annotation('arrow',X11,Y11)
% %----------------------------