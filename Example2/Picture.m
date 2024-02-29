
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
%---------------------------Controller1--------------------------
k10=[-3.8464   -9.8010;-3.8464   -9.8010];k11=[-3.8853   -9.9000;-3.8853   -9.9000];
k12=[-3.9245  -10.0000;-3.9245  -10.0000];k20=[-4.2004   -3.9097;-4.2004   -3.9097];
k21=[-4.2429   -3.9492;-4.2429   -3.9492];k22=[ -4.2857  -3.9891;-4.2857   -3.9891];
%---------------------------Controller2--------------------------
k100=[ -3.8179   -9.8010 ; -3.8179   -9.8010 ];k111=[-3.8565  -9.9000   ;-3.8565   -9.9000 ];
k122=[-3.8955  -10.0000;-3.8955  -10.0000];
k200=[-4.2004   -3.9204;-4.2004   -3.9204];
k211=[-4.2429   -3.9600;-4.2429   -3.9600];k222=[ -4.2857   -4.0000; -4.2857   -4.0000];
%---------------------------Controller3--------------------------
 k01=[ -4.0000   -9.3714
   -4.0000   -9.3714];
 k02=[ -4.2536   -3.9711
   -4.2536   -3.9711];

x_0=[3;1]; % the initial state
x=[x_0,zeros(2,N)];
x0=[x_0,zeros(2,N)];
x1=[x_0,zeros(2,N)];
t=0:stepsize:T;
E=eye(2);
%本文扰动F1*exp(-0.2*k)
%---------------------------Controller1--------------------------
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
%---------------------------Controller2--------------------------
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
%---------------------------Controller3--------------------------
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

%-----------------------------系统状态-------------------------------------------
figure;
plot(t,x(1,1:end),'r',t,x(2,1:end),'--r','LineWidth',1.5);
hold on
plot(t,x0(1,1:end),'-k',t,x0(2,1:end),'--k','LineWidth',1.5);
hold on
plot(t,x1(1,1:end),'-b',t,x1(2,1:end),'--b','LineWidth',1.5);

grid off
axis([0 10 -0 3]),
LE=legend('$x_1$ under Controller 1','$x_2$ under  Controller 1','$x_1$  under Controller 2','$x_2$ under Controller 2','$x_1$ under Controller 3','$x_2$ under Controller 3')
set(LE,'Interpreter','Latex','FontName','Times New Roman','FontSize',14),
xx=xlabel('Time(s)'),
set(xx,'Interpreter','Latex','FontName','Times New Roman','FontSize',12),
yy=ylabel('States $x(t)$'),
set(yy,'Interpreter','Latex','FontName','Times New Roman','FontSize',12),
grid off
%-----------------------------系统输出-------------------------------------------
figure;
plot(t,y(1,1:end),'r','LineWidth',1.5);
hold on;
plot(t,y0(1,1:end),'--k','LineWidth',1.5);
hold on;grid on;
plot(t,y1(1,1:end),'-.b','LineWidth',1.5);

grid off
axis([0 10 -0 0.7]),
LE=legend('$z(t)$ under Controller 1','$z(t)$ under Controller 2','$z(t)$ under Controller 3')
set(LE,'Interpreter','Latex','FontName','Times New Roman','FontSize',14),
xx=xlabel('Time(s)'),
set(xx,'Interpreter','Latex','FontName','Times New Roman','FontSize',12),
yy=ylabel('Output z(t)'),
set(yy,'Interpreter','Latex','FontName','Times New Roman','FontSize',12),

%--------------------控制输入----------------------------------------------------
figure;
plot(t,u(1,1:end),'r','LineWidth',1.5);
hold on;
plot(t,u0(1,1:end),'--k','LineWidth',1.5);
hold on;
plot(t,u1(1,1:end),'-.b','LineWidth',1.5);
% title('控制输入u')
% legend('本文u','无无扰切换u','恒定增益u')
grid off
axis([0 10 -25 0]),
LE=legend('$u(t)$ under Controller 1','$u(t)$ under Controller 2','$u(t)$ under Controller 3')
set(LE,'Interpreter','Latex','FontName','Times New Roman','FontSize',12),
xx=xlabel('Time(s)'),
set(xx,'Interpreter','Latex','FontName','Times New Roman','FontSize',12),
yy=ylabel('Control Input u(t)'),
set(yy,'Interpreter','Latex','FontName','Times New Roman','FontSize',12),
% %-------------------局部放大图-1--------------------------------
axes('Position',[0.35,0.25,0.25,0.2]); % 生成子图  % 生成子图，第一个位置控制左右移动，第二个位置控制上下移动，第三个位置控制 宽窄
plot(t,u(1,1:end),'r','LineWidth',1.5);
hold on;grid on;
plot(t,u0(1,1:end),'--k','LineWidth',1.5);
hold on;
plot(t,u1(1,1:end),'-.b','LineWidth',1.5);
grid off
xlim([0.69,0.74]);
ylim([-10.8,-9]);
X11=[0.36,0.2];
Y11=[0.48,0.62];
annotation('arrow',X11,Y11)
% %-------------------局部放大图-2--------------------------------
axes('Position',[0.38,0.55,0.25,0.2]); % 生成子图  % 生成子图，第一个位置控制左右移动，第二个位置控制上下移动，第三个位置控制 宽窄
plot(t,u(1,1:end),'r','LineWidth',1.5);
hold on;grid on;
plot(t,u0(1,1:end),'--k','LineWidth',1.5);
hold on;
plot(t,u1(1,1:end),'-.b','LineWidth',1.5);
%plot(t,u1(1,1:end),);
grid off
xlim([1.295,1.33]);
ylim([-2.27,-2.1]);
X11=[0.35,0.24];
Y11=[0.65,0.86];
annotation('arrow',X11,Y11)
%----------------------------