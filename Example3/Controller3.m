clc;
clear all;
E=eye(2);
Im=[1;1];
I=1;
I1=1;
I2=1;
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
L1=5991;L2=9000;
b=0.97;
 k01=[-0.0006,-0.0070];
 k02=[-0.0006,-0.0070];

Qau1=(log(k1)+3*log(xi))/g1;
Qau2=(log(k2)+3*log(xi))/g2;

A11=A1';
A22=A2';
vecA11=A11(:);
vecA22=A22(:);
vecE=E(:);
vecI=I(:);
F11=F1';
F22=F2';
Xi=[xi,0;0,xi];
k011=k01';
k022=k02';
veck011=k011(:);
veck022=k022(:);

k0111=k01'*I;
k0222=k02'*I;
veck0111=k0111(:);
veck0222=k0222(:);

H11=[A1'+g1*eye(2),E,zeros(2,4),zeros(2,4),zeros(2,2)];
H14=[zeros(2,4),A2'+g2*eye(2),E,zeros(2,4),zeros(2,2)];

H21=[F11(1,1),F11(1,2),zeros(1,10),-1,zeros(1,1)];
H24=[zeros(1,4),F22(1,1),F22(1,2),zeros(1,6),-1,zeros(1,1)];

H41=[(E+Q1)',zeros(2,2),-k1*eye(2),zeros(2,6),zeros(2,2)];
H44=[-k2*eye(2),zeros(2,2),(E+Q2)',zeros(2,6),zeros(2,2)];

H51=[veck011*I'*B1',zeros(2,6),-kron(I1,E),zeros(2,2),zeros(2,2)];
H54=[zeros(2,4),veck022*I'*B2',zeros(2,2),zeros(2,2),-kron(I2,E),zeros(2,2)];

H61=[-veck0111*I'*B1',zeros(2,6),kron(I'*I1,E),zeros(2,2),zeros(2,2)];
H64=[zeros(2,4),-veck0222*I'*B2',zeros(2,2),zeros(2,2),kron(I'*I2,E),zeros(2,2)];

H71=[-I'*B1',zeros(1,10),zeros(1,1),1];
H74=[zeros(1,4),-I'*B2',zeros(1,6),zeros(1,1),1];

H81=[-vecA11*I'*((B1)'),zeros(4,6),-kron(B1*I1,E),zeros(4,2),zeros(4,2)];
H84=[zeros(4,4),-vecA22*I'*B2',zeros(4,4),-kron(B2*I1,E),zeros(4,2)];

H91=[zeros(2,2),-E,zeros(2,4),E,zeros(2,2),zeros(2,2)];
H93=[zeros(2,6),-E,zeros(2,2),E,zeros(2,2)];

H=[H11;H14;H21;H24;H41;H44;H51;H54;H61;H64;H71;H74;H81;H84;H91;H93];
c1=[-C1';-C2';-D1;-D2;zeros(8,1);b*Im;b*Im;zeros(2,1);L1*vecE;L2*vecE;zeros(4,1)];;
f=[0;0;0;0;0;0;0;0;0;0;0;0;1;-1];
[x,fval,exitflag]=linprog(f,H,c1,[],[]);

P10=x(1:2);
n1=x(3:4);
P20=x(5:6);
n2=x(7:8);
n11=x(9:10);
n21=x(11:12);
r=x(13);
c=x(14);
k10=[vpa((I1*n11'),8)/vpa((P10'*B1*I),8),vpa((I1*n21'),8)/vpa((P20'*B2*I),8)]
