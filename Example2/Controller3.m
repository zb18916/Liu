clc;
clear all;
E=eye(2);
I=[1;1];
I0=[1,1];
I1=[1;0];
I2=[0;1];

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
G=2;
k1=2.4;k2=3.2;
xi=0.99;
g1=1.5;g2=3.7; 
L1=0.6;L2=1.5;
b=1.;a=4.8;r=1;
k01=[-4.00,-10.00;-4.00,-10.00];k02=[-4.3,-4;-4.3,-4];

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

H11=[A1'+g1*eye(2),E,zeros(2,4),zeros(2,8)];
H14=[zeros(2,4),A2'+g2*eye(2),E,zeros(2,8)];

H21=[F11(1,1),F11(1,2),zeros(1,6),zeros(1,8)];
H24=[zeros(1,4),F22(1,1),F22(1,2),zeros(1,10)];

H41=[(E+Q1)',zeros(2,2),-k1*eye(2),zeros(2,10)];
H44=[-k2*eye(2),zeros(2,2),(E+Q2)',zeros(2,10)];

H51=[veck011*I'*B1',zeros(4,6),-kron(I1,E),-kron(I2,E),zeros(4,4)];
H54=[zeros(4,4),veck022*I'*B2',zeros(4,6),-kron(I1,E),-kron(I2,E)];

H61=[-veck0111*I'*B1',zeros(2,6),kron(I'*I1,E),kron(I'*I2,E),zeros(2,4)];
H64=[zeros(2,4),-veck0222*I'*B2',zeros(2,6),kron(I'*I1,E),kron(I'*I2,E)];

H71=[-I'*B1',zeros(1,14)];
H74=[zeros(1,4),-I'*B2',zeros(1,10)];

H81=[-vecA11*I'*((B1)'),zeros(4,6),-kron(((B1)')*I1,E),-kron(((B1)')*I2,E),zeros(4,4)];
H84=[zeros(4,4),-vecA22*I'*B2',zeros(4,6),-kron(B2'*I1,E),-kron(B2'*I2,E)];

H91=[zeros(2,2),-E,zeros(2,4),E,zeros(2,6)];
H92=[zeros(2,2),-E,zeros(2,6),E,zeros(2,4)];
H93=[zeros(2,6),-E,zeros(2,4),E,zeros(2,2)];
H94=[zeros(2,6),-E,zeros(2,6),E];

H=[H11;H14;H21;H24;H41;H44
    H51; H54;H61;H64;H71;H74;H81;H84;H91;H92;H93;H94];

c1=[-C1';-C2';r-D1;r-D2;zeros(12,1);b*vecI;b*vecI;-b/a;-b/a;L1*vecE;L2*vecE;zeros(8,1)];
f=zeros(16,1);
[x,fval,exitflag]=linprog(f,H,c1,[],[]);
P10=x(1:2);
n1=x(3:4);
P20=x(5:6);
n2=x(7:8);
n11=x(9:10);
n12=x(11:12);
n21=x(13:14);
n22=x(15:16);
k10=[(I1*n11'+I2*n12')/(P10'*B1*I);(I1*n21'+I2*n22')/(P20'*B2*I)]
