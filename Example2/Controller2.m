clc;
clear  all;
E=eye(2);
I=[1;1];
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
k1=2.5;k2=3.2;
xi=0.99;r=1;
g1=1.5;g2=3.7;
L1=0.6;L2=1.5;

A11=A1';
A22=A2';
vecA11=A11(:);
vecA22=A22(:);
vecE=E(:);
F11=F1';
F22=F2';

Qau1=(log(k1)+log(xi))/g1;
Qau2=(log(k2)+log(xi))/g2;

H11=[A1'+g1*eye(2),E,zeros(2,12),zeros(2,8)];%p10
H12=[zeros(2,2),E,A1'+g1*eye(2),zeros(2,10),zeros(2,8)];%p11
H13=[zeros(2,2),E,zeros(2,2),A1'+g1*eye(2),zeros(2,8),zeros(2,8)];%p12
H14=[zeros(2,8),A2'+g2*eye(2),E,zeros(2,4),zeros(2,8)];%p20
H15=[zeros(2,10),E,A2'+g2*eye(2),zeros(2,2),zeros(2,8)];%p21
H16=[zeros(2,10),E,zeros(2,2),A2'+g2*eye(2),zeros(2,8)];%p22

H21=[F11(1,1),F11(1,2),zeros(1,14),zeros(1,8)];
H22=[zeros(1,4),F11(1,1),F11(1,2),zeros(1,10),zeros(1,8)];
H23=[zeros(1,6),F11(1,1),F11(1,2),zeros(1,8),zeros(1,8)];
H24=[zeros(1,8),F22(1,1),F22(1,2),zeros(1,6),zeros(1,8)];
H25=[zeros(1,12),F22(1,1),F22(1,2),zeros(1,2),zeros(1,8)];
H26=[zeros(1,14),F22(1,1),F22(1,2),zeros(1,8)];

H31=[-xi*eye(2),zeros(2,2),eye(2),zeros(2,10),zeros(2,8)];%p11-xip10
H32=[zeros(2,4),-xi*eye(2),eye(2),zeros(2,8),zeros(2,8)];%p12-xip11
H33=[zeros(2,8),-xi*eye(2),zeros(2,2),eye(2),zeros(2,2),zeros(2,8)];%p21-xip20
H34=[zeros(2,12),-xi*eye(2),eye(2),zeros(2,8)];%p22-xip21

H41=[(E+Q1)',zeros(2,12),-k1*eye(2),zeros(2,8)];%p10-p22;
H44=[zeros(2,6),-k2*eye(2),(E+Q2)',zeros(2,6),zeros(2,8)];%p20-p12;

H81=[-vecA11*I'*B1',zeros(4,14),-kron(B1*I1,E),-kron(B1*I2,E),zeros(4,4)];
H82=[zeros(4,4),-vecA11*I'*B1',zeros(4,10),-kron(B1*I1,E),-kron(B1*I2,E),zeros(4,4)];
H83=[zeros(4,6),-vecA11*I'*B1',zeros(4,8),-kron(B1*I1,E),-kron(B1*I2,E),zeros(4,4)];
H84=[zeros(4,8),-vecA22*I'*B2',zeros(4,10),-kron(B2*I1,E),-kron(B2*I2,E)];
H85=[zeros(4,12),-vecA22*I'*B2',zeros(4,6),-kron(B2*I1,E),-kron(B2*I2,E)];
H86=[zeros(4,14),-vecA22*I'*B2',zeros(4,4),-kron(B2*I1,E),-kron(B2*I2,E)];

H91=[zeros(2,2),-E,zeros(2,12),E,zeros(2,6)];
H92=[zeros(2,2),-E,zeros(2,14),E,zeros(2,4)];
H93=[zeros(2,10),-E,zeros(2,8),E,zeros(2,2)];
H94=[zeros(2,10),-E,zeros(2,10),E];

H=[H11;H12;H13;H14;H15;H16;H21;H22;H23;H24;H25;H26;H31;H32;H33;H34;H41;H44;
    H81;H82;H83;H84;H85;H86;H91;H92;H93;H94];
c1=[-C1';-C1';-C1';-C2';-C2';-C2';r-D1;r-D1;r-D1;r-D2;r-D2;r-D2;zeros(12,1);
   L1*vecE;L1*vecE;L1*vecE;L2*vecE;L2*vecE;L2*vecE;zeros(8,1)];
f=zeros(24,1);
[x,fval,exitflag]=linprog(f,H,c1,[],[]);
P10=x(1:2);
n1=x(3:4);
P11=x(5:6);
P12=x(7:8);
P20=x(9:10);
n2=x(11:12);
P21=x(13:14);
P22=x(15:16);
n11=x(17:18);
n12=x(19:20);
n21=x(21:22);
n22=x(23:24);
PP=[P10,P11,P12,P20,P21,P22];
k10=[(I1*n11'+I2*n12')/(P10'*B1*I),(I1*n11'+I2*n12')/(P11'*B1*I),(I1*n11'+I2*n12')/(P12'*B1*I);
(I1*n21'+I2*n22')/(P20'*B2*I),(I1*n21'+I2*n22')/(P21'*B2*I),(I1*n21'+I2*n22')/(P22'*B2*I)]