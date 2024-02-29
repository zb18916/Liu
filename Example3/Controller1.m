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
b=0.50;
k01=[-0.0006,-0.0070];
k02=[-0.0006,-0.0070];

Qau1=(log(k1)+3*log(xi))/g1
Qau2=(log(k2)+3*log(xi))/g2

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

H11=[A1'+g1*eye(2),E,zeros(2,12),zeros(2,4),zeros(2,2)];%p10
H12=[zeros(2,2),E,A1'+g1*eye(2),zeros(2,10),zeros(2,4),zeros(2,2)];%p11
H13=[zeros(2,2),E,zeros(2,2),A1'+g1*eye(2),zeros(2,8),zeros(2,4),zeros(2,2)];%p12
H14=[zeros(2,8),A2'+g2*eye(2),E,zeros(2,4),zeros(2,4),zeros(2,2)];%p20
H15=[zeros(2,10),E,A2'+g2*eye(2),zeros(2,2),zeros(2,4),zeros(2,2)];%p21
H16=[zeros(2,10),E,zeros(2,2),A2'+g2*eye(2),zeros(2,4),zeros(2,2)];%p22

H21=[F11(1,1),F11(1,2),zeros(1,14),zeros(1,4),-1,zeros(1,1)];
H22=[zeros(1,4),F11(1,1),F11(1,2),zeros(1,10),zeros(1,4),-1,zeros(1,1)];
H23=[zeros(1,6),F11(1,1),F11(1,2),zeros(1,8),zeros(1,4),-1,zeros(1,1)];
H24=[zeros(1,8),F22(1,1),F22(1,2),zeros(1,6),zeros(1,4),-1,zeros(1,1)];
H25=[zeros(1,12),F22(1,1),F22(1,2),zeros(1,2),zeros(1,4),-1,zeros(1,1)];
H26=[zeros(1,14),F22(1,1),F22(1,2),zeros(1,4),-1,zeros(1,1)];

H31=[-xi*eye(2),zeros(2,2),eye(2),zeros(2,10),zeros(2,4),zeros(2,2)];%p11-xip10
H32=[zeros(2,4),-xi*eye(2),eye(2),zeros(2,8),zeros(2,4),zeros(2,2)];%p12-xip11
H33=[zeros(2,8),-xi*eye(2),zeros(2,2),eye(2),zeros(2,2),zeros(2,4),zeros(2,2)];%p21-xip20
H34=[zeros(2,12),-xi*eye(2),eye(2),zeros(2,4),zeros(2,2)];%p22-xip21

H41=[(E+Q1)',zeros(2,12),-k1*eye(2),zeros(2,4),zeros(2,2)];%p10-p22;
H44=[zeros(2,6),-k2*eye(2),(E+Q2)',zeros(2,6),zeros(2,4),zeros(2,2)];%p20-p12;

H51=[veck011*I'*B1',zeros(2,14),-kron(I1,E),zeros(2,2),zeros(2,2)];
H52=[zeros(2,4),veck011*I'*B1',zeros(2,10),-kron(I1,E),zeros(2,2),zeros(2,2)];
H53=[zeros(2,6),veck011*I'*B1',zeros(2,8),-kron(I1,E),zeros(2,2),zeros(2,2)];
H54=[zeros(2,8),veck022*I'*B2',zeros(2,6),zeros(2,2),-kron(I2,E),zeros(2,2)];
H55=[zeros(2,12),veck022*I'*B2',zeros(2,2),zeros(2,2),-kron(I2,E),zeros(2,2)];
H56=[zeros(2,14),veck022*I'*B2',zeros(2,2),-kron(I2,E),zeros(2,2)];

H61=[-veck0111*I'*B1',zeros(2,14),kron(I'*I1,E),zeros(2,2),zeros(2,2)];
H62=[zeros(2,4),-veck0111*I'*B1',zeros(2,10),kron(I'*I1,E),zeros(2,2),zeros(2,2)];
H63=[zeros(2,6),-veck0111*I'*B1',zeros(2,8),kron(I'*I1,E),zeros(2,2),zeros(2,2)];
H64=[zeros(2,8),-veck0222*I'*B2',zeros(2,6),zeros(2,2),kron(I'*I2,E),zeros(2,2)];
H65=[zeros(2,12),-veck0222*I'*B2',zeros(2,2),zeros(2,2),kron(I'*I2,E),zeros(2,2)];
H66=[zeros(2,14),-veck0222*I'*B2',zeros(2,2),kron(I'*I2,E),zeros(2,2)];

H71=[-I'*B1',zeros(1,18),zeros(1,1),1];
H72=[zeros(1,4),-I'*B1',zeros(1,14),zeros(1,1),1];
H73=[zeros(1,6),-I'*B1',zeros(1,12),zeros(1,1),1];
H74=[zeros(1,8),-I'*B2',zeros(1,10),zeros(1,1),1];
H75=[zeros(1,12),-I'*B2',zeros(1,6),zeros(1,1),1];
H76=[zeros(1,14),-I'*B2',zeros(1,4),zeros(1,1),1];

H81=[-vecA11*I'*((B1)'),zeros(4,14),-kron(B1*I1,E),zeros(4,2),zeros(4,2)];
H82=[zeros(4,4),-vecA11*I'*B1',zeros(4,10),-kron(B1*I1,E),zeros(4,2),zeros(4,2)];
H83=[zeros(4,6),-vecA11*I'*B1',zeros(4,8),-kron(B1*I1,E),zeros(4,2),zeros(4,2)];
H84=[zeros(4,8),-vecA22*I'*B2',zeros(4,8),-kron(B2*I1,E),zeros(4,2)];
H85=[zeros(4,12),-vecA22*I'*B2',zeros(4,4),-kron(B2*I1,E),zeros(4,2)];
H86=[zeros(4,14),-vecA22*I'*B2',zeros(4,2),-kron(B2*I1,E),zeros(4,2)];

H91=[zeros(2,2),-E,zeros(2,12),E,zeros(2,2),zeros(2,2)];
H93=[zeros(2,10),-E,zeros(2,6),E,zeros(2,2)];

H=[H11;H12;H13;H14;H15;H16;H21;H22;H23;H24;H25;H26;H31;H32;H33;H34;H41;H44
    H51; H52; H53;H54;H55;H56;H61;H62;H63;H64;H65;H66;H71;H72;H73;H74;H75;H76;
    H81;H82;H83;H84;H85;H86;H91;H93;];
c1=[-C1';-C1';-C1';-C2';-C2';-C2';-D1;-D1;-D1;-D2;-D2;-D2;zeros(12,1);zeros(12,1);
    b*Im; b*Im; b*Im;b*Im; b*Im; b*Im;zeros(6,1);
    L1*vecE;L1*vecE;L1*vecE;L2*vecE;L2*vecE;L2*vecE;zeros(4,1)];
f=[zeros(20,1);1;-1];
[x,fval,exitflag]=linprog(f,H,c1,[],[]);
P10=x(1:2);n1=x(3:4);P11=x(5:6);P12=x(7:8);P20=x(9:10);n2=x(11:12);P21=x(13:14);
P22=x(15:16);n11=x(17:18);n21=x(19:20);
r=x(21);
c=x(22);
k10=[vpa((I1*n11'),8)/vpa((P10'*B1*I),8),vpa((I1*n11'),8)/vpa((P11'*B1*I),8),vpa((I1*n11'),8)/vpa((P12'*B1*I),8);
vpa((I1*n21'),8)/vpa((P20'*B2*I),8),vpa((I1*n21'),8)/vpa((P21'*B2*I),8),vpa((I1*n21'),8)/vpa((P22'*B2*I),8)]
