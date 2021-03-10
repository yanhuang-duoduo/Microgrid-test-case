clc;
clear all;
close all;
format long

A1=[1,1,1;
   1*exp(-2/3*pi*j),1*exp(2/3*pi*j),1;
   1*exp(2/3*pi*j),1*exp(-2/3*pi*j),1];

A2=1/3*[1,1,1;
   1,1*exp(2/3*pi*j),1*exp(-2/3*pi*j);
   1,1*exp(-2/3*pi*j),1*exp(2/3*pi*j)];
% branch 1
Z12s=0.003608000000000+0.003189000000000*j;
Z1_2 = [Z12s;Z12s;3*Z12s];
Z1_2s = A1*(diag([Z1_2(1,1),Z1_2(2,1),Z1_2(3,1)],0))*inv(A1);
Y1_2 = inv(Z1_2s);  % Three phase branch admittance
gsh12=[0,0,0;
       0,0,0;
       0,0,0];
bsh12=[0,0,0;
       0,0,0;
       0,0,0];
Ysh12=gsh12+bsh12*j;

Y12_11=Y1_2+Ysh12; 
Y12_12=-Y1_2;
Y12_21=-Y1_2;
Y12_22=Y1_2+Ysh12;

% branch 2
Z13s=0.005683750000000+0.004634062500000*j;
Z1_3 = [Z13s;Z13s;3*Z13s];
Z1_3s = A1*(diag([Z1_3(1,1),Z1_3(2,1),Z1_3(3,1)],0))*inv(A1);
Y1_3 = inv(Z1_3s);
gsh13=[0,0,0;
       0,0,0;
       0,0,0];
bsh13=[0,0,0;
       0,0,0;
       0,0,0];
Ysh13=gsh13+bsh13*j;
Y13_11=Y1_3+Ysh13;
Y13_12=-Y1_3;
Y13_21=-Y1_3;
Y13_22=Y1_3+Ysh13;
% branch 3
Z24s=0.052912875000000+0.039472500000000*j;
Z2_4 = [Z24s;Z24s;3*Z24s];
Z2_4s = A1*(diag([Z2_4(1,1),Z2_4(2,1),Z2_4(3,1)],0))*inv(A1);
Y2_4 = inv(Z2_4s);
gsh24=[0,0,0;
       0,0,0;
       0,0,0];
bsh24=[0,0,0;
       0,0,0;
       0,0,0];
Ysh24=gsh24+bsh24*j;
Y24_11=Y2_4+Ysh24;
Y24_12=-Y2_4;
Y24_21=-Y2_4;
Y24_22=Y2_4+Ysh24;

% branch 4
Z35s=0.052912875000000+0.039472500000000*j;
Z3_5 = [Z35s;Z35s;3*Z35s];
Z3_5s = A1*(diag([Z3_5(1,1),Z3_5(2,1),Z3_5(3,1)],0))*inv(A1);
Y3_5 = inv(Z3_5s);
gsh35=[0,0,0;
       0,0,0;
       0,0,0];
bsh35=[0,0,0;
       0,0,0;
       0,0,0];
Ysh35=gsh35+bsh35*j;
Y35_11=Y3_5+Ysh35;
Y35_12=-Y3_5;
Y35_21=-Y3_5;
Y35_22=Y3_5+Ysh35;
% branch 5
Z45s=0.005982400000000+0.005156750000000*j;
Z4_5 = [Z45s;Z45s;3*Z45s];
Z4_5s = A1*(diag([Z4_5(1,1),Z4_5(2,1),Z4_5(3,1)],0))*inv(A1);
Y4_5 = inv(Z4_5s);
gsh45=[0,0,0;
       0,0,0;
       0,0,0];
bsh45=[0,0,0;
       0,0,0;
       0,0,0];
Ysh45=gsh45+bsh45*j;
Y45_11=Y4_5+Ysh45;
Y45_12=-Y4_5;
Y45_21=-Y4_5;
Y45_22=Y4_5+Ysh45;

% branch 6
Z23s=0.005821062500000+0.004874362500000*j;
Z2_3 = [Z23s;Z23s;3*Z23s];
Z2_3s = A1*(diag([Z2_3(1,1),Z2_3(2,1),Z2_3(3,1)],0))*inv(A1);
Y2_3 = inv(Z2_3s);
gsh23=[0,0,0;
       0,0,0;
       0,0,0];
bsh23=[0,0,0;
       0,0,0;
       0,0,0];
Ysh23=gsh23+bsh23*j;
Y23_11=Y2_3+Ysh23;
Y23_12=-Y2_3;
Y23_21=-Y2_3;
Y23_22=Y2_3+Ysh23;

% branch 7
Z26s=0.001+0.0008*j;
Z2_6 = [Z26s;Z26s;3*Z26s];
Z2_6s = A1*(diag([Z2_6(1,1),Z2_6(2,1),Z2_6(3,1)],0))*inv(A1);
Y2_6 = inv(Z2_6s);
gsh26=[0,0,0;
       0,0,0;
       0,0,0];
bsh26=[0,0,0;
       0,0,0;
       0,0,0];
Ysh26=gsh26+bsh26*j;
Y26_11=Y2_6+Ysh26;
Y26_12=-Y2_6;
Y26_21=-Y2_6;
Y26_22=Y2_6+Ysh26;

Zdg=0.1+0.1*j;
Qdg=0.3;
Pdg=0.3;

pL2=[1.6;1.1;1.5];
qL2=[0.6;0.2;0.4];
pL3=[1.8;1.2;1.6];
qL3=[0.6;0.3;0.9];
pL4=[0.9;0.3;0.6];
qL4=[0.2;0.1;0.15];

Vbase = 1;
Ibase = 1;

V_lb = 0.95 * Vbase;
V_ub = 1.05 * Vbase;
v_lb = V_lb * V_lb;
v_ub = V_ub * V_ub;

su=[25;25;25];

pGl=0;
pGu=50;
qGl=-0.3;
qGu=18;
Transl = 0.95 * 0.95;
Transu = 1.1 * 1.1;
pDGl=0;
pDGu=0.1;
qDGl=0;
qDGu=0.1;


cvx_begin

    variable vG(3,3) hermitian
    variable v1(3,3) hermitian
    variables pT1(3) qT1(3) Trans(3)
    variables pf(21) qf(21) pt(21) qt(21)
    variables pG(6) qG(6)
    variable v(18,18) hermitian 
    variable v2seq(3,3) hermitian
    variable Idg2seq(3,3) hermitian
    variables pDG(3) qDG(3)
    minimize ( 100*(4*(pG(1)+pG(2)+pG(3))+(pG(4)+pG(5)+pG(6))) )
    subject to
        % SDR
        vG == hermitian_semidefinite(3)
        v1 == hermitian_semidefinite(3)
        v == hermitian_semidefinite(18)
        v1==[v(1,1) v(1,2) v(1,3);v(2,1) v(2,2) v(2,3);v(3,1) v(3,2) v(3,3)]


        % Vg = T*V1
Transl.*diag(v1) + [v_lb;v_lb;v_lb].*Trans-Transl.*[v_lb;v_lb;v_lb] <= diag(vG)
Transu.*diag(v1) + [v_ub;v_ub;v_ub].*Trans-Transu.*[v_ub;v_ub;v_ub] <= diag(vG)
diag(vG) <= Transl.*diag(v1) + [v_ub;v_ub;v_ub].*Trans-Transl.*[v_ub;v_ub;v_ub]
diag(vG) <= Transu.*diag(v1) + [v_lb;v_lb;v_lb].*Trans-Transu.*[v_lb;v_lb;v_lb]

% bus G
pG(1)==-pT1(1)
pG(2)==-pT1(2)
pG(3)==-pT1(3)

qG(1)==-qT1(1)
qG(2)==-qT1(2)
qG(3)==-qT1(3)

% bus 1
-pT1(1)==pf(1)+pf(4);
-pT1(2)==pf(2)+pf(5);
-pT1(3)==pf(3)+pf(6);

-qT1(1)==qf(1)+qf(4);
-qT1(2)==qf(2)+qf(5);
-qT1(3)==qf(3)+qf(6);

% bus 2
pt(1)==pf(7)+pf(16)+pf(19)+pL2(1);
pt(2)==pf(8)+pf(17)+pf(20)+pL2(2);
pt(3)==pf(9)+pf(18)+pf(21)+pL2(3);

qt(1)==qf(7)+qf(16)+qf(19)+qL2(1);
qt(2)==qf(8)+qf(17)+qf(20)+qL2(2);
qt(3)==qf(9)+qf(18)+qf(21)+qL2(3);

% bus 3
pt(4)+pt(16)==pf(10)+pL3(1);
pt(5)+pt(17)==pf(11)+pL3(2);
pt(6)+pt(18)==pf(12)+pL3(3);

qt(4)+qt(16)==qf(10)+qL3(1);
qt(5)+qt(17)==qf(11)+qL3(2);
qt(6)+qt(18)==qf(12)+qL3(3);

% bus 4
pt(7)==pf(13)+pL4(1);
pt(8)==pf(14)+pL4(2);
pt(9)==pf(15)+pL4(3);

qt(7)==qf(13)+qL4(1);
qt(8)==qf(14)+qL4(2);
qt(9)==qf(15)+qL4(3);

% bus 5
pG(4)==-pt(10)-pt(13);
pG(5)==-pt(11)-pt(14);
pG(6)==-pt(12)-pt(15);

qG(4)==-qt(10)-qt(13);
qG(5)==-qt(11)-qt(14);
qG(6)==-qt(12)-qt(15);

%bus 6
pDG(1)==-pt(19)
pDG(2)==-pt(20)
pDG(3)==-pt(21)

qDG(1)==-qt(19)
qDG(2)==-qt(20)
qDG(3)==-qt(21)

pDGl<=pDG(1)<=pDGu;
pDGl<=pDG(2)<=pDGu;
pDGl<=pDG(3)<=pDGu;

qDGl<=qDG(1)<=qDGu;
qDGl<=qDG(2)<=qDGu;
qDGl<=qDG(3)<=qDGu;

v_lb <= v(1,1) <= v_ub;
v_lb <= v(2,2) <= v_ub;
v_lb <= v(3,3) <= v_ub;
v_lb <= v(4,4) <= v_ub;
v_lb <= v(5,5) <= v_ub;
v_lb <= v(6,6) <= v_ub;
v_lb <= v(7,7) <= v_ub;
v_lb <= v(8,8) <= v_ub;
v_lb <= v(9,9) <= v_ub;
v_lb <= v(10,10) <= v_ub;
v_lb <= v(11,11) <= v_ub;
v_lb <= v(12,12) <= v_ub;
v_lb <= v(13,13) <= v_ub;
v_lb <= v(14,14) <= v_ub;
v_lb <= v(15,15) <= v_ub;
v_lb <= v(16,16) <= v_ub;
v_lb <= v(17,17) <= v_ub;
v_lb <= v(18,18) <= v_ub;

v_lb <= vG(1,1) <= v_ub;
v_lb <= vG(2,2) <= v_ub;
v_lb <= vG(3,3) <= v_ub;

pGl <= pG(1,1) <= pGu;
pGl <= pG(2,1) <= pGu;
pGl <= pG(3,1) <= pGu;
pGl <= pG(4,1) <= pGu;
pGl <= pG(5,1) <= pGu;
pGl <= pG(6,1) <= pGu;

qGl <= qG(1,1) <= qGu;
qGl <= qG(2,1) <= qGu;
qGl <= qG(3,1) <= qGu;
qGl <= qG(4,1) <= qGu;
qGl <= qG(5,1) <= qGu;
qGl <= qG(6,1) <= qGu;


% branch 1
pf(1) + qf(1)*j == conj(Y12_11(1,1))*v(1,1)+conj(Y12_11(1,2))*v(1,2)+conj(Y12_11(1,3))*v(1,3)+conj(Y12_12(1,1))*v(1,4)+conj(Y12_12(1,2))*v(1,5)+conj(Y12_12(1,3))*v(1,6);
pf(2) + qf(2)*j == conj(Y12_11(2,1))*v(2,1)+conj(Y12_11(2,2))*v(2,2)+conj(Y12_11(2,3))*v(2,3)+conj(Y12_12(2,1))*v(2,4)+conj(Y12_12(2,2))*v(2,5)+conj(Y12_12(2,3))*v(2,6);
pf(3) + qf(3)*j == conj(Y12_11(3,1))*v(3,1)+conj(Y12_11(3,2))*v(3,2)+conj(Y12_11(3,3))*v(3,3)+conj(Y12_12(3,1))*v(3,4)+conj(Y12_12(3,2))*v(3,5)+conj(Y12_12(3,3))*v(3,6);

-pt(1,1) - qt(1,1)*j == conj(Y12_21(1,1))*v(4,1)+conj(Y12_21(1,2))*v(4,2)+conj(Y12_21(1,3))*v(4,3)+conj(Y12_22(1,1))*v(4,4)+conj(Y12_22(1,2))*v(4,5)+conj(Y12_22(1,3))*v(4,6);
-pt(2,1) - qt(2,1)*j == conj(Y12_21(2,1))*v(5,1)+conj(Y12_21(2,2))*v(5,2)+conj(Y12_21(2,3))*v(5,3)+conj(Y12_22(2,1))*v(5,4)+conj(Y12_22(2,2))*v(5,5)+conj(Y12_22(2,3))*v(5,6);
-pt(3,1) - qt(3,1)*j == conj(Y12_21(3,1))*v(6,1)+conj(Y12_21(3,2))*v(6,2)+conj(Y12_21(3,3))*v(6,3)+conj(Y12_22(3,1))*v(6,4)+conj(Y12_22(3,2))*v(6,5)+conj(Y12_22(3,3))*v(6,6);

% branch 2
pf(4) + qf(4)*j == conj(Y13_11(1,1))*v(1,1)+conj(Y13_11(1,2))*v(1,2)+conj(Y13_11(1,3))*v(1,3)+conj(Y13_12(1,1))*v(1,7)+conj(Y13_12(1,2))*v(1,8)+conj(Y13_12(1,3))*v(1,9);
pf(5) + qf(5)*j == conj(Y13_11(2,1))*v(2,1)+conj(Y13_11(2,2))*v(2,2)+conj(Y13_11(2,3))*v(2,3)+conj(Y13_12(2,1))*v(2,7)+conj(Y13_12(2,2))*v(2,8)+conj(Y13_12(2,3))*v(2,9);
pf(6) + qf(6)*j == conj(Y13_11(3,1))*v(3,1)+conj(Y13_11(3,2))*v(3,2)+conj(Y13_11(3,3))*v(3,3)+conj(Y13_12(3,1))*v(3,7)+conj(Y13_12(3,2))*v(3,8)+conj(Y13_12(3,3))*v(3,9);

-pt(4,1) - qt(4,1)*j == conj(Y13_21(1,1))*v(7,1)+conj(Y13_21(1,2))*v(7,2)+conj(Y13_21(1,3))*v(7,3)+conj(Y13_22(1,1))*v(7,7)+conj(Y13_22(1,2))*v(7,8)+conj(Y13_22(1,3))*v(7,9);
-pt(5,1) - qt(5,1)*j == conj(Y13_21(2,1))*v(8,1)+conj(Y13_21(2,2))*v(8,2)+conj(Y13_21(2,3))*v(8,3)+conj(Y13_22(2,1))*v(8,7)+conj(Y13_22(2,2))*v(8,8)+conj(Y13_22(2,3))*v(8,9);
-pt(6,1) - qt(6,1)*j == conj(Y13_21(3,1))*v(9,1)+conj(Y13_21(3,2))*v(9,2)+conj(Y13_21(3,3))*v(9,3)+conj(Y13_22(3,1))*v(9,7)+conj(Y13_22(3,2))*v(9,8)+conj(Y13_22(3,3))*v(9,9);

% branch 3
pf(7) + qf(7)*j == conj(Y24_11(1,1))*v(4,4)+conj(Y24_11(1,2))*v(4,5)+conj(Y24_11(1,3))*v(4,6)+conj(Y24_12(1,1))*v(4,10)+conj(Y24_12(1,2))*v(4,11)+conj(Y24_12(1,3))*v(4,12);
pf(8) + qf(8)*j == conj(Y24_11(2,1))*v(5,4)+conj(Y24_11(2,2))*v(5,5)+conj(Y24_11(2,3))*v(5,6)+conj(Y24_12(2,1))*v(5,10)+conj(Y24_12(2,2))*v(5,11)+conj(Y24_12(2,3))*v(5,12);
pf(9) + qf(9)*j == conj(Y24_11(3,1))*v(6,4)+conj(Y24_11(3,2))*v(6,5)+conj(Y24_11(3,3))*v(6,6)+conj(Y24_12(3,1))*v(6,10)+conj(Y24_12(3,2))*v(6,11)+conj(Y24_12(3,3))*v(6,12);

-pt(7,1) - qt(7,1)*j == conj(Y24_21(1,1))*v(10,4)+conj(Y24_21(1,2))*v(10,5)+conj(Y24_21(1,3))*v(10,6)+conj(Y24_22(1,1))*v(10,10)+conj(Y24_22(1,2))*v(10,11)+conj(Y24_22(1,3))*v(10,12);
-pt(8,1) - qt(8,1)*j == conj(Y24_21(2,1))*v(11,4)+conj(Y24_21(2,2))*v(11,5)+conj(Y24_21(2,3))*v(11,6)+conj(Y24_22(2,1))*v(11,10)+conj(Y24_22(2,2))*v(11,11)+conj(Y24_22(2,3))*v(11,12);
-pt(9,1) - qt(9,1)*j == conj(Y24_21(3,1))*v(12,4)+conj(Y24_21(3,2))*v(12,5)+conj(Y24_21(3,3))*v(12,6)+conj(Y24_22(3,1))*v(12,10)+conj(Y24_22(3,2))*v(12,11)+conj(Y24_22(3,3))*v(12,12);

% branch 4
pf(10) + qf(10)*j == conj(Y35_11(1,1))*v(7,7)+conj(Y35_11(1,2))*v(7,8)+conj(Y35_11(1,3))*v(7,9)+conj(Y35_12(1,1))*v(7,13)+conj(Y35_12(1,2))*v(7,14)+conj(Y35_12(1,3))*v(7,15);
pf(11) + qf(11)*j == conj(Y35_11(2,1))*v(8,7)+conj(Y35_11(2,2))*v(8,8)+conj(Y35_11(2,3))*v(8,9)+conj(Y35_12(2,1))*v(8,13)+conj(Y35_12(2,2))*v(8,14)+conj(Y35_12(2,3))*v(8,15);
pf(12) + qf(12)*j == conj(Y35_11(3,1))*v(9,7)+conj(Y35_11(3,2))*v(9,8)+conj(Y35_11(3,3))*v(9,9)+conj(Y35_12(3,1))*v(9,13)+conj(Y35_12(3,2))*v(9,14)+conj(Y35_12(3,3))*v(9,15);

-pt(10,1) - qt(10,1)*j == conj(Y35_21(1,1))*v(13,7)+conj(Y35_21(1,2))*v(13,8)+conj(Y35_21(1,3))*v(13,9)+conj(Y35_22(1,1))*v(13,13)+conj(Y35_22(1,2))*v(13,14)+conj(Y35_22(1,3))*v(13,15);
-pt(11,1) - qt(11,1)*j == conj(Y35_21(2,1))*v(14,7)+conj(Y35_21(2,2))*v(14,8)+conj(Y35_21(2,3))*v(14,9)+conj(Y35_22(2,1))*v(14,13)+conj(Y35_22(2,2))*v(14,14)+conj(Y35_22(2,3))*v(14,15);
-pt(12,1) - qt(12,1)*j == conj(Y35_21(3,1))*v(15,7)+conj(Y35_21(3,2))*v(15,8)+conj(Y35_21(3,3))*v(15,9)+conj(Y35_22(3,1))*v(15,13)+conj(Y35_22(3,2))*v(15,14)+conj(Y35_22(3,3))*v(15,15);

% branch 5
pf(13) + qf(13)*j == conj(Y45_11(1,1))*v(10,10)+conj(Y45_11(1,2))*v(10,11)+conj(Y45_11(1,3))*v(10,12)+conj(Y45_12(1,1))*v(10,13)+conj(Y45_12(1,2))*v(10,14)+conj(Y45_12(1,3))*v(10,15);
pf(14) + qf(14)*j == conj(Y45_11(2,1))*v(11,10)+conj(Y45_11(2,2))*v(11,11)+conj(Y45_11(2,3))*v(11,12)+conj(Y45_12(2,1))*v(11,13)+conj(Y45_12(2,2))*v(11,14)+conj(Y45_12(2,3))*v(11,15);
pf(15) + qf(15)*j == conj(Y45_11(3,1))*v(12,10)+conj(Y45_11(3,2))*v(12,11)+conj(Y45_11(3,3))*v(12,12)+conj(Y45_12(3,1))*v(12,13)+conj(Y45_12(3,2))*v(12,14)+conj(Y45_12(3,3))*v(12,15);

-pt(13,1) - qt(13,1)*j == conj(Y45_21(1,1))*v(13,10)+conj(Y45_21(1,2))*v(13,11)+conj(Y45_21(1,3))*v(13,12)+conj(Y45_22(1,1))*v(13,13)+conj(Y45_22(1,2))*v(13,14)+conj(Y45_22(1,3))*v(13,15);
-pt(14,1) - qt(14,1)*j == conj(Y45_21(2,1))*v(14,10)+conj(Y45_21(2,2))*v(14,11)+conj(Y45_21(2,3))*v(14,12)+conj(Y45_22(2,1))*v(14,13)+conj(Y45_22(2,2))*v(14,14)+conj(Y45_22(2,3))*v(14,15);
-pt(15,1) - qt(15,1)*j == conj(Y45_21(3,1))*v(15,10)+conj(Y45_21(3,2))*v(15,11)+conj(Y45_21(3,3))*v(15,12)+conj(Y45_22(3,1))*v(15,13)+conj(Y45_22(3,2))*v(15,14)+conj(Y45_22(3,3))*v(15,15);
% branch 6
pf(16) + qf(16)*j == conj(Y23_11(1,1))*v(4,4)+conj(Y23_11(1,2))*v(4,5)+conj(Y23_11(1,3))*v(4,6)+conj(Y23_12(1,1))*v(4,7)+conj(Y23_12(1,2))*v(4,8)+conj(Y23_12(1,3))*v(4,9);
pf(17) + qf(17)*j == conj(Y23_11(2,1))*v(5,4)+conj(Y23_11(2,2))*v(5,5)+conj(Y23_11(2,3))*v(5,6)+conj(Y23_12(2,1))*v(5,7)+conj(Y23_12(2,2))*v(5,8)+conj(Y23_12(2,3))*v(5,9);
pf(18) + qf(18)*j == conj(Y23_11(3,1))*v(6,4)+conj(Y23_11(3,2))*v(6,5)+conj(Y23_11(3,3))*v(6,6)+conj(Y23_12(3,1))*v(6,7)+conj(Y23_12(3,2))*v(6,8)+conj(Y23_12(3,3))*v(6,9);

-pt(16,1) - qt(16,1)*j == conj(Y23_21(1,1))*v(7,4)+conj(Y23_21(1,2))*v(7,5)+conj(Y23_21(1,3))*v(7,6)+conj(Y23_22(1,1))*v(7,7)+conj(Y23_22(1,2))*v(7,8)+conj(Y23_22(1,3))*v(7,9);
-pt(17,1) - qt(17,1)*j == conj(Y23_21(2,1))*v(8,4)+conj(Y23_21(2,2))*v(8,5)+conj(Y23_21(2,3))*v(8,6)+conj(Y23_22(2,1))*v(8,7)+conj(Y23_22(2,2))*v(8,8)+conj(Y23_22(2,3))*v(8,9);
-pt(18,1) - qt(18,1)*j == conj(Y23_21(3,1))*v(9,4)+conj(Y23_21(3,2))*v(9,5)+conj(Y23_21(3,3))*v(9,6)+conj(Y23_22(3,1))*v(9,7)+conj(Y23_22(3,2))*v(9,8)+conj(Y23_22(3,3))*v(9,9);

% branch 7
pf(19) + qf(19)*j == conj(Y26_11(1,1))*v(4,4)+conj(Y26_11(1,2))*v(4,5)+conj(Y26_11(1,3))*v(4,6)+conj(Y26_12(1,1))*v(4,16)+conj(Y26_12(1,2))*v(4,17)+conj(Y26_12(1,3))*v(4,18);
pf(20) + qf(20)*j == conj(Y26_11(2,1))*v(5,4)+conj(Y26_11(2,2))*v(5,5)+conj(Y26_11(2,3))*v(5,6)+conj(Y26_12(2,1))*v(5,16)+conj(Y26_12(2,2))*v(5,17)+conj(Y26_12(2,3))*v(5,18);
pf(21) + qf(21)*j == conj(Y26_11(3,1))*v(6,4)+conj(Y26_11(3,2))*v(6,5)+conj(Y26_11(3,3))*v(6,6)+conj(Y26_12(3,1))*v(6,16)+conj(Y26_12(3,2))*v(6,17)+conj(Y26_12(3,3))*v(6,18);

-pt(19,1) - qt(19,1)*j == conj(Y26_21(1,1))*v(16,4)+conj(Y26_21(1,2))*v(16,5)+conj(Y26_21(1,3))*v(16,6)+conj(Y26_22(1,1))*v(16,16)+conj(Y26_22(1,2))*v(16,17)+conj(Y26_22(1,3))*v(16,18);
-pt(20,1) - qt(20,1)*j == conj(Y26_21(2,1))*v(17,4)+conj(Y26_21(2,2))*v(17,5)+conj(Y26_21(2,3))*v(17,6)+conj(Y26_22(2,1))*v(17,16)+conj(Y26_22(2,2))*v(17,17)+conj(Y26_22(2,3))*v(17,18);
-pt(21,1) - qt(21,1)*j == conj(Y26_21(3,1))*v(18,4)+conj(Y26_21(3,2))*v(18,5)+conj(Y26_21(3,3))*v(18,6)+conj(Y26_22(3,1))*v(18,16)+conj(Y26_22(3,2))*v(18,17)+conj(Y26_22(3,3))*v(18,18);

[su(1) pf(1) + 1j*qf(1); pf(1) - 1j*qf(1) su(1)] == hermitian_semidefinite(2)
[su(1) pt(1) + 1j*qt(1); pt(1) - 1j*qt(1) su(1)] == hermitian_semidefinite(2)
[su(2) pf(2) + 1j*qf(2); pf(2) - 1j*qf(2) su(2)] == hermitian_semidefinite(2)
[su(2) pt(2) + 1j*qt(2); pt(2) - 1j*qt(2) su(2)] == hermitian_semidefinite(2)
[su(3) pf(3) + 1j*qf(3); pf(3) - 1j*qf(3) su(3)] == hermitian_semidefinite(2)
[su(3) pt(3) + 1j*qt(3); pt(3) - 1j*qt(3) su(3)] == hermitian_semidefinite(2)

[su(1) pf(4) + 1j*qf(4); pf(4) - 1j*qf(4) su(1)] == hermitian_semidefinite(2)
[su(1) pt(4) + 1j*qt(4); pt(4) - 1j*qt(4) su(1)] == hermitian_semidefinite(2)
[su(2) pf(5) + 1j*qf(5); pf(5) - 1j*qf(5) su(2)] == hermitian_semidefinite(2)
[su(2) pt(5) + 1j*qt(5); pt(5) - 1j*qt(5) su(2)] == hermitian_semidefinite(2)
[su(3) pf(6) + 1j*qf(6); pf(6) - 1j*qf(6) su(3)] == hermitian_semidefinite(2)
[su(3) pt(6) + 1j*qt(6); pt(6) - 1j*qt(6) su(3)] == hermitian_semidefinite(2)

[su(1) pf(7) + 1j*qf(7); pf(7) - 1j*qf(7) su(1)] == hermitian_semidefinite(2)
[su(1) pt(7) + 1j*qt(7); pt(7) - 1j*qt(7) su(1)] == hermitian_semidefinite(2)
[su(2) pf(8) + 1j*qf(8); pf(8) - 1j*qf(8) su(2)] == hermitian_semidefinite(2)
[su(2) pt(8) + 1j*qt(8); pt(8) - 1j*qt(8) su(2)] == hermitian_semidefinite(2)
[su(3) pf(9) + 1j*qf(9); pf(9) - 1j*qf(9) su(3)] == hermitian_semidefinite(2)
[su(3) pt(9) + 1j*qt(9); pt(9) - 1j*qt(9) su(3)] == hermitian_semidefinite(2)

[su(1) pf(10) + 1j*qf(10); pf(10) - 1j*qf(10) su(1)] == hermitian_semidefinite(2)
[su(1) pt(10) + 1j*qt(10); pt(10) - 1j*qt(10) su(1)] == hermitian_semidefinite(2)
[su(2) pf(11) + 1j*qf(11); pf(11) - 1j*qf(11) su(2)] == hermitian_semidefinite(2)
[su(2) pt(11) + 1j*qt(11); pt(11) - 1j*qt(11) su(2)] == hermitian_semidefinite(2)
[su(3) pf(12) + 1j*qf(12); pf(12) - 1j*qf(12) su(3)] == hermitian_semidefinite(2)
[su(3) pt(12) + 1j*qt(12); pt(12) - 1j*qt(12) su(3)] == hermitian_semidefinite(2)

[su(1) pf(13) + 1j*qf(13); pf(13) - 1j*qf(13) su(1)] == hermitian_semidefinite(2)
[su(1) pt(13) + 1j*qt(13); pt(13) - 1j*qt(13) su(1)] == hermitian_semidefinite(2)
[su(2) pf(14) + 1j*qf(14); pf(14) - 1j*qf(14) su(2)] == hermitian_semidefinite(2)
[su(2) pt(14) + 1j*qt(14); pt(14) - 1j*qt(14) su(2)] == hermitian_semidefinite(2)
[su(3) pf(15) + 1j*qf(15); pf(15) - 1j*qf(15) su(3)] == hermitian_semidefinite(2)
[su(3) pt(15) + 1j*qt(15); pt(15) - 1j*qt(15) su(3)] == hermitian_semidefinite(2)

[su(1) pf(16) + 1j*qf(16); pf(16) - 1j*qf(16) su(1)] == hermitian_semidefinite(2)
[su(1) pt(16) + 1j*qt(16); pt(16) - 1j*qt(16) su(1)] == hermitian_semidefinite(2)
[su(2) pf(17) + 1j*qf(17); pf(17) - 1j*qf(17) su(2)] == hermitian_semidefinite(2)
[su(2) pt(17) + 1j*qt(17); pt(17) - 1j*qt(17) su(2)] == hermitian_semidefinite(2)
[su(3) pf(18) + 1j*qf(18); pf(18) - 1j*qf(18) su(3)] == hermitian_semidefinite(2)
[su(3) pt(18) + 1j*qt(18); pt(18) - 1j*qt(18) su(3)] == hermitian_semidefinite(2)

[su(1) pf(19) + 1j*qf(19); pf(19) - 1j*qf(19) su(1)] == hermitian_semidefinite(2)
[su(1) pt(19) + 1j*qt(19); pt(19) - 1j*qt(19) su(1)] == hermitian_semidefinite(2)
[su(2) pf(20) + 1j*qf(20); pf(20) - 1j*qf(20) su(2)] == hermitian_semidefinite(2)
[su(2) pt(20) + 1j*qt(20); pt(20) - 1j*qt(20) su(2)] == hermitian_semidefinite(2)
[su(3) pf(21) + 1j*qf(21); pf(21) - 1j*qf(21) su(3)] == hermitian_semidefinite(2)
[su(3) pt(21) + 1j*qt(21); pt(21) - 1j*qt(21) su(3)] == hermitian_semidefinite(2)

cvx_end

% Optimal solution
cpusdr = cvx_cputime;
optsdr = cvx_optval;
statussdr = cvx_status;
Vsdr = v
Vmag=sqrt(diag(Vsdr))
angle1A=0
angle1B=-2.094395102393195
angle1C=2.094395102393195
angle2A=angle1A+atan2(imag(v(4,1)),real(v(4,1)))
angle2B=angle1B+atan2(imag(v(5,2)),real(v(5,2)))
angle2C=angle1C+atan2(imag(v(6,3)),real(v(6,3)))
angle3A=angle1A+atan2(imag(v(7,1)),real(v(7,1)))
angle3B=angle1B+atan2(imag(v(8,2)),real(v(8,2)))
angle3C=angle1C+atan2(imag(v(9,3)),real(v(9,3)))
angle4A=angle2A+atan2(imag(v(10,4)),real(v(10,4)))
angle4B=angle2B+atan2(imag(v(11,5)),real(v(11,5)))
angle4C=angle2C+atan2(imag(v(12,6)),real(v(12,6)))
angle5A=angle3A+atan2(imag(v(13,7)),real(v(13,7)))
angle5B=angle3B+atan2(imag(v(14,8)),real(v(14,8)))
angle5C=angle3C+atan2(imag(v(15,9)),real(v(15,9)))
angle6A=angle2A+atan2(imag(v(16,4)),real(v(16,4)))
angle6B=angle2B+atan2(imag(v(17,5)),real(v(17,5)))
angle6C=angle2C+atan2(imag(v(18,6)),real(v(18,6)))
angle7A=angle1A;
angle7B=angle1B;
angle7C=angle1C;

% vsdr = approx_volt_profile(Adj,Vsdr,slack,angslack)
pG=value(pG)
qG=value(qG)
pf=value(pf)
qf=value(qf)
pt=value(pt)
qt=value(qt)
pDG = value(pDG)
qDG = value(qDG)
vG=sqrt(diag(vG))
sqrt(Trans)

A=[angle1A;angle1B;angle1C;angle2A;angle2B;angle2C;angle3A;angle3B;angle3C;angle4A;angle4B;angle4C;angle5A;angle5B;angle5C;angle6A;angle6B;angle6C;angle7A;angle7B;angle7C];