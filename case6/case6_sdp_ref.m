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
Z12s=0.002608000000000+0.003189000000000*j;
Z1_2 = [Z12s;Z12s;3*Z12s];
Z12 = A1*(diag([Z1_2(1,1),Z1_2(2,1),Z1_2(3,1)],0))*inv(A1);
Y1_2 = inv(Z12);  % 三相支路导纳
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
Z13s=0.003683750000000+0.004634062500000*j;
Z1_3 = [Z13s;Z13s;3*Z13s];
Z13 = A1*(diag([Z1_3(1,1),Z1_3(2,1),Z1_3(3,1)],0))*inv(A1);
Y1_3 = inv(Z13);
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
Z24s=0.022912875000000+0.039472500000000*j;
Z2_4 = [Z24s;Z24s;3*Z24s];
Z24 = A1*(diag([Z2_4(1,1),Z2_4(2,1),Z2_4(3,1)],0))*inv(A1);
Y2_4 = inv(Z24);
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
Z35s=0.022912875000000+0.039472500000000*j;
Z3_5 = [Z35s;Z35s;3*Z35s];
Z35 = A1*(diag([Z3_5(1,1),Z3_5(2,1),Z3_5(3,1)],0))*inv(A1);
Y3_5 = inv(Z35);
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

pL2=[1.3;1.3;1.3];
qL2=[0.2;0.2;0.2];
pL3=[1.3;1.3;1.3];
qL3=[0.2;0.2;0.2];
pL4=[0.65;0.65;0.65];
qL4=[0.1;0.1;0.1];

v0=[
     1.005582856229995+0j;
     1.005582856230481*exp(-2/3*pi*j);
     1.005582856226576*exp(2/3*pi*j);
     ];


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
tic
cvx_begin

    variable vG(3,3) hermitian 
    variable v11(3,3) hermitian
    variable v12(3,3) complex
    variable v21(3,3) complex
    variable v22(3,3) hermitian
    variable v13(3,3) complex 
    variable v31(3,3) complex
    variable v33(3,3) hermitian
    variable v24(3,3) complex
    variable v42(3,3) complex 
    variable v44(3,3) hermitian
    variable v35(3,3) complex 
    variable v53(3,3) complex
    variable v55(3,3) hermitian
    
    variable I12(3,3) hermitian
    variable I13(3,3) hermitian
    variable I24(3,3) hermitian
    variable I35(3,3) hermitian

    variable S12(3,3) complex
    variable S13(3,3) complex
    variable S24(3,3) complex
    variable S35(3,3) complex

    variable DeltaS12(3,3) complex
    variable DeltaS13(3,3) complex
    variable DeltaS24(3,3) complex
    variable DeltaS35(3,3) complex
    
    variable S21(3,3) complex
    variable S31(3,3) complex
    variable S42(3,3) complex
    variable S53(3,3) complex
    
    variable Sg1(3,1) complex
    variable Sg2(3,1) complex
    variable Sg3(3,1) complex
    variable Sg4(3,1) complex
    variable Sg5(3,1) complex

%     variables Trans(3)
    Trans = diag([1.029711386463642;1.029711386464766;1.029711386459019],0);

    minimize ( 100*(4*(real(Sg1(1,1))+real(Sg1(2,1))+real(Sg1(3,1)))+(real(Sg5(1,1))+real(Sg5(2,1))+real(Sg5(3,1)))))
    subject to
        % SDR
        vG == hermitian_semidefinite(3)
        v11 == hermitian_semidefinite(3)
        v22 == hermitian_semidefinite(3)
        v33 == hermitian_semidefinite(3)
        v44 == hermitian_semidefinite(3)
        v55 == hermitian_semidefinite(3)
        
        I12 == hermitian_semidefinite(3)
        I13 == hermitian_semidefinite(3)
        I24 == hermitian_semidefinite(3)
        I35 == hermitian_semidefinite(3)
        

        % Vg = T*V1        

vG == v0 * ctranspose(v0); 

vG == Trans * v11 * ctranspose(Trans);

Sg1 == diag(S12) + diag(S13);
Sg2 == diag(S21) + diag(S24) +(pL2 + qL2*j);
Sg3 == diag(S31) + diag(S35) +(pL3 + qL3*j);
Sg4 == diag(S42) + (pL4 + qL4*j);
Sg5 == diag(S53);

% 1-2

DeltaS12 == Z12*I12;
diag(S12) - diag(DeltaS12)==-diag(S21);
v22 == v11 - S12*ctranspose(Z12) - Z12*ctranspose(S12) + Z12*I12*ctranspose(Z12);
v11 == v22 - S21*ctranspose(Z12) - Z12*ctranspose(S21) + Z12*I12*ctranspose(Z12);

S12(1,1) == conj(Y12_11(1,1))*v11(1,1)+conj(Y12_11(1,2))*v11(1,2)+conj(Y12_11(1,3))*v11(1,3)+conj(Y12_12(1,1))*v12(1,1)+conj(Y12_12(1,2))*v12(1,2)+conj(Y12_12(1,3))*v12(1,3);
S12(2,2) == conj(Y12_11(2,1))*v11(2,1)+conj(Y12_11(2,2))*v11(2,2)+conj(Y12_11(2,3))*v11(2,3)+conj(Y12_12(2,1))*v12(2,1)+conj(Y12_12(2,2))*v12(2,2)+conj(Y12_12(2,3))*v12(2,3);
S12(3,3) == conj(Y12_11(3,1))*v11(3,1)+conj(Y12_11(3,2))*v11(3,2)+conj(Y12_11(3,3))*v11(3,3)+conj(Y12_12(3,1))*v12(3,1)+conj(Y12_12(3,2))*v12(3,2)+conj(Y12_12(3,3))*v12(3,3);

S21(1,1) == conj(Y12_21(1,1))*v21(1,1)+conj(Y12_21(1,2))*v21(1,2)+conj(Y12_21(1,3))*v21(1,3)+conj(Y12_22(1,1))*v22(1,1)+conj(Y12_22(1,2))*v22(1,2)+conj(Y12_22(1,3))*v22(1,3);
S21(2,2) == conj(Y12_21(2,1))*v21(2,1)+conj(Y12_21(2,2))*v21(2,2)+conj(Y12_21(2,3))*v21(2,3)+conj(Y12_22(2,1))*v22(2,1)+conj(Y12_22(2,2))*v22(2,2)+conj(Y12_22(2,3))*v22(2,3);
S21(3,3) == conj(Y12_21(3,1))*v21(3,1)+conj(Y12_21(3,2))*v21(3,2)+conj(Y12_21(3,3))*v21(3,3)+conj(Y12_22(3,1))*v22(3,1)+conj(Y12_22(3,2))*v22(3,2)+conj(Y12_22(3,3))*v22(3,3);

% 2-4
DeltaS24 == Z24*I24;
diag(S24) - diag(DeltaS24)== -diag(S42);

v44 == v22 - S24*ctranspose(Z24) - Z24*ctranspose(S24) + Z24*I24*ctranspose(Z24);
v22 == v44 - S42*ctranspose(Z24) - Z24*ctranspose(S42) + Z24*I24*ctranspose(Z24);

S24(1,1) == conj(Y24_11(1,1))*v22(1,1)+conj(Y24_11(1,2))*v22(1,2)+conj(Y24_11(1,3))*v22(1,3)+conj(Y24_12(1,1))*v24(1,1)+conj(Y24_12(1,2))*v24(1,2)+conj(Y24_12(1,3))*v24(1,3);
S24(2,1) == conj(Y24_11(2,1))*v22(2,1)+conj(Y24_11(2,2))*v22(2,2)+conj(Y24_11(2,3))*v22(2,3)+conj(Y24_12(2,1))*v24(2,1)+conj(Y24_12(2,2))*v24(2,2)+conj(Y24_12(2,3))*v24(2,3);
S24(3,1) == conj(Y24_11(3,1))*v22(3,1)+conj(Y24_11(3,2))*v22(3,2)+conj(Y24_11(3,3))*v22(3,3)+conj(Y24_12(3,1))*v24(3,1)+conj(Y24_12(3,2))*v24(3,2)+conj(Y24_12(3,3))*v24(3,3);

S42(1,1) == conj(Y24_21(1,1))*v42(1,1)+conj(Y24_21(1,2))*v42(1,2)+conj(Y24_21(1,3))*v42(1,3)+conj(Y24_22(1,1))*v44(1,1)+conj(Y24_22(1,2))*v44(1,2)+conj(Y24_22(1,3))*v44(1,3);
S42(2,1) == conj(Y24_21(2,1))*v42(2,1)+conj(Y24_21(2,2))*v42(2,2)+conj(Y24_21(2,3))*v42(2,3)+conj(Y24_22(2,1))*v44(2,1)+conj(Y24_22(2,2))*v44(2,2)+conj(Y24_22(2,3))*v44(2,3);
S42(3,1) == conj(Y24_21(3,1))*v42(3,1)+conj(Y24_21(3,2))*v42(3,2)+conj(Y24_21(3,3))*v42(3,3)+conj(Y24_22(3,1))*v44(3,1)+conj(Y24_22(3,2))*v44(3,2)+conj(Y24_22(3,3))*v44(3,3);

% 1-3
DeltaS13 == Z13*I13;
diag(S13) - diag(DeltaS13)== -diag(S31);
v33 == v11 - S13*ctranspose(Z13) - Z13*ctranspose(S13) + Z13*I13*ctranspose(Z13);
v11 == v33 - S31*ctranspose(Z13) - Z13*ctranspose(S31) + Z13*I13*ctranspose(Z13);

S13(1,1) == conj(Y13_11(1,1))*v11(1,1)+conj(Y13_11(1,2))*v11(1,2)+conj(Y13_11(1,3))*v11(1,3)+conj(Y13_12(1,1))*v13(1,1)+conj(Y13_12(1,2))*v13(1,2)+conj(Y13_12(1,3))*v13(1,3);
S13(2,1) == conj(Y13_11(2,1))*v11(2,1)+conj(Y13_11(2,2))*v11(2,2)+conj(Y13_11(2,3))*v11(2,3)+conj(Y13_12(2,1))*v13(2,1)+conj(Y13_12(2,2))*v13(2,2)+conj(Y13_12(2,3))*v13(2,3);
S13(3,1) == conj(Y13_11(3,1))*v11(3,1)+conj(Y13_11(3,2))*v11(3,2)+conj(Y13_11(3,3))*v11(3,3)+conj(Y13_12(3,1))*v13(3,1)+conj(Y13_12(3,2))*v13(3,2)+conj(Y13_12(3,3))*v13(3,3);

S31(1,1) == conj(Y13_21(1,1))*v31(1,1)+conj(Y13_21(1,2))*v31(1,2)+conj(Y13_21(1,3))*v31(1,3)+conj(Y13_22(1,1))*v33(1,1)+conj(Y13_22(1,2))*v33(1,2)+conj(Y13_22(1,3))*v33(1,3);
S31(2,1) == conj(Y13_21(2,1))*v31(2,1)+conj(Y13_21(2,2))*v31(2,2)+conj(Y13_21(2,3))*v31(2,3)+conj(Y13_22(2,1))*v33(2,1)+conj(Y13_22(2,2))*v33(2,2)+conj(Y13_22(2,3))*v33(2,3);
S31(3,1) == conj(Y13_21(3,1))*v31(3,1)+conj(Y13_21(3,2))*v31(3,2)+conj(Y13_21(3,3))*v31(3,3)+conj(Y13_22(3,1))*v33(3,1)+conj(Y13_22(3,2))*v33(3,2)+conj(Y13_22(3,3))*v33(3,3);

% 3-5
DeltaS35 == Z35*I35;
diag(S35) - diag(DeltaS35)== -diag(S53);
v55 == v33 - S35*ctranspose(Z35) - Z35*ctranspose(S35) + Z35*I35*ctranspose(Z35);
v33 == v55 - S53*ctranspose(Z35) - Z35*ctranspose(S53) + Z35*I35*ctranspose(Z35);

S35(1,1) == conj(Y35_11(1,1))*v33(1,1)+conj(Y35_11(1,2))*v33(1,2)+conj(Y35_11(1,3))*v33(1,3)+conj(Y35_12(1,1))*v35(1,1)+conj(Y35_12(1,2))*v35(1,2)+conj(Y35_12(1,3))*v35(1,3);
S35(2,1) == conj(Y35_11(2,1))*v33(2,1)+conj(Y35_11(2,2))*v33(2,2)+conj(Y35_11(2,3))*v33(2,3)+conj(Y35_12(2,1))*v35(2,1)+conj(Y35_12(2,2))*v35(2,2)+conj(Y35_12(2,3))*v35(2,3);
S35(3,1) == conj(Y35_11(3,1))*v33(3,1)+conj(Y35_11(3,2))*v33(3,2)+conj(Y35_11(3,3))*v33(3,3)+conj(Y35_12(3,1))*v35(3,1)+conj(Y35_12(3,2))*v35(3,2)+conj(Y35_12(3,3))*v35(3,3);

S53(1,1) == conj(Y35_21(1,1))*v53(1,1)+conj(Y35_21(1,2))*v53(1,2)+conj(Y35_21(1,3))*v53(1,3)+conj(Y35_22(1,1))*v55(1,1)+conj(Y35_22(1,2))*v55(1,2)+conj(Y35_22(1,3))*v55(1,3);
S53(2,1) == conj(Y35_21(2,1))*v53(2,1)+conj(Y35_21(2,2))*v53(2,2)+conj(Y35_21(2,3))*v53(2,3)+conj(Y35_22(2,1))*v55(2,1)+conj(Y35_22(2,2))*v55(2,2)+conj(Y35_22(2,3))*v55(2,3);
S53(3,1) == conj(Y35_21(3,1))*v53(3,1)+conj(Y35_21(3,2))*v53(3,2)+conj(Y35_21(3,3))*v53(3,3)+conj(Y35_22(3,1))*v55(3,1)+conj(Y35_22(3,2))*v55(3,2)+conj(Y35_22(3,3))*v55(3,3);

v_lb <= v11(1,1) <= v_ub;
v_lb <= v11(2,2) <= v_ub;
v_lb <= v11(3,3) <= v_ub;
v_lb <= v22(1,1) <= v_ub;
v_lb <= v22(2,2) <= v_ub;
v_lb <= v22(3,3) <= v_ub;
v_lb <= v33(1,1) <= v_ub;
v_lb <= v33(2,2) <= v_ub;
v_lb <= v33(3,3) <= v_ub;
v_lb <= v44(1,1) <= v_ub;
v_lb <= v44(2,2) <= v_ub;
v_lb <= v44(3,3) <= v_ub;
v_lb <= v55(1,1) <= v_ub;
v_lb <= v55(2,2) <= v_ub;
v_lb <= v55(3,3) <= v_ub;

v_lb <= vG(1,1) <= v_ub;
v_lb <= vG(2,2) <= v_ub;
v_lb <= vG(3,3) <= v_ub;

pGl <= real(Sg1(1,1)) <= pGu;
pGl <= real(Sg1(2,1)) <= pGu;
pGl <= real(Sg1(3,1)) <= pGu;

diag(Sg2,0) == zeros(3);
diag(Sg3,0) == zeros(3);
diag(Sg4,0) == zeros(3);

pGl <= real(Sg5(1,1)) <= pGu;
pGl <= real(Sg5(2,1)) <= pGu;
pGl <= real(Sg5(3,1)) <= pGu;

qGl <= imag(Sg1(1,1)) <= qGu;
qGl <= imag(Sg1(2,1)) <= qGu;
qGl <= imag(Sg1(3,1)) <= qGu;

qGl <= imag(Sg5(1,1)) <= qGu;
qGl <= imag(Sg5(2,1)) <= qGu;
qGl <= imag(Sg5(3,1)) <= qGu;

[v11, S12; ctranspose(S12), I12] == hermitian_semidefinite(6)
[v11, S13; ctranspose(S13), I13] == hermitian_semidefinite(6)
[v22, S24; ctranspose(S24), I24] == hermitian_semidefinite(6)
[v33, S35; ctranspose(S35), I35] == hermitian_semidefinite(6)

cvx_end
toc
% Optimal solution
cpusdr = cvx_cputime;
optsdr = cvx_optval;
statussdr = cvx_status;

vG=sqrt(diag(vG))
sqrt(Trans)

Sg1 = value(Sg1)
Sg2 = value(Sg2)
Sg3 = value(Sg3)
Sg4 = value(Sg4)
Sg5 = value(Sg5)

v1 = value(v11);
v2 = value(v22);
v3 = value(v33);
v4 = value(v44);
v5 = value(v55);

I12 = value(I12);
I13 = value(I13);
I24 = value(I24);
I35 = value(I35);

S12 = value(S12);
S13 = value(S13);
S24 = value(S24);
S35 = value(S35);

S21 = value(S21);
S31 = value(S31);
S42 = value(S42);
S53 = value(S53);

pG=[real(Sg1);real(Sg5)]
qG=[imag(Sg1);imag(Sg5)]

pf=[diag(real(S12));diag(real(S13));diag(real(S24));diag(real(S35))]
pf=[diag(imag(S12));diag(imag(S13));diag(imag(S24));diag(imag(S35))]

pt=(-1).*[diag(real(S21));diag(real(S31));diag(real(S42));diag(real(S53))]
qt=(-1).*[diag(imag(S21));diag(imag(S31));diag(imag(S42));diag(imag(S53))]

Vmag1 = diag(sqrt(value(v11)))
Vmag2 = diag(sqrt(value(v22)))
Vmag3 = diag(sqrt(value(v33)))
Vmag4 = diag(sqrt(value(v44)))
Vmag5 = diag(sqrt(value(v55)))

Vg = v0;
V1 = Vg ./ diag(Trans);
Angle1 = atan2(imag(V1),real(V1))
I12 = 1 ./ ( v1([1, 2, 3],[1, 2, 3]) ) .* ctranspose(S12) .* V1([1, 2, 3]);
V2 =V1([1, 2, 3]) -  Z12 *I12;
V2 = diag (V2)
Angle2 = atan2(imag(V2),real(V2))
I13 = 1 ./ ( v1([1, 2, 3],[1, 2, 3]) ) .* ctranspose(S13) .* V1([1, 2, 3]);
V3 =V1([1, 2, 3]) -  Z13 *I13;
V3 = diag (V3)
Angle3 = atan2(imag(V3),real(V3))
I24 = 1 ./ ( v2([1, 2, 3],[1, 2, 3]) ) .* ctranspose(S24) .* V2([1, 2, 3]);
V4 =V2([1, 2, 3]) -  Z24 *I24;
V4 = diag (V4)
Angle4 = atan2(imag(V4),real(V4))
I35 = 1 ./ ( v3([1, 2, 3],[1, 2, 3]) ) .* ctranspose(S35) .* V3([1, 2, 3]);
V5 =V3([1, 2, 3]) -  Z35 *I35;
V5 = diag (V5)
Angle5 = atan2(imag(V5),real(V5))
Angle6=Angle1;
VG=Vg
V = [V1;V2;V3;V4;V5;VG]
Vmag = abs(V)
Angle = [Angle1;Angle2;Angle3;Angle4;Angle5;Angle6];