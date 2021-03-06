clc;
clear;


Vgams=[0.97647999
0.97647999
0.97647999
0.97014405
0.97014405
0.97014405
0.98088792
0.98088792
0.98088792
0.95000000
0.95000000
0.95000000
1.05000000
1.05000000
1.05000000
1.00549257
1.00549257
1.00549257
];
Vsdp=[0.97648001
0.97648001
0.97648001
0.97014407
0.97014407
0.97014407
0.98088793
0.98088793
0.98088793
0.95000000
0.95000000
0.95000000
1.05000000
1.05000000
1.05000000
1.00558286
1.00558286
1.00558286
];
Vsdp_ref=[0.97656768
0.97656768
0.97656768
0.97023233
0.97023233
0.97023233
0.98096602
0.98096602
0.98096602
0.95009024
0.95009024
0.95009024
1.05000000
1.05000000
1.05000000
1.00558286
1.00558286
1.00558286
];

Agams=[0.00000000;
    -2.09439510;
    2.09439510;
    -0.00572333;
    -2.10011844;
    2.08867177;
    0.01205493;
    -2.08234018;
    2.10645003;
    -0.03107859;
    -2.12547369;
    2.06331651;
    0.14813822;
    -1.94625688;
    2.24253332;
    0.00000000;
    -2.09439510;
    2.09439510;
    ];

Asdp =[            0
  -2.094395102393195
   2.094395102393195
  -0.005723331870225
  -2.100118434263298
   2.088671770522867
   0.012054932422535
  -2.082340169971222
   2.106450034815151
  -0.031078590689566
  -2.125473693082874
   2.063316511703699
   0.148138268780510
  -1.946256833618453
   2.242533371168783
                   0
  -2.094395102393195
   2.094395102393195];

Asdp_ref=[0
  -2.094395102393196
   2.094395102393195
  -0.005722300299969
  -2.100117402693194
   2.088672802093210
   0.012060366702704
  -2.082334735690521
   2.106455469095953
  -0.031072843472301
  -2.125467945865730
   2.063322258920798
   0.148176481064039
  -1.946218621328940
   2.242571583457952
                   0
  -2.094395102393196
   2.094395102393195];
figure(1)
subplot(2,1,1) 
plot(1:18,Vgams,'-*',1:18,Vsdp,'-o',1:18,Vsdp_ref,'-+')
xlabel('Bus');
ylabel('Voltage Magnitude/p.u.');
xticks([2 5 8 11 14 17])
set(gca, 'xticklabel',{'V1';'V2';'V3';'V4';'V5';'V6'});
legend('V-nonconvex','V-sdp','V-sdp-ref');

subplot(2,1,2) 

plot(1:18,Agams,'-*',1:18,Asdp,'-o',1:18,Asdp_ref,'-+')
xlabel('Bus');
ylabel('Phase Angle/rad');
xticks([2 5 8 11 14 17])
 set(gca, 'xticklabel',{'V1';'V2';'V3';'V4';'V5';'V6'});
legend('\theta-nonconvex','\theta-sdp','\theta-sdp-ref');

Pg_gams=[0;-6.6613E-16;0;3.55060046;3.55060046;3.55060046];
Qg_gams=[0.98601274;0.98601274;0.98601274;0.01862115;0.01862115;0.01862115];
Pg_sdp=[0;0;0;3.55060064;3.55060064;3.55060064];
Qg_sdp=[0.98601387;0.98601387;0.98601387;0.01862016;0.01862016;0.01862016];
Pg_sdp_ref=[0;0;0;3.55060213;3.55060213;3.55060213];
Qg_sdp_ref=[0.98795703;0.98795703;0.98795703;0.01667742;0.01667742;0.01667742];

Pg_sdp_MAE=max(abs(Pg_sdp-Pg_gams))
Qg_sdp_MAE=max(abs(Qg_sdp-Qg_gams))
Pg_sdp_RMSE=sqrt((sum(abs(Pg_sdp-Pg_gams).^2))/15)
Qg_sdp_RMSE=sqrt((sum(abs(Qg_sdp-Qg_gams).^2))/15)
V_sdp_MAE=max(abs(Vsdp-Vgams))
V_sdp_RMSE=sqrt((sum(abs(Vsdp-Vgams).^2))/15)

Pg_sdp_ref_MAE=max(abs(Pg_sdp_ref-Pg_gams))
Qg_sdp_ref_MAE=max(abs(Qg_sdp_ref-Qg_gams))
Pg_sdp_ref_RMSE=sqrt((sum(abs(Pg_sdp_ref-Pg_gams).^2))/15)
Qg_sdp_ref_RMSE=sqrt((sum(abs(Qg_sdp_ref-Qg_gams).^2))/15)
V_sdp_ref_MAE=max(abs(Vsdp_ref-Vgams))
V_sdp_ref_RMSE=sqrt((sum(abs(Vsdp_ref-Vgams).^2))/15)

figure(2)
S_rmse=[sqrt((sum(abs(Pg_sdp-Pg_gams).^2))/6);sqrt((sum(abs(Pg_sdp_ref-Pg_gams).^2))/6);sqrt((sum(abs(Qg_sdp-Qg_gams).^2))/6);sqrt((sum(abs(Qg_sdp_ref-Qg_gams).^2))/6)];
bar(S_rmse)
xlabel('DG Bus');
ylabel('Power RMSE');
xticks([1 2 3 4])
set(gca, 'xticklabel',{'Pdg1-sdp';'Pdg5-sdp-ref';'Qdg1-sdp';'Qdg5-sdp-ref'});

