%% EoM quadruple Pendulum


clc; clear;
disp('Calculating the equations of motion symbolically...')
syms q1 q2 q3 q4 q5 dq1 dq2 dq3 dq4 dq5 ddq1 ddq2 ddq3 ddq4 ddq5 x dx ddx y dy ddy T1 T2 T3 T4 T5  'real'; 
syms q1_plus q2_plus q3_plus q4_plus q5_plus dq1_plus dq2_plus dq3_plus dq4_plus dq5_plus
syms q1_min q2_min q3_min q4_min q5_min dq1_min dq2_min dq3_min dq4_min dq5_min

% th_A == hook of angle to world axisses
% th_K == 
% th_H == 
% dth_A == first time derivative of th_A
% dth_K == first time derivative of th_K
% ddth_A == second time derivative of th_A
% ddth_K == second time derivative of th_K

syms m1 m2 m3 m4 m5 g l1 l2 l3 l4 l5 I1 I2 I3 I4 I5 lc1 lc2 lc3 lc4 lc5 'real';
% m1, m2, m3 == mass of the pendulum links
% g == acceleration due to gravity
% l1, l2, l3== length of the pendulum's first and second links
% I1, I2, I3 == moment of inertias of the first and second links
% lc1, lc2, lc3 == distance between link's parent joint and CoM


% Note on the ANGLE DEFINITIONS:
% The angles are the positive angles with respect to the direction of the
% previous segment. For the first angle (th_A) - this is the positive angle
% with respect to the y-axis.

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%                            vector stuff                                 %
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
i = [1; 0; 0];
j = [0; 1; 0];
k = [0; 0; 1];

% Positions of joints and link COMs
P_0 = [0;0;0];

G_1 = P_0 + lc1*[cos(q1 + pi/2); sin(q1 + pi/2);0];
P_1 = P_0 + l1*[cos(q1 + pi/2); sin(q1 + pi/2);0];

G_2 = P_1 + lc2*[cos(q2 + pi/2); sin(q2 + pi/2);0];
P_2 = P_1 + l2*[cos(q2 + pi/2); sin(q2 + pi/2);0];

G_3 = P_2 + lc3*[cos(q3 + pi/2); sin(q3 + pi/2);0];
P_3 = P_2 + l3*[cos(q3 + pi/2); sin(q3 + pi/2);0];

G_4 = P_2 - lc4*[cos(q4 + pi/2); sin(q4 + pi/2);0];
P_4 = P_2 - l4*[cos(q4 + pi/2); sin(q4 + pi/2);0];

G_5 = P_4 - lc5*[cos(q5 + pi/2); sin(q5 + pi/2);0];
P_5 = P_4 - l5*[cos(q5 + pi/2); sin(q5 + pi/2);0];

% Collect in column vectors
P = [P_0; P_1; P_2; P_3; P_4; P_5];
G = [G_1; G_2; G_3; G_4; G_5];
q = [q1; q2; q3; q4; q5]; q_plus = [q1_plus; q2_plus; q3_plus; q4_plus; q5_plus]; q_min = [q1_min; q2_min; q3_min; q4_min; q5_min]; 
dq = [dq1; dq2; dq3; dq4; dq5]; dq_plus = [dq1_plus; dq2_plus; dq3_plus; dq4_plus; dq5_plus]; dq_min = [dq1_min; dq2_min; dq3_min; dq4_min; dq5_min]; 
ddq = [ddq1; ddq2; ddq3; ddq4; ddq5];

% Get velocities
dP_1 = jacobian(P_1,q)*dq; dP_2 = jacobian(P_2,q)*dq; dP_3 = jacobian(P_3,q)*dq; dP_4 = jacobian(P_4,q)*dq; dP_5 = jacobian(P_5,q)*dq;
dG_1 = jacobian(G_1,q)*dq; dG_2 = jacobian(G_2,q)*dq; dG_3 = jacobian(G_3,q)*dq; dG_4 = jacobian(G_4,q)*dq; dG_5 = jacobian(G_5,q)*dq; 

% Get accelerations
ddP_1 = jacobian(dP_1,[q;dq])*[dq; ddq]; ddP_2 = jacobian(dP_2,[q;dq])*[dq; ddq]; ddP_3 = jacobian(dP_3,[q;dq])*[dq; ddq]; ddP_4 = jacobian(dP_4,[q;dq])*[dq; ddq]; ddP_5 = jacobian(dP_5,[q;dq])*[dq; ddq];
ddG_1 = jacobian(dG_1,[q;dq])*[dq; ddq]; ddG_2 = jacobian(dG_2,[q;dq])*[dq; ddq]; ddG_3 = jacobian(dG_3,[q;dq])*[dq; ddq]; ddG_4 = jacobian(dG_4,[q;dq])*[dq; ddq]; ddG_5 = jacobian(dG_5,[q;dq])*[dq; ddq];


% EQUATIONS OF MOTION
% Write equations for angular momentum balance about each joint o get five
% independent equations in order to construct EoM

% 1.
eq1 = [0; 0; T1] + k.*(cross((G_1 - P_0),(-m1*g*j)) + cross((G_2 - P_0),(-m2*g*j)) + cross((G_3 - P_0),(-m3*g*j)) + cross((G_4 - P_0),(-m4*g*j)) + cross((G_5 - P_0),(-m5*g*j))) - ...
    k.*(cross((G_1 - P_0),(m1*ddG_1)) + cross((G_2 - P_0),(m2*ddG_2)) + cross((G_3 - P_0),(m3*ddG_3)) + cross((G_4 - P_0),(m4*ddG_4)) + cross((G_5 - P_0),(m5*ddG_5)) + ...
    ddq1*I1*k + ddq2*I2*k + ddq3*I3*k + ddq4*I4*k + ddq5*I5*k);
eq1 = eq1(3);


% 2.
eq2 = [0; 0; T2] + k.*(cross((G_2 - P_1),(-m2*g*j)) + cross((G_3 - P_1),(-m3*g*j)) + cross((G_4 - P_1),(-m4*g*j)) + cross((G_5 - P_1),(-m5*g*j))) - ...
    k.*(cross((G_2 - P_1),(m2*ddG_2)) + cross((G_3 - P_1),(m3*ddG_3)) + cross((G_4 - P_1),(m4*ddG_4)) + cross((G_5 - P_1),(m5*ddG_5)) + ...
    ddq2*I2*k + ddq3*I3*k + ddq4*I4*k + ddq5*I5*k);
eq2 = eq2(3);

% 3.
eq3 = [0; 0; T3] + k.*(cross((G_3 - P_2),(-m3*g*j)) + cross((G_4 - P_2),(-m4*g*j)) + cross((G_5 - P_2),(-m5*g*j))) - ...
    k.*(cross((G_3 - P_2),(m3*ddG_3)) + cross((G_4 - P_2),(m4*ddG_4)) + cross((G_5 - P_2),(m5*ddG_5)) + ...
    ddq3*I3*k + ddq4*I4*k + ddq5*I5*k);
eq3 = eq3(3);

% 4.
eq4 = [0; 0; T4] + k.*(cross((G_4 - P_2),(-m4*g*j)) + cross((G_5 - P_2),(-m5*g*j))) - ...
    k.*(cross((G_4 - P_2),(m4*ddG_4)) + cross((G_5 - P_2),(m5*ddG_5)) + ...
    ddq4*I4*k + ddq5*I5*k);
eq4 = eq4(3);

% 5.
eq5 = [0; 0; T5] + k.*(cross((G_5 - P_4),(-m5*g*j))) - ...
    k.*(cross((G_5 - P_4),(m5*ddG_5)) + ...
    ddq5*I5*k);
eq5 = eq5(3);

eq_SysDyn = simplify([eq1; eq2; eq3; eq4; eq5]);
eq_SysDyn_fcn = matlabFunction(eq_SysDyn,'File','eq_SysDyn.m');

% HEEL STRIKE MAP

% Switch positions map
eq_h1 = q1_plus - q5_min;
eq_h2 = q2_plus - q4_min;
eq_h3 = q3_plus - q3_min;
eq_h4 = q4_plus - q2_min;
eq_h5 = q5_plus - q1_min;

% Impulsive collision equations (conserve angular momentum around collision
% point and all joints
G_1_plus = subs(G_1,q,q_plus); G_2_plus = subs(G_2,q,q_plus); G_3_plus = subs(G_3,q,q_plus); G_4_plus = subs(G_4,q,q_plus); G_5_plus = subs(G_5,q,q_plus); 
dG_1_plus = subs(dG_1,[q; dq],[q_plus; dq_plus]); dG_2_plus = subs(dG_2,[q; dq],[q_plus; dq_plus]); dG_3_plus = subs(dG_3,[q; dq],[q_plus; dq_plus]); dG_4_plus = subs(dG_4,[q; dq],[q_plus; dq_plus]); dG_5_plus = subs(dG_5,[q; dq],[q_plus; dq_plus]); 
G_1_min = subs(G_1,q,q_min); G_2_min = subs(G_2,q,q_min); G_3_min = subs(G_3,q,q_min); G_4_min = subs(G_4,q,q_min); G_5_min = subs(G_5,q,q_min); 
dG_1_min = subs(dG_1,[q; dq],[q_min; dq_min]); dG_2_min = subs(dG_2,[q; dq],[q_min; dq_min]); dG_3_min = subs(dG_3,[q; dq],[q_min; dq_min]); dG_4_min = subs(dG_4,[q; dq],[q_min; dq_min]); dG_5_min = subs(dG_5,[q; dq],[q_min; dq_min]); 

P_0_plus = [0; 0; 0]; P_1_plus = subs(P_1,q,q_plus); P_2_plus = subs(P_2,q,q_plus); P_3_plus = subs(P_3,q,q_plus); P_4_plus = subs(P_4,q,q_plus); P_5_plus = subs(P_5,q,q_plus); 
P_1_min = subs(P_1,q,q_min); P_2_min = subs(P_2,q,q_min); P_3_min = subs(P_3,q,q_min); P_4_min = subs(P_4,q,q_min); P_5_min = subs(P_5,q,q_min); 

% 6.
eq_h6 = k.*( cross((G_1_min-P_5_min),m1*dG_1_min)+dq1_min*I1*k + cross((G_2_min-P_5_min),m2*dG_2_min)+dq2_min*I2*k + cross((G_3_min-P_5_min),m3*dG_3_min)+dq3_min*I3*k + cross((G_4_min-P_5_min),m4*dG_4_min)+dq4_min*I4*k + cross((G_5_min-P_5_min),m5*dG_5_min)+dq5_min*I5*k) - ...
        k.*( cross((G_1_plus-P_0_plus),m1*dG_1_plus)+dq1_plus*I1*k + cross((G_2_plus-P_0_plus),m2*dG_2_plus)+dq2_plus*I2*k + cross((G_3_plus-P_0_plus),m3*dG_3_plus)+dq3_plus*I3*k + cross((G_4_plus-P_0_plus),m4*dG_4_plus)+dq4_plus*I4*k + cross((G_5_plus-P_0_plus),m5*dG_5_plus)+dq5_plus*I5*k)
eq_h6 = simplify(eq_h6(3));

% 7.
eq_h7 = k.*( cross((G_1_min-P_4_min),m1*dG_1_min)+dq1_min*I1*k + cross((G_2_min-P_4_min),m2*dG_2_min)+dq2_min*I2*k + cross((G_3_min-P_4_min),m3*dG_3_min)+dq3_min*I3*k + cross((G_4_min-P_4_min),m4*dG_4_min)+dq4_min*I4*k) - ...
        k.*( cross((G_2_plus-P_1_plus),m2*dG_2_plus)+dq2_plus*I2*k + cross((G_3_plus-P_1_plus),m3*dG_3_plus)+dq3_plus*I3*k + cross((G_4_plus-P_1_plus),m4*dG_4_plus)+dq4_plus*I4*k + cross((G_5_plus-P_1_plus),m5*dG_5_plus)+dq5_plus*I5*k)
eq_h7 = simplify(eq_h7(3));

% 8.
eq_h8 = k.*( cross((G_1_min-P_2_min),m1*dG_1_min)+dq1_min*I1*k + cross((G_2_min-P_2_min),m2*dG_2_min)+dq2_min*I2*k + cross((G_3_min-P_2_min),m3*dG_3_min)+dq3_min*I3*k) - ...
        k.*( cross((G_3_plus-P_2_plus),m3*dG_3_plus)+dq3_plus*I3*k + cross((G_4_plus-P_2_plus),m4*dG_4_plus)+dq4_plus*I4*k + cross((G_5_plus-P_2_plus),m5*dG_5_plus)+dq5_plus*I5*k)
eq_h8 = simplify(eq_h8(3));

% 9.
eq_h9 = k.*( cross((G_1_min-P_2_min),m1*dG_1_min)+dq1_min*I1*k + cross((G_2_min-P_2_min),m2*dG_2_min)+dq2_min*I2*k) - ...
        k.*( cross((G_4_plus-P_2_plus),m4*dG_4_plus)+dq4_plus*I4*k + cross((G_5_plus-P_2_plus),m5*dG_5_plus)+dq5_plus*I5*k)
eq_h9 = simplify(eq_h9(3));

% 10.
eq_h10 = k.*( cross((G_1_min-P_1_min),m1*dG_1_min)+dq1_min*I1*k) - ...
         k.*( cross((G_5_plus-P_4_plus),m5*dG_5_plus)+dq5_plus*I5*k)
eq_h10 = simplify(eq_h10(3));

eq_h = [eq_h1; eq_h2; eq_h3; eq_h4; eq_h5; eq_h6; eq_h7; eq_h8; eq_h9; eq_h10];
eq_h_fcn = matlabFunction(eq_h,'File','eq_HeelStrike.m');

JointPos = P([4 5 7 8 10 11 13 14 16 17],:);
JointVel = [dP_1(1:2); dP_2(1:2); dP_3(1:2); dP_4(1:2); dP_5(1:2)];

P_fcn = matlabFunction(JointPos,'File','JointPos.m');
dP_fcn = matlabFunction(JointVel,'File','JointVel.m');

% Get relative joint angles
q_ANK = q1;
q_stanceKNEE = q1 - q2;
q_stanceHIP = q2 - q3;
q_swingHIP = q4 - q3;
q_swingKNEE = q5 - q4;
relativeJointPos = [q_ANK; q_stanceKNEE; q_stanceHIP; q_swingHIP; q_swingKNEE];
Prel_fcn = matlabFunction(relativeJointPos,'File','relativeJointPos.m');

% Get relative joint velocities
dq_ANK = dq1;
dq_stanceKNEE = dq1 - dq2;
dq_stanceHIP = dq2 - dq3;
dq_swingHIP = dq4 - dq3;
dq_swingKNEE = dq5 - dq4;
relativeJointVel = [dq_ANK; dq_stanceKNEE; dq_stanceHIP; dq_swingHIP; dq_swingKNEE];
dPrel_fcn = matlabFunction(relativeJointVel,'File','relativeJointVel.m');

