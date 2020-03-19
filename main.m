clear all; close all; clc;
import casadi.*

% Parameters of our walker (masses, lengths inertias, ...)
m1 = 3.2; m5 = 3.2;
m2 = 6.8; m4 = 6.8;
m3 = 20;
I1 = 0.93; I5 = 0.93;
I2 = 1.08; I4 = 1.08;
I3 = 2.22;
l1 = 0.4; l5 = 0.4;
l2 = 0.4; l4 = 0.4;
l3 = 0.625;
lc1 = l1 - 0.128; lc5 = lc1;
lc2 = l2-0.163; lc4 = lc2;
lc3 = 0.2;
g = 9.81;

% Declare CasADi variables that will be used to create CasADi function
q1_MX = MX.sym('q1_MX',1); q2_MX = MX.sym('q2_MX',1); q3_MX = MX.sym('q3_MX',1); q4_MX = MX.sym('q4_MX',1); q5_MX = MX.sym('q5_MX',1);
dq1_MX = MX.sym('dq1_MX',1); dq2_MX = MX.sym('dq2_MX',1); dq3_MX = MX.sym('dq3_MX',1); dq4_MX = MX.sym('dq4_MX',1); dq5_MX = MX.sym('dq5_MX',1);
ddq1_MX = MX.sym('ddq1_MX',1); ddq2_MX = MX.sym('ddq2_MX',1); ddq3_MX = MX.sym('ddq3_MX',1); ddq4_MX = MX.sym('ddq4_MX',1); ddq5_MX = MX.sym('ddq5_MX',1);
T1_MX = MX.sym('T1_MX',1); T2_MX = MX.sym('T2_MX',1); T3_MX = MX.sym('T3_MX',1); T4_MX = MX.sym('T4_MX',1); T5_MX = MX.sym('T5_MX',1);

% Generate a CasADi function for the implicit constraint in order to satisfy the system dynamics --> f(T,q,dq,ddq) == 0
eq_SysDyn_Error = eq_SysDyn(I1,I2,I3,I4,I5,T1_MX,T2_MX,T3_MX,T4_MX,T5_MX,ddq1_MX,ddq2_MX,ddq3_MX,ddq4_MX,ddq5_MX,dq1_MX,dq2_MX,dq3_MX,dq4_MX,dq5_MX,g,l1,l2,l4,lc1,lc2,lc3,lc4,lc5,m1,m2,m3,m4,m5,q1_MX,q2_MX,q3_MX,q4_MX,q5_MX);
f_eq_SysDyn_Error_Nominal = Function('f_eq_SysDyn_Error_Nominal',{T1_MX,T2_MX,T3_MX,T4_MX,T5_MX,ddq1_MX,ddq2_MX,ddq3_MX,ddq4_MX,ddq5_MX,dq1_MX,dq2_MX,dq3_MX,dq4_MX,dq5_MX,q1_MX,q2_MX,q3_MX,q4_MX,q5_MX},{eq_SysDyn_Error});

% Time horizon and mesh size for simulation
dt = 0.01;  % Mesh size
T = 0.8;    % Stride time
N = T/dt;   % Nr mesh intervals
time = 0:dt:T;

opti = casadi.Opti(); % Create opti instance

% Create optimization variables
q1 = opti.variable(1,N+1);   q2 = opti.variable(1,N+1);   q3 = opti.variable(1,N+1);   q4 = opti.variable(1,N+1); , q5 = opti.variable(1,N+1);
dq1 = opti.variable(1,N+1);  dq2 = opti.variable(1,N+1);  dq3 = opti.variable(1,N+1);  dq4 = opti.variable(1,N+1);  dq5 = opti.variable(1,N+1);
ddq1 = opti.variable(1,N);   ddq2 = opti.variable(1,N);   ddq3 = opti.variable(1,N);   ddq4 = opti.variable(1,N);   ddq5 = opti.variable(1,N);
T1 = opti.variable(1,N);     T2 = opti.variable(1,N);     T3 = opti.variable(1,N);     T4 = opti.variable(1,N);     T5 = opti.variable(1,N);

% Crude bounds on the segment orientations
opti.subject_to(-pi/2 < q1 < pi/2);
opti.subject_to(-pi/2 < q2 < pi/2);
opti.subject_to(-pi/3 < q3 < pi/3);
opti.subject_to(-pi/2 < q4 < pi/2);
opti.subject_to(-pi/2 < q5 < pi/2);

% Physiological joint limits
opti.subject_to(-pi < q1 - q2 < 0); % Knee joint limit (no hyperflexion)
opti.subject_to(-pi < q5 - q4 < 0); % Knee joint limit (no hyperflexion)


% Generate an initial guess for the positions (done very naively!)
q1_init = -pi/8; q2_init = pi/6; q3_init = -pi/6; q4_init = -pi/8; q5_init = -pi/6;
q1_final = -pi/6; q2_final = -pi/8; q3_final = 0; q4_final = pi/6; q5_final = -pi/8;
q1guess = linspace( q1_init, q1_final, N+1);
q2guess = linspace( q2_init, q2_final, N+1);
q3guess = linspace( q3_init, q3_final, N+1);
q4guess = linspace( q4_init, q4_final, N+1);
q5guess = linspace( q5_init, q5_final, N+1);

opti.set_initial(q1, q1guess);
opti.set_initial(q2, q2guess);
opti.set_initial(q3, q3guess);
opti.set_initial(q4, q4guess);
opti.set_initial(q5, q5guess);


% Generate heel-strike map (periodicity for swing vs stance leg + impulsive
% collision constraint)
% Impulsive collision :: joint configuration does not change at collision,
% joint velocities do change, in such a way that angular momentum is
% conserved around the collision point and all joints (giving 5 constraint equations).
q1_min = q1(:,end); q2_min = q2(:,end); q3_min = q3(:,end); q4_min = q4(:,end); q5_min = q5(:,end);
q1_plus = q1(:,1); q2_plus = q2(:,1); q3_plus = q3(:,1); q4_plus = q4(:,1); q5_plus = q5(:,1);
dq1_min = dq1(:,end); dq2_min = dq2(:,end); dq3_min = dq3(:,end); dq4_min = dq4(:,end); dq5_min = dq5(:,end);
dq1_plus = dq1(:,1); dq2_plus = dq2(:,1); dq3_plus = dq3(:,1); dq4_plus = dq4(:,1); dq5_plus = dq5(:,1);
heelStrike_error = eq_HeelStrike(I1,I2,I3,I4,I5,dq1_min,dq2_min,dq3_min,dq4_min,dq5_min,dq1_plus,dq2_plus,dq3_plus,dq4_plus,dq5_plus,l1,l2,l4,l5,lc1,lc2,lc3,lc4,lc5,m1,m2,m3,m4,m5,q1_min,q2_min,q3_min,q4_min,q5_min,q1_plus,q2_plus,q3_plus,q4_plus,q5_plus)
opti.subject_to(heelStrike_error == 0);


J = 0;
for k=1:N
    % State at mesh point k
    q1k = q1(:,k);     q2k = q2(:,k);     q3k = q3(:,k);     q4k = q4(:,k);     q5k = q5(:,k);
    dq1k = dq1(:,k);   dq2k = dq2(:,k);   dq3k = dq3(:,k);   dq4k = dq4(:,k);   dq5k = dq5(:,k);
    
    % Control/lifted state variable for mesh k
    ddq1k = ddq1(:,k); ddq2k = ddq2(:,k); ddq3k = ddq3(:,k); ddq4k = ddq4(:,k); ddq5k = ddq5(:,k);
    T1k = T1(:,k);     T2k = T2(:,k);     T3k = T3(:,k);     T4k = T4(:,k);     T5k = T5(:,k);
    
    % State at mesh point k+1
    q1k_plus = q1(:,k+1);     q2k_plus = q2(:,k+1);     q3k_plus = q3(:,k+1);     q4k_plus = q4(:,k+1);     q5k_plus = q5(:,k+1);
    dq1k_plus = dq1(:,k+1);   dq2k_plus = dq2(:,k+1);   dq3k_plus = dq3(:,k+1);   dq4k_plus = dq4(:,k+1);   dq5k_plus = dq5(:,k+1);
       
    % Collect state at k and k+1
    Xk = [q1k; q2k; q3k; q4k; q5k; dq1k; dq2k; dq3k; dq4k; dq5k];
    Xk_next = [q1k_plus; q2k_plus; q3k_plus; q4k_plus; q5k_plus; dq1k_plus; dq2k_plus; dq3k_plus; dq4k_plus; dq5k_plus];
    
    % Collect state derivative (we use backward (not forward) Euler!)
    Uk = [dq1k_plus; dq2k_plus; dq3k_plus; dq4k_plus; dq5k_plus; ddq1k; ddq2k; ddq3k; ddq4k; ddq5k];
    
    % Integration
    opti.subject_to(eulerIntegrator(Xk,Xk_next,Uk,dt) == 0);
       
    % Dynamics error (backward Euler - derivative at k+1)
    error = f_eq_SysDyn_Error_Nominal(T1k,T2k,T3k,T4k,T5k,ddq1k,ddq2k,ddq3k,ddq4k,ddq5k,dq1k_plus,dq2k_plus,dq3k_plus,dq4k_plus,dq5k_plus,q1k_plus,q2k_plus,q3k_plus,q4k_plus,q5k_plus);
    opti.subject_to(error == 0);
    
    % Cost function contributions
    % Main term is torque minimization, we add minimization of
    % accelerations for regularization (contribution is confirmed to be
    % small)
    J = J + (T1k.^2 + T2k.^2 + T3k.^2 + T4k.^2 + T5k.^2)*dt + 1e-1*(ddq1k.^2 + ddq2k.^2 + ddq3k.^2 + ddq4k.^2 + ddq5k.^2)*dt;
    
    % Joint locations in x-y plane
    P_J = JointPos(l1,l2,l3,l4,l5,q1k,q2k,q3k,q4k,q5k);
    
    % Impose that swing foot does not penetrate ground
    opti.subject_to(P_J(10) > -1e-4);
    
    
    % OPTIONAL CONSTRAINTS TO CHANGE GAIT PATTERN
    
    % Avoid with swing foot a 45cm circle around stance foot
%     opti.subject_to(P_J(9)^2 + P_J(10)^2 > 0.45^2);
    
    % Impose crouching by keeping pelvis below 0.6m
%     opti.subject_to(P_J(4)<0.6);

    % Walk with 'pin-boots' (no ankle muscles)
%     opti.subject_to(T1k  == 0.0);

    % Walk without muscles in stance knee
%     opti.subject_to(T2k  == 0.0);

    % Walk without hip muscles (This is a nice one :-)!!)
%     opti.subject_to(T3k  == 0.0);
%     opti.subject_to(T4k  == 0.0);

    % Walk with only ankle muscles
%     opti.subject_to(T2k  == 0.0);
%     opti.subject_to(T3k  == 0.0);
%     opti.subject_to(T4k  == 0.0);
%     opti.subject_to(T5k  == 0.0);

    % Walk with only knee muscles
%     opti.subject_to(T1k  == 0.0);
%     opti.subject_to(T3k  == 0.0);
%     opti.subject_to(T4k  == 0.0);
    
    % Walk with only hip muscles
%     opti.subject_to(T1k  == 0.0);
%     opti.subject_to(T2k  == 0.0);
%     opti.subject_to(T5k  == 0.0);
end


% Periodicity constraint (part is in the heel strike map)
JointPosFinal = JointPos(l1,l2,l3,l4,l5,q1(:,end),q2(:,end),q3(:,end),q4(:,end),q5(:,end));
opti.subject_to(JointPosFinal(9:10,:) == [0.5;0]); %  Impose step length of 1 meter

% Impose that simulation start and end is at 'toe-off' and 'heel-strike'.
% We do this by saying that the swing foot should have a positive
% y-velocity at the start (toe-off) and negative at the end (heel strike).
JointVelInit = JointVel(dq1(:,1),dq2(:,1),dq3(:,1),dq4(:,1),dq5(:,1),l1,l2,l3,l4,l5,q1(:,1),q2(:,1),q3(:,1),q4(:,1),q5(:,1));
JointVelFinal = JointVel(dq1(:,end),dq2(:,end),dq3(:,end),dq4(:,end),dq5(:,end),l1,l2,l3,l4,l5,q1(:,end),q2(:,end),q3(:,end),q4(:,end),q5(:,end));
opti.subject_to(JointVelInit(10) > 0);
opti.subject_to(JointVelFinal(10) < 0);

% Define cost
opti.minimize(J);

% Create an NLP solver
optionssol.ipopt.linear_solver = 'mumps';
optionssol.ipopt.tol = 1e-2;
optionssol.ipopt.constr_viol_tol = 1e-10;
optionssol.ipopt.dual_inf_tol = 1e-8;
optionssol.ipopt.compl_inf_tol = 1e-6;
optionssol.ipopt.max_iter = 10000;
% optionssol.ipopt.hessian_approximation = 'limited-memory';
optionssol.ipopt.mu_strategy = 'adaptive';
opti.solver('ipopt',optionssol);

diary('2D-gait_diary.txt');
sol = opti.solve();
diary off

% Plot the torques, positions and velocities
T1_sol = sol.value(T1);
T2_sol = sol.value(T2);
T3_sol = sol.value(T3);
T4_sol = sol.value(T4);
T5_sol = sol.value(T5);
figure(1)
plot(time(1:end-1),T1_sol,time(1:end-1),T2_sol,time(1:end-1),T3_sol,time(1:end-1),T4_sol,time(1:end-1),T5_sol);
legend('stance ankle','stance knee','stance hip','swing hip','swing knee');
xlabel(['time [s]'])
ylabel(['Torque [Nm]'])

q1_sol = sol.value(q1);
q2_sol = sol.value(q2);
q3_sol = sol.value(q3);
q4_sol = sol.value(q4);
q5_sol = sol.value(q5);
figure(2)
plot(time,180/pi*q1_sol,time,180/pi*q2_sol,time,180/pi*q3_sol,time,180/pi*q4_sol,time,180/pi*q5_sol);
legend('stance tib','stance fem','trunk','swing fem','swing tib');
xlabel(['time [s]'])
ylabel(['segment pos [°]'])

dq1_sol = sol.value(dq1);
dq2_sol = sol.value(dq2);
dq3_sol = sol.value(dq3);
dq4_sol = sol.value(dq4);
dq5_sol = sol.value(dq5);
figure(3)
plot(time,180/pi*dq1_sol,time,180/pi*dq2_sol,time,180/pi*dq3_sol,time,180/pi*dq4_sol,time,180/pi*dq5_sol);
legend('stance tib','stance fem','trunk','swing fem','swing tib');
xlabel(['time'])
ylabel(['segment vel [°]'])

ddq1_sol = sol.value(ddq1);
ddq2_sol = sol.value(ddq2);
ddq3_sol = sol.value(ddq3);
ddq4_sol = sol.value(ddq4);
ddq5_sol = sol.value(ddq5);

% Reconstruct cost function
J_torque = (sumsqr(T1_sol) + sumsqr(T2_sol) + sumsqr(T3_sol) + sumsqr(T4_sol) + sumsqr(T5_sol))*dt;
J_ddq = 1e-1*(sumsqr(ddq1_sol) + sumsqr(ddq2_sol) + sumsqr(ddq3_sol) + sumsqr(ddq4_sol) + sumsqr(ddq5_sol))*dt;

% Animate the gait pattern
P_J = JointPos(l1,l2,l3,l4,l5,q1_sol,q2_sol,q3_sol,q4_sol,q5_sol)';
animation(P_J,dt);