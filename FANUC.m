%% FANUC R-2000iC/165F

syms q1 q2 q3 q4 q5 q6 d1 d2 d3 d4 d5 d6 real

%% FK

T = Rz(q1)*Tz(d1)*Ty(d2)*Rx(q2)*Tz(d3)*Rx(q3)*Tz(d4)*Ty(d5)*Ry(q4)*Rx(q5)*Ry(q6)*Ty(d6);
T = simplify(T);

disp('Forward kinematics');
disp(T);

%% Scew theory

% Origins
T0 = eye(4);
T1 = Rz(q1)*Tz(d1)*Ty(d2);
T2 = T1*Rx(q2)*Tz(d3);
T3 = T2*Rx(q3)*Tz(d4)*Ty(d5);
T4 = T3*Ry(q4);
T5 = T4*Rx(q5);
T6 = T5*Ry(q6)*Ty(d6);

T1 = simplify(T1);
T2 = simplify(T2);
T3 = simplify(T3);
T4 = simplify(T4);
T5 = simplify(T5);
T6 = simplify(T6);

O0 = T0(1:3,4);
O1 = T1(1:3,4);
O2 = T2(1:3,4);
O3 = T3(1:3,4);
O4 = O3;
O5 = O3;
O6 = T6(1:3,4);

% Rotation axis
Z0 = T0(1:3,3);
X1 = T1(1:3,1);
X2 = T2(1:3,1);
Y3 = T3(1:3,2);
X4 = T4(1:3,1);
Y5 = T5(1:3,2);

% Jacobians
J1 = [cross(Z0,O6-O0);Z0];
J2 = [cross(X1,O6-O1);X1];
J3 = [cross(X2,O6-O2);X2];
J4 = [cross(Y3,O6-O3);Y3];
J5 = [cross(X4,O6-O4);X4];
J6 = [cross(Y5,O6-O5);Y5];

% Jsc = [J1, J2, J3, J4, J5, J6];
Jsc = [simplify(J1), simplify(J2), simplify(J3), simplify(J4), simplify(J5), simplify(J6)];

disp('Jacobian obtained by Scew theory');
disp(Jsc);

%% Numerical method

% Inverse transpose
Ti = T;
Ti = inv(Ti(1:3,1:3));
Ti=[Ti,zeros(3,1);0 0 0 1];
T1 = simplify(Ti);

% Jacobians
Td = Rzd(q1)*Tz(d1)*Ty(d2)* Rx(q2)*Tz(d3)* Rx(q3)*Tz(d4)*Ty(d5)* Ry(q4)* Rx(q5)* Ry(q6)*Ty(d6)*Ti;
Td = simplify(Td);
J1 = Jcol(Td);
Td =  Rz(q1)*Tz(d1)*Ty(d2)*Rxd(q2)*Tz(d3)* Rx(q3)*Tz(d4)*Ty(d5)* Ry(q4)* Rx(q5)* Ry(q6)*Ty(d6)*Ti;
Td = simplify(Td);
J2 = Jcol(Td);
Td =  Rz(q1)*Tz(d1)*Ty(d2)* Rx(q2)*Tz(d3)*Rxd(q3)*Tz(d4)*Ty(d5)* Ry(q4)* Rx(q5)* Ry(q6)*Ty(d6)*Ti;
Td = simplify(Td);
J3 = Jcol(Td);
Td =  Rz(q1)*Tz(d1)*Ty(d2)* Rx(q2)*Tz(d3)* Rx(q3)*Tz(d4)*Ty(d5)*Ryd(q4)* Rx(q5)* Ry(q6)*Ty(d6)*Ti;
Td = simplify(Td);
J4 = Jcol(Td);
Td =  Rz(q1)*Tz(d1)*Ty(d2)* Rx(q2)*Tz(d3)* Rx(q3)*Tz(d4)*Ty(d5)* Ry(q4)*Rxd(q5)* Ry(q6)*Ty(d6)*Ti;
Td = simplify(Td);
J5 = Jcol(Td);
Td =  Rz(q1)*Tz(d1)*Ty(d2)* Rx(q2)*Tz(d3)* Rx(q3)*Tz(d4)*Ty(d5)* Ry(q4)* Rx(q5)*Ryd(q6)*Ty(d6)*Ti;
Td = simplify(Td);
J6 = Jcol(Td);

% Jnum = [J1, J2, J3, J4, J5, J6];
Jnum = [simplify(J1), simplify(J2), simplify(J3), simplify(J4), simplify(J5), simplify(J6)];

disp('Jacobian obtained by Numerical method');
disp(Jnum);

%% Comparing two solutions

disp('J numerical - J scew');
dif = Jnum - Jsc;
dif = simplify(dif);
disp(dif);

%% Kinematic Singularity analysis

q1 = 0;
q2 = 0;
q3 = 0;
q4 = 0;
q5 = 0;
q6 = 0;
d1 = 670;
d2 = 312;
d3 = 1075;
d4 = 225;
d5 = 1280;
d6 = 215;

J = Jsc;
J = subs(J, q1);
J = subs(J, q2);
J = subs(J, q3);
J = subs(J, q4);
J = subs(J, q5);
J = subs(J, q6);
J = subs(J, d1);
J = subs(J, d2);
J = subs(J, d3);
J = subs(J, d4);
J = subs(J, d5);
J = subs(J, d6);

disp('Jacobian with the given joint and link positions');
disp(J);

A = J'*J;
sigma = svd(A);

disp('Singular values');
disp(sigma);
disp('Multiplication of singular values');
disp(sigma(1)*sigma(2)*sigma(3)*sigma(4)*sigma(5)*sigma(6));