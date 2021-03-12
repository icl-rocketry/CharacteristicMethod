clear
clc
close all

%% solve for limits of P and Q as described in sheffield paper

syms V1 V2 V3 th1 th2 th3 mu1 mu2 mu3 x1 x2 x3 r1 r2 r3

eq1 = V3 == 1/(cot(mu1)/V1+cot(mu2)/V2)*(cot(mu1)*(1+th3/r3*(x3-x1))+cot(mu2)*(1+th3/r3*(x3-x2))+th2-th1); %implicit equation for velocity for first iteration near centreline
eq2 = th3 == th1 + cot(mu1)*((V3-V1)/V1)-th3/r3*(x3-x1); %implicit equation for theta for first iteration near centreline

eqs = [eq1 ;eq2];

S = solve(eqs,[V3 th3],'ReturnConditions',true); %solutions of the two equations

Veq = S.V3; %velocity for first iteration near centreline
thetaeq = S.th3; %theta for first iteration near centreline

syms V3 th3
syms P1 P2 P3 Q1 Q2 Q3 r4 th4

Veq = simplify(Veq,'All',true,'IgnoreAnalyticConstraints',true);
thetaeq = simplify(thetaeq,'All',true,'IgnoreAnalyticConstraints',true);

%% solve for alternative MOC scheme

syms r1 r2 r3 x1 x2 x3 th1 th2 th3 mu1 mu2 mu3 V1 V2 V3 

eq1 = ((r3 - r1)/(x3-x1)) == tan(th1+mu1);
eq2 = ((r2-r3)/(x2-x3)) == tan(th2-mu2);

S = solve ([eq1;eq2],[x3 r3]);

x3eq = S.x3;

S2 = solve(eq1,r3);

r3eq = S2;

%

eq3 = (V2-V3)/V3 + (th2 - th3)*tan(mu2) - (tan(mu2)*sin(mu2)*sin(th2)*(x2-x3))/(cos(th2-mu2)*r3) ==0;
eq4 = (V3-V1)/V1 - (th3 - th1)*tan(mu1) - (tan(mu1)*sin(mu1)*sin(th1)*(x3-x1))/(cos(th1+mu1)*r1) ==0;

S = solve([eq3;eq4],[V3 th3]);

V3eq = S.V3(1);

S2 = solve(eq3,th3);

th3eq = S2;

% 
% Pn = tan(mu1)*sin(mu1)*sin(th1)/(r1*cos(th1+mu1));
% Qn = tan(mu2)*sin(mu2)*sin(th2)/(r2*cos(th2+mu2));
% 
% simplify(V3eq)

P = tan(mu1)*sin(mu1)*sin(th1)/(r3*cos(th1+mu1));



