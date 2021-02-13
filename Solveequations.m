clear
clc
close all

syms V1 V2 V3 th1 th2 th3 mu1 mu2 mu3 x1 x2 x3 r1 r2 r3

eq1 = V3 == 1/(cot(mu1)/V1+cot(mu2)/V2)*(cot(mu1)*(1+th3/r3*(x3-x1))+cot(mu2)*(1+th3/r3*(x3-x2))+th2-th1); %implicit equation for velocity for first iteration near centreline
eq2 = th3 == th1 + cot(mu1)*((V3-V1)/V1)-th3/r3*(x3-x1); %implicit equation for theta for first iteration near centreline

eqs = [eq1 ;eq2];

S = solve(eqs,[V3 th3],'ReturnConditions',true); %solutions of the two equations

Veq = S.V3; %velocity for first iteration near centreline
thetaeq = S.th3; %theta for first iteration near centreline

syms V3 th3
syms P1 P2 P3 Q1 Q2 Q3 r4 th4

simplify(Veq,'All',true,'IgnoreAnalyticConstraints',true)
simplify(thetaeq,'All',true,'IgnoreAnalyticConstraints',true)