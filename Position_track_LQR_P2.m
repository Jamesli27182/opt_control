clc;close all;clear all
m = 1;k = 1;b = 0.5;
%setup state-space-c
A = [0 1;-k/m -b/m];
B = [0;1/m];
C = [1 0];
D = 0;
n = size(A,1);
%%setup lqr========================
R= 0.1;S = [1 0;0 1];Q = [1 0;0 1];
d_T = 0.1;
ss_d = c2d(ss(A,B,C,D),d_T);
%%setup state-space-d
A = ss_d.A;
B = ss_d.B;
C = ss_d.C;
D = ss_d.D;
%Inital set
x_0 = [0;0];
u_0 = 1;
u = u_0 ;
%Goal position
x_d = [1;0];
x = x_0;
Xa = [x;x_d];
k_steps = 200+1;
N = k_steps;
x_record = zeros(n(1),k_steps);
t = 0:d_T:(k_steps-1)*d_T;
u_record = zeros(1,k_steps);
u(1)=u;
x_record(:,1)=x_0;
u_record(1) = u_0;
[p,~] = size(B');
%% Augmented Matrix
Aa = [A,zeros(n);zeros(n),eye(n)];
Ba = [B;0;0];
Ca = [eye(n) -eye(n)];
Qa = Ca'*Q*Ca;
Sa = Ca'*S*Ca;

%% lqr Control
F_N=LQR_Gain(Aa,Ba,Qa,R,Sa,N)
F = F_N(1,:)
x = x_0;
for k = 2:k_steps;
   u = -F*Xa;
   x = A*x+B*u;
   %x_d = x_d + [0.5;0.5*sin(k)]
   Xa = [x;x_d];
   x_record(:,k)=x;
   u_record(:,k)=u;
end
subplot(2,1,1);
plot(t,x_record);
legend('x','v');
subplot(2,1,2);
plot(t,u_record);
legend('input u');
%% local function
function [F_N] = LQR_Gain(A,B,Q,R,S,N)
P_k = S;
for k = 1:N;
F = inv(B'*P_k*B+R)*B'*P_k*A;
P_k = (A-B*F)'*P_k*(A-B*F)+(F)'*R*F+Q;
if k == 1;
    F_N =F;
else
    F_N = [F;F_N];
end
end
end