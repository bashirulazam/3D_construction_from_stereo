function [W,R,T] = calc_W_R_T(P)
P3 = P(:,1:3);
P4 = P(:,4);
[A,B] = rq_mine(P3);
R = B;
W = A;
T = W\P4;
