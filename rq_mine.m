function [R,Q] = rq_mine(A)

P = flipud(eye(size(A)));
A_reversed = P*A;
[Q_reversed,R_reversed] = qr(A_reversed');
Q = P*Q_reversed';
R = P*R_reversed'*P;