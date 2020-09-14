function [R,T,F] = form_F(Wl,Wr,Rl,Rr,Tl,Tr)

R = Rl*Rr';
T = Tl - R*Tr;
S = [ 0 -T(3) T(2);
    T(3)  0   -T(1);
    -T(2) T(1) 0];
E = S'*R;
F = inv(Wl)'*E*inv(Wr);