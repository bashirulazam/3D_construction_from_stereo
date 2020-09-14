function [a,b,c] = epiline(cl,rl,F)

Ul = [cl; rl ; 1];
right_param = F'*Ul;
a = right_param(1);
b = right_param(2);
c = right_param(3);