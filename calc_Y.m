function [Y] = calc_Y(B,C,Z)

Y = -inv(B'*B)*B'*C*Z;