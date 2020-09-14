function [B,C] = form_B_C(points_3D,points_2D)

N = size(points_3D,2);
for i = 1:1:N
    M(:,i) = points_3D(:,i);
    c(i) = points_2D(1,i);
    r(i) = points_2D(2,i);
    B(2*i-1,:) = [M(:,i)' 1 zeros(1,3) 0 -c(i)];
    B(2*i,:) = [zeros(1,3) 0 M(:,i)' 1 -r(i)];
    C(2*i-1,:) = -c(i)*M(:,i)';
    C(2*i,:) = -r(i)*M(:,i)';
end