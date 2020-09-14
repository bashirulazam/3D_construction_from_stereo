clear all
close all
clc

points_2D_left = read_points('project2_data\calibration\pts_2D_left.txt');
points_2D_right = read_points('project2_data\calibration\pts_2D_right.txt');
points_3D = read_points('project2_data\calibration\pts_3D.txt');

N = size(points_2D_left,2);

% avg_l = sum(points_2D_left,2)/N;
% avg_r = sum(points_2D_right,2)/N;
% var_l = var(points_2D_left,1,2);
% var_r = var(points_2D_right,1,2);
% 
% xl = avg_l(1); 
% yl = avg_l(2);
% dl = sqrt(var_l(1)^2 + var_l(2)^2);
% 
% xr = avg_r(1); 
% yr = avg_r(2);
% dr = sqrt(var_r(1)^2 + var_r(2)^2);
% 
% Tl = [1/dl 0 -xl/dl;
%     0 1/dl -yl/dl;
%     0 0 1];
% Tr = [1/dr 0 -xr/dr;
%     0 1/dr -yr/dr;
%     0 0 1];


for i = 1:N
    Ul(:,i) = [points_2D_left(:,i) ; 1];
    Ur(:,i) = [points_2D_right(:,i) ; 1];
    cr(i) = Ur(1,i);
    rr(i) = Ur(2,i);
    A(i,:) = [Ul(:,i)'*cr(i) Ul(:,i)'*rr(i) Ul(:,i)'];
end

[V,D] = eig(A'*A,'vector');
D_min_index = find(D==min(D));
F_solve = V(:,D_min_index);
F_est = [F_solve(1:3) F_solve(4:6) F_solve(7:9)];
[U,S,V] = svd(F_est);
S(3,3) = 0 ;
F_ref = U*S*V';