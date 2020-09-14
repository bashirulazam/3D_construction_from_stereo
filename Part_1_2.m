clear all
close all
clc

points_2D_left = read_points('project2_data\calibration\pts_2D_left.txt');
points_2D_right = read_points('project2_data\calibration\pts_2D_right.txt');
points_3D = read_points('project2_data\calibration\pts_3D.txt');

P_left = compute_P(points_3D,points_2D_left);
[Wl,Rl,Tl] = calc_W_R_T(P_left);

P_right = compute_P(points_3D,points_2D_right);
[Wr,Rr,Tr] = calc_W_R_T(P_right);
[R,T,F] = form_F(Wl,Wr,Rl,Rr,Tl,Tr);



