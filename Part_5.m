clear all
close all
clc

points_2D_left_cal = read_points('project2_data\calibration\pts_2D_left.txt');
points_2D_right_cal = read_points('project2_data\calibration\pts_2D_right.txt');
points_3D_cal = read_points('project2_data\calibration\pts_3D.txt');
P_left = compute_P(points_3D_cal,points_2D_left_cal);
[Wl,Rl,Tl] = calc_W_R_T(P_left);
P_right = compute_P(points_3D_cal,points_2D_right_cal);
[Wr,Rr,Tr] = calc_W_R_T(P_right);
R =Rl*Rr';
T = Tl - R*Tr;

points_2D_left_face = read_points('project2_data\faceimage\pts_left.txt');
points_2D_right_face = read_points('project2_data\faceimage\pts_right.txt');

P = [Wr*R Wr*T];

for i = 1:size(points_2D_left_face,2)
    cl = points_2D_left_face(1,i);
    rl = points_2D_left_face(2,i);
    cr = points_2D_right_face(1,i);
    rr = points_2D_right_face(2,i);

    A = [Wl(1,1) Wl(1,2) Wl(1,3) -cl 0;
        Wl(2,1) Wl(2,2)  Wl(2,3) -rl 0;
        Wl(3,1) Wl(3,2)  Wl(3,3) -1 0;
        P(1,1) P(1,2)   P(1,3)  0   -cr;
        P(2,1)  P(2,2) P(2,3)   0   -rr;
        P(3,1)  P(3,2) P(3,3)   0   -1];
    b = [0; 0; 0; -P(1,4); -P(2,4); -P(3,4)];
    X = inv(A'*A)*A'*b;
    points_3D_face(1,i) = X(1);
    points_3D_face(2,i) = X(2);
    points_3D_face(3,i) = X(3);
end

plot3(points_3D_face(1,:),points_3D_face(2,:),points_3D_face(3,:),'x') 
fileID = fopen('3D_reconstructed_pts.txt','w');
fprintf(fileID,'%6f %6f %6f\n',points_3D_face');
fclose(fileID);


