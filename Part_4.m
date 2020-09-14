clear all
close all
clc

points_2D_left = read_points('project2_data\calibration\pts_2D_left.txt');
points_2D_right = read_points('project2_data\calibration\pts_2D_right.txt');
points_3D = read_points('project2_data\calibration\pts_3D.txt');
left_face = imread('project2_data\faceimage\left_face.jpg');
right_face = imread('project2_data\faceimage\right_face.jpg');
P_left = compute_P(points_3D,points_2D_left);
[Wl,Rl,Tl] = calc_W_R_T(P_left);

P_right = compute_P(points_3D,points_2D_right);
[Wr,Rr,Tr] = calc_W_R_T(P_right);
[R,T,F] = form_F(Wl,Wr,Rl,Rr,Tl,Tr);

rL1 = T'/norm(T);
rL2 = [rL1(2) -rL1(1) 0]/sqrt(rL1(2)^2 + rL1(1)^2); 
rL3 = cross(rL1,rL2);

RL = [rL1 ; rL2 ; rL3];

left_face_rect = zeros(size(left_face));
right_face_rect = zeros(size(right_face));

W = 0.5*(Wl+Wr);
RR = RL*R;
Wl(1,1) = W(2,2);
Wl(2,2) = W(2,2);
Wr(1,1) = W(2,2);
Wr(2,2) = W(2,2);
W(1,1)  = W(2,2);

Wl(1,2) = 0;
Wl(2,1) = 0;
Wr(1,2) = 0;
Wr(2,1) = 0;
W(1,2) = 0;
W(2,1) = 0;

 for cl_prime = 1:size(left_face,2)
    for rl_prime = 1:size(left_face,1)
        %rectifying left image
        scaled_2D_left_org = inv(W*RL*inv(Wl))*[cl_prime; rl_prime ; 1];
        scaled_2D_right_org = inv(W*RR*inv(Wr))*[cl_prime; rl_prime ; 1];
        lambda_l = scaled_2D_left_org(3);
        lambda_r = scaled_2D_right_org(3);
        cl = scaled_2D_left_org(1)/lambda_l;
        rl = scaled_2D_left_org(2)/lambda_l;
        cr = scaled_2D_right_org(1)/lambda_r;
        rr = scaled_2D_right_org(2)/lambda_r;
        if (cl >= 1 && cl < size(left_face,2) && rl>=1 && rl <= size(left_face,1))
            intensity_app = left_face(floor(rl),floor(cl),:);
            left_face_rect(rl_prime,cl_prime,:) = intensity_app;
        end
        if (cr >= 1 && cr < size(right_face,2) && rr>=1 && rr <= size(right_face,1))
            intensity_app = right_face(floor(rr),floor(cr),:);
            right_face_rect(rl_prime,cl_prime,:) = intensity_app;
        end
            
    end
end


[Vr,Dr] = eig(F,'vector');
indr = find(Dr == min(abs(Dr)));
[Vl,Dl] = eig(F','vector');
indl = find(Dl == min(abs(Dl)));
%Epipoles
er = Vr(:,indr);
el = Vl(:,indl);

point_1 = [805, 648];
point_2 = [871, 637];
point_3 = [758, 793];
point_4 = [817, 791];
cl = [point_1(1) point_2(1) point_3(1) point_4(1)];
rl = [point_1(2) point_2(2) point_3(2) point_4(2)];
cl_prime = zeros(size(cl));
rl_prime = zeros(size(rl));
for i = 1:4
    scaled_2D_left_org = W*RL*inv(Wl)*[cl(i); rl(i) ; 1];
    lambda_l = scaled_2D_left_org(3);
    cl_prime(i) = scaled_2D_left_org(1)/lambda_l;
    rl_prime(i) = scaled_2D_left_org(2)/lambda_l;
    [a(i), b(i), c(i)] = epiline(cl(i), rl(i),F);
end

cr = 1:20:size(right_face,2);
rr1 = (1/b(1))*(-a(1)*cr -c(1)); 
rr2 = (1/b(2))*(-a(2)*cr -c(2)); 
rr3 = (1/b(3))*(-a(3)*cr -c(3)); 
rr4 = (1/b(4))*(-a(4)*cr -c(4)); 

cr_prime = 1:20:size(right_face_rect,2);
rr1_prime = rl_prime(1)*ones(size(cr_prime)); 
rr2_prime = rl_prime(2)*ones(size(cr_prime));
rr3_prime = rl_prime(3)*ones(size(cr_prime));
rr4_prime = rl_prime(4)*ones(size(cr_prime));
figure
subplot(1,2,1)
imshow(uint8(left_face));
title('left face before rectification')
hold on
scatter(cl,rl,36);
subplot(1,2,2)
imshow(uint8(right_face));
hold on 
plot(cr,rr1,'x',cr,rr2,'o',cr,rr3,'+',cr,rr4,'r')
title('right face before rectification')
figure
subplot(1,2,1)
imshow(uint8(left_face_rect));
title('left face after rectification')
hold on
scatter(cl_prime,rl_prime,36);
subplot(1,2,2)
imshow(uint8(right_face_rect));
hold on 
plot(cr_prime,rr1_prime,'x',cr_prime,rr2_prime,'o',cr_prime,rr3_prime,'+',cr_prime,rr4_prime,'r')
title('right face after rectification')


