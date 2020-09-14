function [Z] = calc_Z(B,C)

temp1 = B*(inv(B'*B))*B';
deri_Z = C'*(eye(size(temp1)) - temp1)*C;
[V_Z,D_Z] = eig(deri_Z,'vector');
[d,ind] = sort(abs(D_Z));    % sorting eigenvalues from minimum to maximum 
Z = V_Z(:,ind(1));      % eigenvector corresponding to minimum eigenvalue
