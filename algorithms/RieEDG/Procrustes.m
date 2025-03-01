function [P2_align,T,P2_trans] = Procrustes(P1,P2)
% the goal of this is to align a point matrix P2 with a given matrix P1

n = size(P1,1);
d = size(P1,2);

J = eye(n)-1/n*ones(n);
Psi = P1'*J*P2;
[U,~,V] = svds(Psi,3);
T = V*U';
P2_rot = P2*T;
P2_trans = (1/n)*(P1-P2_rot)'*ones(n,1);
P2_align = P2_rot + ones(n,1)*P2_trans';

return