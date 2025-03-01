function [X_l,P_approx,norm_diffs] = distgeo_rgrad_Pomega(X_0,samples,d,max_iter,X_true,rate,result_print,threshing,norm_thresh)
% P_omega_X is masked distance matrix
% d is dimension of points that generate P and gram matrix
% P is returned matrix of points, n x d

 
% first step, build the mask from the masked distance matrix
n = size(X_0,1);
I = find(X_0 ~= 0);
norm_diffs = zeros(max_iter,1); 
Pomega_true = P_omega(X_true,samples);

r = d+2; %rank of the EDM
[X_l,U_l,D_l,V_l]=hard_thresh(X_0,r,2);

for l=1:max_iter
    % form residual
    G_l = zeros(n);
    G_l(I) = (1/rate)*(Pomega_true(I) - X_l(I));

    %form optimal stepsize
    PU_Gl = U_l*(U_l'*G_l);
    PV_Gl = (G_l*V_l)*V_l';
    Pt_Gl = PU_Gl + PV_Gl -(PU_Gl*V_l)*V_l';
    Pomega_Pt_Gl = zeros(n);
    Pomega_Pt_Gl(I) =  (1/rate)*Pt_Gl(I);
    alpha_l=(norm(Pt_Gl,'fro')^2)/sum(sum(Pt_Gl.*Pomega_Pt_Gl));

    %following Wei et al 2016
    [Q1,R1] = qr(alpha_l*(eye(n)-V_l*V_l')*G_l*U_l,'econ');
    [Q2,R2] = qr(alpha_l*(eye(n)-U_l*U_l')*G_l*V_l,'econ');
    M_l = [D_l + alpha_l*U_l'*G_l*V_l,    R1' ;
           R2            ,            zeros(r)];

    %hard thresholding step
    [U,D_l,V] = svd(M_l);
    D_l = diag(D_l);
    D_l = D_l(1:r);
    D_l = diag(D_l);

    U_l = [U_l Q2]*U(:,1:r);
    V_l = [V_l Q1]*V(:,1:r);
    
    X_l1 = X_l;
    X_l = U_l*D_l*V_l';

    if result_print == 1
        norm_diffs(l) = norm(X_true-X_l,'fro')/norm(X_true,'fro');
    end
% 
    if threshing == 1
        if norm(X_true-X_l,'fro')/norm(X_true,'fro') < norm_thresh
            break
        end
    end
    
end

P_approx = classical_edg(X_l,d,0);

return






