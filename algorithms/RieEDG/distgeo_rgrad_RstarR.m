function [X_l,P_approx,iterate_alignment,time] = distgeo_rgrad_RstarR(X_0,samples,d,max_iter,Dist_sampled,X_true,result_print,diff_threshing,norm_thresh,rel_thresh)
% This algorithm implements the (X-M,R^*_omega R_omega(X-M)) minimization for the
% operator R^*_omega R_omega.
% X_0 is the initialization. The standard initialization is just 
% X_0 =R^*_omega R_omega(M)
% samples is the array with each row corresponding to the (i,j)-th
% coefficients used in the masking procedure. This is needed to evaluate 
% d is the dimension of the datapoints
% max_iter is the desired number of iterations for the algorithm
% X_true, while not something we have access to in practice, is the true
% underlying Gram matrix. This is useful for comparison when wanting some
% illustration of convergence rates
% result_print is either 0 or 1, and when its 1 the relative difference
% between X_true and X_l in Frobenius norm is returned.
% first step, build the mask from the masked distance matrix

r = d; %rank of the gram matrix
n = size(X_0,1);
L = n*(n-1)/2;
m = size(samples,1)/2;
rate = m/L;
norm_diffs = zeros(max_iter,1);

%initialize the iteration using the hard thresholding operator
%change this so that initialization happens outside of this algorithm
%entirely
[X_l,~,D_l,U_l]=hard_thresh(X_0,r,1);

D_hat = Dist_sampled;
JXJ = -1/2*(D_hat - 1/n*ones(n,1)*sum(D_hat,2)' - 1/n*sum(D_hat,1)'*ones(1,n) + 1/n^2*sum(sum(D_hat))*ones(n));
P_JXJ = P_omega(JXJ,samples);
RstarR_Xtrue  = P_JXJ - diag(sum(P_JXJ,1));
time = zeros(1,max_iter);
tic;
for l=1:max_iter
    %compute the gradient
    G_l = (1/rate)^2*(RstarR_Xtrue - RstarR_omega(X_l,samples));

    %build the proper step size
    PU_Gl = U_l*(U_l'*G_l);
    Pt_Gl = PU_Gl + PU_Gl' -(PU_Gl*U_l)*U_l';
    RstarR_Pt_Gl = (1/rate)^2*RstarR_omega(Pt_Gl,samples);
    alpha_l=(norm(Pt_Gl,'fro')^2)/sum(sum(Pt_Gl.*RstarR_Pt_Gl));

    %this step is a way to more efficiently compute the next iteration
    %following the Wei et al 2016 paper
    [Q,R] = qr(alpha_l*(eye(n)-U_l*U_l')*G_l*U_l,'econ');
    M_l = [D_l + alpha_l*U_l'*G_l*U_l,    R' ;
          R            ,            zeros(r)];
    
    
    %[U,D_l] = eigs(M_l,d,'largestreal');
    [U,D_l] = eig(M_l);
    D_l = diag(D_l);
    D_l = real(D_l);
    [D_l,idx] = sort(D_l,'descend');
    D_l = D_l(1:r);
    U = U(:,idx(1:r));
    D_l = diag(D_l);
    %project back to the SPD cone if there are any negative eigenvalues
    for i=1:r
        if D_l(i,i)<0
            D_l(i,i)=0;
        end
    end

    %build new unitary matrix
    U_l = [U_l Q]*U(:,1:r);
    %construct new iterate
    X_l1 = X_l;
    X_l = U_l*D_l*U_l';
    iterate_alignment = U_l*sqrt(D_l);

    if result_print == 1 && not(isempty(X_true))
        norm_diffs(l) = norm(X_true-X_l,'fro')/norm(X_true,'fro');
    end
    
    time(l) = toc;
    if diff_threshing == 1 && not(isempty(X_true))
        if norm(X_true-X_l,'fro')/norm(X_true,'fro') < norm_thresh
            break
        end
    else
        if norm(X_l1-X_l,'fro')/norm(X_l,'fro') < rel_thresh
            break
        end
    end

end

P_approx = classical_edg(X_l,d,0);
time = time(1:l);
return
