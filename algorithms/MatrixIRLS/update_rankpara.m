function [r_c,weight_op,eps]=update_rankpara(XXX,weight_op,eps,R,r_c,...
    N_SVD,tracking,opts,varargin)
%update_rankpara This functions determines how many singular values
%are larger than eps, and updates the rank r_c accordingly. It further
%outputs the first r_c singular values (sing) and first r_c singular vectors
%(U_X and V_X).
%
%       Input:  XXX = (1 x 4) cell with function handles of X and X', 
%                     respectively, and the dimensions d1 and d2.
%                 R = maximal number of singular values to be computed.
%               r_c = Number of singular values of last iterate larger than
%                     the eps of last iterate.
%             N_SVD = Max. number of iterations for bksvd
%                     
tol_precision = 1e-5; % precision parameter

decomposition = 'svd';
if ~isempty(varargin)
    for tt = 1:2:length(varargin)
        switch lower(varargin{tt})
            case 'decomposition'
                decomposition = varargin{tt+1};
        end
    end
end
U = weight_op.U;
V = weight_op.V;
sing = weight_op.sing;

d1=size(U,1);
d2=size(V,1);
d=min(d1,d2);
r_SVD=max(r_c,1);
r_correct = 0;
while not(r_correct)
%     if length(find(sing>eps+tol_precision)) <= r_SVD
% %         r_c = length(find(sing>eps+tol_precision));
% %         r_correct = 1;
%     else
        r_SVD = min(2*r_SVD,R);
        if ~tracking
            [U,singval_mat_c,V] = update_matrixdec_from_handle(XXX,...
            min(d,r_SVD+1),N_SVD,decomposition,opts);
%             if strcmp(decomposition,'svd')
%                 [U_X, singval_mat_c, V_X]=bksvd_mod(XXX(1:4), min(d,r_SVD+1), N_SVD);
%             elseif strcmp(decomposition,'svdexact')
%                 [U_X, singval_mat_c, V_X]=svds(XXX{5},min(d,r_SVD+1));
%             elseif strcmp(decomposition,'eig_abs')
%                 OPTS.isreal = true;
%                 TARGET = 'LA';
%                 [U_X, singval_mat_c]= primme_eigs(XXX{1},d,min(d,r_SVD+1),TARGET,OPTS);
%                 V_X = U_X;
%             elseif strcmp(decomposition,'eigexact_abs')
%                 [U_X, singval_mat_c] = eigs(XXX{5},min(d,r_SVD+1),'largestabs');
%                 singval_mat_c = abs(singval_mat_c);
%                 V_X = U_X;
%             elseif strcmp(decomposition,'eigexact')
%                 [U_X, singval_mat_c] = eigs(XXX{5},min(d,r_SVD+1),'largestreal');
%                 V_X = U_X;
%             else 
%                 error('Please indiciate whether svd or eig is used in the weight updates.')
%             end
            sing = diag(singval_mat_c);
        end
        if length(find(sing>eps+tol_precision)) <= r_SVD+1
            r_c = length(find(sing>eps+tol_precision));
            r_correct = 1;
        end
%             else
%                 error('Implement loop here to calculate even more singular values.')
%             end
%     end
    if r_SVD == R
        r_c=min(r_SVD,r_c);%min(r_c,r_SVD);
        r_correct = 1;
    end
end
if r_c == 0
    eps = 0.99*eps;
    r_c = length(find(sing>eps+tol_precision));
%     weight_op.sing = [];
%     weight_op.U = zeros(d1,r_c);
%     weight_op.V = zeros(d2,r_c);
    weight_op.sing = sing(1:r_c);
    weight_op.U = U(:,1:r_c);
    weight_op.V = V(:,1:r_c);
else
    weight_op.sing = sing(1:r_c);
    weight_op.U = U(:,1:r_c);
    weight_op.V = V(:,1:r_c);
end
if isfield(opts,'tracking') && opts.tracking
    weight_op.Uperp = null(weight_op.U');
    weight_op.Vperp = null(weight_op.V');
else
    weight_op.Uperp = [];
    weight_op.Vperp = [];
end
if isfield(weight_op,'symmetricflag') && weight_op.symmetricflag
    weight_op.symTmat = sparse(ones(r_c,r_c));
    weight_op.symTmat = triu(weight_op.symTmat,0);
    weight_op.symInd = find(weight_op.symTmat);
else
    weight_op.symTmat = [];
    weight_op.symInd = [];
end
end

