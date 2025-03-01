function val = norm_frob_compact(M,meas_op,weight_op,mode)
% Calculates the Frobenius norm ||M||_F of a matrix M given in the the
% compact format
% M_mat =  U^{(k)}*(\tilde{Gamma}_1*V^{(k)'}+\tilde{Gamma}_2)
%          + \tilde{Gamma}_3*V^{(k)'} + P_{\Omega}^*(r_{k+1}).
% =========================================================================
% Parameters
% ----------
% M:    struct. Represents a matrix M_mat, is given such that,
%       M.Gam1 = \tilde{Gamma}_1  [(r_M x r_M) matrix]
%       M.Gam2 = \tilde{Gamma}_2  [(r_M x d_2) matrix]
%       M.Gam3 = \tilde{Gamma}_3  [(d_1 x r_M) matrix]
%       M.res_range = r_{k+1}     [(m x 1) vector].
%       M.U    = U^{(k)}          [(d_1 x r_M) matrix w/
%                                   orthonormal columns]
%       V    = V^{(k)}          [(d_2 x r_M) matrix w/
%                                   orthonormal columns
% meas_op: struct. Contains data of measurement operator. Usually (for
%   completion problems, contains field:
%  sps_plc: (d_1 x d_2) sparse matrix with support of size m. Non-zeros 
%       correspond to indices \Omega (values are irrevelant).
% mode: string of characters. If
%        = 'standard': As above.
%        = 'lsqr': (to be updated): M = {M1,M2,M3,M4} such that
%           Mat = [M1         ,  M2*V_perp;
%                  U_perp'*M3 ,  U_perp'*\mathcal{A}'(M4)*V_perp],
%          and U_perp resp. V_perp are orthonormal bases of the
%          spaces that are orthogonal on the columns of U and
%          V, respectively.
% Returns
% ----------
% err_fro: double. Frobenius norm ||M||_F.

%              mode: Indicates if the format of M is such that:
%                    = 'standard': M = {M1,M2,M3,M4} such that
%                                  Mat = U*M1*V' +
%                                  U*M2'+M3*V'+\mathcal{A}'(M4),
%                      and \mathcal{A} is the operator that maps to the
%                      entries represented by sps_plc.
%                    = 'lsqr': M = {M1,M2,M3,M4} such that
%                       Mat = [M1         ,  M2*V_perp;
%                              U_perp'*M3 ,  U_perp'*\mathcal{A}'(M4)*V_perp],
%                      and U_perp resp. V_perp are orthonormal bases of the
%                      spaces that are orthogonal on the columns of U and
%                      V, respectively.
% =========================================================================
% Author: Christian Kuemmerle, 2019-2021.
U = weight_op.U;
d1 = size(U,1);
if not(weight_op.symmetricflag)
    V      = weight_op.V;
    d2      =size(V,1);
else
    V      = U;
    M.Gam2 = M.Gam3';
end
d2      =size(V,1);
if nargin == 3
    mode = 'standard';
end
if strcmp(mode,'standard')
    if strcmp(meas_op.problem_type,'PhaseRetrieval')
        val = sc_prod_IRLS_PR_Astr(M,M,weight_op.AU,meas_op.A_corr);
        val = sqrt(val);
    else
        if strcmp(meas_op.problem_type,'MultidimensionalScaling')
            m=meas_op.m;%length(M.res_range);
        else
            m=length(M.res_range);
        end
        if strcmp(meas_op.problem_type,'MatrixCompletion')
            sps_plc = meas_op.sps_plc;
        else
            if strcmp(meas_op.problem_type,'MultidimensionalScaling')
                if isfield(meas_op,'efficiencymode') && strcmp(meas_op.efficiencymode,'fast')
                n = d1;
                m = length(M.res_range)-n;
                row_diag = zeros(n,1);
                col_diag = zeros(n,1);
                for l = 1:n
    %                     rowind_set = find(meas_op.rowind==l)';
    %                     colind_set = find(meas_op.colind==l)';
                    row_diag(l) = sum(M.res_range(meas_op.rowind_diag_set{l}));
                    col_diag(l) = sum(M.res_range(meas_op.colind_diag_set{l}));
                end

                Sps1 = sparse(1:n,1:n,row_diag,n,n);%meas_op.sps_diag_rows;
                Sps2 = sparse(1:n,1:n,col_diag,n,n);%meas_op.sps_diag_cols;
                Sps3 = meas_op.sps_upper;
    %                 Sps4 = meas_op.sps_lower;
                %setSval(meas_op.sps_diag_rows,M.res_range(1:m),m);%meas_op.sps_diag_rows;
    %                 setSval(meas_op.sps_diag_cols,M.res_range(1:m),m);%meas_op.sps_diag_cols;
    %                 setSval(meas_op.sps_upper,M.res_range(1:m),m);%meas_op.sps_upper;
                setSval(Sps3,M.res_range(1:m),m);
                Sps4 = Sps3';
    %                 setSval(meas_op.sps_lower,M.res_range(1:m),m);%meas_op.sps_lower;
                sps_plc = Sps1 + Sps2 - Sps3 - Sps4; 
                else
                    if contains(meas_op.mode_linsolve,'tangspace')
                        sps_plc = reshape(applyMeasopBackward((M.res_range).',meas_op),[d1,d1]);
                        %meas_op.kappa(meas_op.kappainv(meas_op.kappainv(reshape(meas_op.Samp_opmat'*(M.res_range).',size(U,1),size(U,1)))));
                    else
                        sps_plc = reshape(meas_op.Phi'*(M.res_range).',[d1,d2]);
                    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    % doesn't work yet
    %                 sps_plc_y = reshape(meas_op.Phi_W'*(M.ytilde).',[d1,d2]);
    %                 sps_plc = reshape(meas_op.Phi'*(M.XcOm).',[d1,d2]);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                end
    %             Sps5 = (M.res_range(m+1:end)'*ones(1,n)+ones(n,1)*M.res_range(m+1:end))./sqrt(2*n); % used for testing
    %             sps_plc_test = Sps1 + Sps2 - Sps3 - Sps4 + Sps5; % used for testing
            else
                if strcmp(meas_op.efficiencymode,'fast')
                    yPhiAd = M.res_range(1:d1).';
                    if strcmp(meas_op.problem_type,'CommunityDetection')
                        yAvg = M.res_range(d1+1);
        %             elseif strcmp(meas_op.problem_type,'MaxCut CD')
        %                 
        %             elseif strcmp(meas_op.problem_type,'MaxCut')
        %                 
                    end
                else
                    if strcmp(meas_op.problem_type,'CommunityDetection') || strcmp(meas_op.problem_type,'MaxCut')
                        sps_plc = reshape(meas_op.Phi'*(M.res_range.'),[d1,d2])+M.Afac.*meas_op.A;
                    elseif strcmp(meas_op.problem_type,'MaxCut CD')
                        sps_plc = reshape(meas_op.Phi'*(M.res_range).',[d1,d2])+...
                            M.Afac.*(2.*meas_op.A+speye(d1)-ones(d1,d1));
                    elseif strcmp(meas_op.problem_type,'MaxCut fixedval')
                        nn = size(U,1);
                        sps_plc_fast = meas_op.sps_plc_fast;
                        setSval(sps_plc_fast,M.res_range(1:nn),nn);
                        sps_plc = sps_plc_fast + (M.res_range(nn+1)/meas_op.normA).*meas_op.A;
                    end
                end
            end
        end
        if (isfield(meas_op,'efficiencymode') && strcmp(meas_op.efficiencymode,'fast') && not(strcmp(meas_op.problem_type,'MatrixCompletion')))...
                || strcmp(meas_op.problem_type,'MultidimensionalScaling')
            if strcmp(meas_op.problem_type,'CommunityDetection') || strcmp(meas_op.problem_type,'MaxCut')
                M4U = conj(yPhiAd).*U;
                M4U = M4U + M.Afac.*meas_op.A*U;

                M4V = yPhiAd.*V;
                M4V = M4V + M.Afac.*meas_op.A*V;
                if strcmp(meas_op.problem_type,'CommunityDetection')
                    M4U = M4U + conj(yAvg).*ones(d2,1)*sum(U,1)./d2;
                    M4V = M4V + yAvg.*ones(d1,1)*sum(V,1)./d1;
                end
            elseif strcmp(meas_op.problem_type,'MaxCut CD')
                M4U = conj(yPhiAd+M.Afac).*U;
                M4U = M4U + conj(M.Afac.*(2.*meas_op.A))*U;
                M4U = M4U - conj(M.Afac).*ones(d1,1).*sum(U,1);
                M4V = (yPhiAd+M.Afac).*V;
                M4V = M4V + M.Afac.*(2.*meas_op.A)*V;
                M4V = M4V - M.Afac.*ones(d1,1).*sum(V,1);
            elseif strcmp(meas_op.problem_type,'MultidimensionalScaling')
    %             setSval(sps_plc,M.res_range,m+n);
                M4U = sps_plc'*U; % d2 x r
                M4V = sps_plc*V;  % d1 x r 
                if isfield(meas_op,'efficiencymode') && strcmp(meas_op.efficiencymode,'fast')
                    M4U = M4U + (ones(n,1)*(M.res_range(m+1:end)*U) +M.res_range(m+1:end).'*sum(U,1))./sqrt(2*n);
                    M4V = M4V + (ones(n,1)*(M.res_range(m+1:end)*V) +M.res_range(m+1:end).'*sum(V,1))./sqrt(2*n);
                else
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    % doesn't work yet
    %                 M4U = M4U + sps_plc_y'*U;
    %                 M4V = M4V + sps_plc_y*V;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                end
            end
        else
            setSval(sps_plc,M.res_range,m);
            M4U=sps_plc'*U; % d2 x r
            M4V=sps_plc*V;  % d1 x r 
        end
        temp=M.Gam1(:);
        val=sum(temp'*temp);
        val=val+2*real(trace(M.Gam1*(V'*M.Gam2')));
        val=val+2*real(trace((M.Gam3'*U)*M.Gam1));
        temp=M.Gam2(:);
        val=val+sum(temp'*temp);
        val=val+2*real(trace(M.Gam2*M4U));
        val=val+2*real(trace(M4V'*M.Gam3));
        temp=M.Gam3(:);
        val=val+sum(temp'*temp);
        val=val+2*real(trace(M.Gam1*(M4V'*U)));
        val=val+2*real(trace((M.Gam3'*U)*(M.Gam2*V)));

    %     if strcmp(meas_op.problem_type,'MultidimensionalScaling') && not(isfield(meas_op,'efficiencymode') && strcmp(meas_op.efficiencymode,'fast'))
    %         temp=M.ytilde(:); % check if this is correct
    %         val=val+sum(temp'*temp); 
    %         temp = M.XcOm(:);
    %     else
            temp=M.res_range(:);
    %     end
        val=val+sum(temp'*temp);
        val=sqrt(val);
    end
elseif strcmp(mode,'lsqr')
    m = length(M.res_range);
    setSval(sps_plc,M.res_range,m);
    UAstM4 = U'*sps_plc;
    UAstM4V = UAstM4*V;
    AstM4V = sps_plc*V; 
    
    temp=M.Gam2(:);
    val2=sum(temp'*temp);
    M2V =M.Gam2*V;
    temp=M2V(:);
    val2=val2-sum(temp'*temp);
    
    temp=M.Gam3(:);
    val3=sum(temp'*temp);
    UM3 =U'*M.Gam3;
    temp=UM3(:);
    val3=val3-sum(temp'*temp);
    
    temp=M.res_range(:);
    val4=sum(temp'*temp);
    temp=UAstM4(:);
    val4=val4-sum(temp'*temp);
    temp=UAstM4V(:);
    val4=val4+sum(temp'*temp);
    temp=AstM4V(:);
    val4=val4-sum(temp'*temp);
    
    temp=M.Gam1(:);
    val=sum(temp'*temp)+val2+val3+val4;
    val=sqrt(val);
end
end

function val = sc_prod_IRLS_PR_Astr(X,Y,weight_op,meas_op)
%sc_prod_IRLS_PR This function calculates the scalar product of
% two matrices (old)
% X = U_T*X_1*U_T' + U_T*X_2*conj(A) + A.'*X_2'*U_T'+ A.'*diag(X_3)*conj(A)
% and 
% Y = U_T*Y_1*U_T' + U_T*Y_2*conj(A) + A.'*Y_2'*U_T'+ A.'*diag(Y_3)*conj(A)
%
% two matrices 
% X = U_T*X.Gam1*U_T' + (U_T*X.Gam2' + X.Gam2*U_T')/sqrt(2) + A.'*diag(X.res_range.')*conj(A)
% and 
% Y =U_T*Y.Gam1*U_T' + (U_T*Y.Gam2' + Y.Gam2*U_T')/sqrt(2) + A.'*diag(Y.res_range.')*conj(A)
% whose coefficients are saved as (old)
% X = cell(1,3) such that X_1 is (R x R) matrix
%                         X_2 is (R x m) matrix
%                         X_3 is (m x 1) matrix         and
% Y = cell(1,3) such that Y_1 is (R x R) matrix
%                         Y_2 is (R x m) matrix
%                         Y_3 is (m x 1) matrix.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% X is struct with fields: X.Gam1 (R x R) matrix (Hermitian/symmetric)
%                          X.Gam2 is (n x R) matrix
%                          X.res_range is (1 x m) matrix         and
% Y is struct with fields: Y.Gam1 (R x R) matrix (Hermitian/symmetric)
%                          Y.Gam2 is (n x R) matrix
%                          Y.res_range is (1 x m) matrix         and
% To be used in HM_IRLS_Phaseretrieval.m.
% The code makes use of precalculated products AU=conj(A)*U_T etc. to save
% calculations.

%       Input:    X: (1 x 3) cell, struture see above
%                 Y: (1 x 3) cell, struture see above
%                AU: (m x R) matrix such that AU=conj(A)*U_T
%            A_corr: (m x m) matrix such that A_corr= conj(A)*A.',
%                                   hopefully sparse
%         A_corr_sq: (m x m) matrix such that 
%                                   A_corr_sq = A_corr.*conj(A_corr)
%      Output:   
U = weight_op.U;
AU = weight_op.AU;


% [m,R]=size(AU);
tmp1= conj(AU)*Y.Gam1.';
% tmp2= A_corr*X.Gam2.';

val=sum(sum(conj(X.Gam1).*Y.Gam1));
% val=val+2.*real(sum(sum((AU*X.Gam1').*Y.Gam2.')));%2.*real(sum(sum(Y.Gam2.'.*(AU))))
val = val + sqrt(2).*real(trace(Y.Gam2'*U*X.Gam1'));
% val=val+sum(sum((X.res_range.'.*A_corr).*(Y.res_range.'.*A_corr).'));
% val=val+sum(sum((conj(AU)*conj(X.Gam1)).*(Y.res_range.'.*AU)));
val = val + sum(sum(conj(AU*X.Gam1).*(Y.res_range.'.*AU)));
val = val + sqrt(2).*real(trace(X.Gam2'*U*Y.Gam1));
val = val + real(trace((X.Gam2'*U)*(Y.Gam2'*U))+trace(X.Gam2'*Y.Gam2));
val = val + sqrt(2).*real(sum(sum(meas_op.A(conj(X.Gam2)).*(Y.res_range.'*(AU)))));
val = val + sum(sum(conj(AU*X.Gam1').*(X.res_range.'.*AU)));
val = val + sqrt(2).*real(sum(sum(meas_op.A(conj(Y.Gam2)).*(X.res_range.'*(AU)))));
val = val  +sum(sum((X.res_range.'.*A_corr).*(Y.res_range.'.*A_corr).'));

% val=val+sum(sum(tmp1.*(X.res_range.'.*AU)));
% val=val+2.*real(sum(sum(tmp1.*X.Gam2')));
% val=val+2.*real(sum(sum(tmp2.*Y.Gam2.')));
% val=val+2.*real(sum(sum((conj(Y.Gam2)*conj(AU)).*(AU.'*X.Gam2.'))));
% val=val+2.*real(sum(sum((Y.res_range.'.*conj(AU)).*tmp2)));
% val=val+2.*real(sum(sum((X.res_range.'.*AU).*(conj(A_corr)*Y.Gam2'))));
end