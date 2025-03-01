function [W,W_long,meas_op,ind_diag,V] = create_Wmatrix_MDS(Phi,n)
%CREATE_WMATRIX_MDS Summary of this function goes here
%   Detailed explanation goes here
meas_op = struct;
meas_op.sps_meas = Phi;
meas_op.Omega = find(meas_op.sps_meas);
% meas_op.Omega = find(meas_op.sps_meas');
meas_op.m = length(meas_op.Omega);
[meas_op.rowind,meas_op.colind] = find(meas_op.sps_meas);
% [rowind_tmp,colind_tmp] = find(meas_op.sps_meas');
% meas_op.colind = rowind_tmp;
% meas_op.rowind = colind_tmp;

meas_op.OmegaT = sub2ind([n.d2,n.d1],meas_op.colind,meas_op.rowind) ;

% meas_op.Omega = sub2ind([n.d2,n.d1],meas_op.rowind,meas_op.colind) ;

meas_op.rowind_diag_set = cell(1,n.d1);
meas_op.colind_diag_set = cell(1,n.d1);
for l = 1:n.d1
    meas_op.rowind_diag_set{l} = find(meas_op.rowind==l)';
    meas_op.colind_diag_set{l} = find(meas_op.colind==l)';
    meas_op.unionind_diag_set{l} = [meas_op.rowind_diag_set{l},meas_op.colind_diag_set{l}];
end
maxx_ind = 1;
for l = 1:n.d1
    maxx_ind = max(length(meas_op.unionind_diag_set{l}),maxx_ind);
end
meas_op.IndMat = (meas_op.m+1).*ones(maxx_ind,n.d1);
for l = 1:n.d1
    meas_op.IndMat(1:length(meas_op.unionind_diag_set{l}),l) = meas_op.unionind_diag_set{l};
end 

meas_op.sps_diag_rows = sparse(meas_op.rowind,meas_op.rowind,ones(meas_op.m,1),n.d1,n.d2);
meas_op.sps_diag_cols = sparse(meas_op.colind,meas_op.colind,ones(meas_op.m,1),n.d1,n.d2);
meas_op.sps_upper = meas_op.sps_meas;
meas_op.sps_lower = meas_op.sps_meas.';



%     diagind = (1:dim + 1:dim*dim);
%     edgeind = I + (J - 1)*dim;
meas_op.linind_test = meas_op.rowind + (meas_op.colind - 1)*n.d1;
meas_op.Omega_transp = meas_op.colind + (meas_op.rowind - 1)*n.d1;
%%% Set up matrix Phi of size (m+n) x (n^2)
meas_op.rowind_diag_lin = meas_op.rowind + (meas_op.rowind - 1)*n.d1;
meas_op.colind_diag_lin = meas_op.colind + (meas_op.colind - 1)*n.d1;

colind_Phi = [meas_op.rowind_diag_lin,meas_op.colind_diag_lin,...
        meas_op.linind_test,meas_op.Omega_transp]';
%     meas_op.Omega,meas_op.Omega_transp]';
colind_Phi = colind_Phi(:);
rowind_Phi = ones(4,meas_op.m).*[1:meas_op.m];
rowind_Phi = rowind_Phi(:);
vals = [ones(2,meas_op.m);-ones(2,meas_op.m)];
vals = vals(:);

vals_short = vals;
colind_Phi_short = colind_Phi;
rowind_Phi_short = rowind_Phi;

colind_Phi_diag = [];
rowind_Phi_diag = [];
vals_diag = [];
for l=1:n.d1
    e_c = sparse(l,1,1,n.d1,1);
%         e_c(l) = 1;
    oness = ones(1,n.d1);
    colind_p2_c = e_c*spones(oness);
    colind_p2_c_1 = reshape(colind_p2_c,n.d1^2,1);
    colind_p2_c_2 = reshape(colind_p2_c',n.d1^2,1);
    colind_p2_c_mat = colind_p2_c_1 + colind_p2_c_2;
    vals_c = nonzeros(colind_p2_c_mat)./2;%(sqrt(2*n.d1));
    colind_p2_c_ind = find(colind_p2_c_mat);
    rolind_p2_c = (meas_op.m+l).*ones(length(vals_c),1);

    colind_Phi_diag = [colind_Phi_diag;colind_p2_c_ind];
    rowind_Phi_diag = [rowind_Phi_diag;rolind_p2_c];
    vals_diag = [vals_diag;vals_c];
end
ind_diag = sub2ind([meas_op.m+n.d1,n.d1^2],rowind_Phi_diag,colind_Phi_diag);

colind_Phi = [colind_Phi_short;colind_Phi_diag];
rowind_Phi = [rowind_Phi_short;rowind_Phi_diag];
vals       = [vals_short;vals_diag];


%     rowind_Phi = [1:n,n+1.*ones(1,s)];
%     colind_Phi = [inds,inds2];
%     vals = [ones(1,n),valsA./meas_op_out.normA];
WT = sparse(rowind_Phi_short,colind_Phi_short,vals_short,meas_op.m,n.d1^2);
W = WT';
WT_long = sparse(rowind_Phi,colind_Phi,vals,meas_op.m+n.d1,n.d1^2);
W_long = WT_long';
meas_op.Phi = W_long';
% Construct V
V = zeros(n.d1*n.d2,n.d1*n.d2);
onev = ones(n.d1,1);
unitvs = cell(1,n.d1);
for i=1:n.d1
    unitvs{i} = zeros(n.d2,1);
    unitvs{i}(i) = 1;
end
for i=1:n.d1
    a = (unitvs{i}-onev./n.d1);
    for j=1:n.d2
        if i == j
            vtmp = unitvs{i}*unitvs{i}'- a*a';
        else
            b = (unitvs{j}-onev./n.d1);
            vtmp = -0.5*(a*b'+b*a');
        end
        V(:,(i-1)*n.d2+j) = vtmp(:);
    end
end

end
