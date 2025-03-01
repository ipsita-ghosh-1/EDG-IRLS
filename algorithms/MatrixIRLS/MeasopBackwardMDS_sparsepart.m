function sps_plc = MeasopBackwardMDS_sparsepart(y,meas_op)
if strcmp(meas_op.problem_type,'MultidimensionalScaling')
    m = meas_op.m;
    n = length(y)-m;

%     row_diag = zeros(n,1);
%     col_diag = zeros(n,1);
%     for l = 1:n
%         row_diag(l) = sum(y(meas_op.rowind_diag_set{l}));
%         col_diag(l) = sum(y(meas_op.colind_diag_set{l}));
%     end
    yy = y(1:m);
    Sps1 = sparse(meas_op.rowind,meas_op.rowind,yy ,n,n);%sparse(1:n,1:n,row_diag,n,n);
    Sps2 = sparse(meas_op.colind,meas_op.colind,yy ,n,n);%sparse(1:n,1:n,col_diag,n,n);
    Sps3 = meas_op.sps_upper;
    setSval(Sps3,y(1:m),m);
    Sps4 = Sps3';
    sps_plc = Sps1 + Sps2 - Sps3 - Sps4;
else
    error('Check measurement operator.')
end
end

