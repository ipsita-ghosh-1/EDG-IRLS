function [X_c_handle,X_c_full] = get_matrix_handle(X_c,meas_op,weight_op,n,varargin)
% This function computes the handles to compute matrix-vector
% multiplications with the matrix X represented by X_c and with its
% (conjugate) transpose X' as well.
% =========================================================================
% Returns
% ----------
% X_c_handle, cell. 
%       X_c_handle{1} = handle such that X_c_handle{1}(v) = X*v;
%       X_c_handle{2} = handle such that X_c_handle{2}(v) = X'*v;
%       X_c_handle{3} = n.d1, number of rows of X.
%       X_c_handle{4} = n.d2, number of columns of X.
%  (if forcefull = true):
%       X_c_handle{5} = dense reprensentation of X_c as (n.d1 x n.d2)
%       matrix.
% =========================================================================
% Author: Christian Kuemmerle, kuemmerle@jhu.edu, 2021.
forcefull = false;
if ~isempty(varargin)
    for tt = 1:2:length(varargin)
        switch lower(varargin{tt})
            case 'forcefull'
                forcefull = varargin{tt+1};
        end
    end
end
U = weight_op.U;
if weight_op.symmetricflag
    V      = weight_op.V;
    X_c.Gam2 = X_c.Gam3';
else
    V      = weight_op.V;
end

if strcmp(meas_op.problem_type,'CommunityDetection')  || strcmp(meas_op.problem_type,'MaxCut CD') ...
        || strcmp(meas_op.problem_type,'MaxCut')
    linobjective = 1;
else
    linobjective = 0;
end
if linobjective
    nn = n.d1;
    if isfield(meas_op,'efficiencymode') && strcmp(meas_op.efficiencymode,'fast')
        yPhiAd = X_c.res_range(1:nn).';
        if strcmp(meas_op.problem_type,'CommunityDetection')
            yAvg = X_c.res_range(nn+1);
            if forcefull
                sps_plc = meas_op.sps_plc_fast;
                setSval(sps_plc,yPhiAd,nn);
                sps_plc = sps_plc + X_c.Afac.*meas_op.A+ yAvg.*(ones(nn,nn)./nn);
            end
        elseif strcmp(meas_op.problem_type,'MaxCut CD')
            if forcefull
                sps_plc = meas_op.sps_plc;
                setSval(sps_plc,yPhiAd,nn);
                sps_plc = sps_plc + X_c.Afac.*(2.*meas_op.A+speye(nn)-ones(nn,nn));
            end
        elseif strcmp(meas_op.problem_type,'MaxCut')
            if forcefull
                sps_plc = meas_op.sps_plc;
                setSval(sps_plc,yPhiAd,nn);
                sps_plc = sps_plc + X_c.Afac.*meas_op.A;
            end
        end
        
    else
        if strcmp(meas_op.problem_type,'CommunityDetection') || strcmp(meas_op.problem_type,'MaxCut')
            sps_plc = reshape(meas_op.Phi'*(X_c.res_range).',[n.d1,n.d2])+X_c.Afac.*meas_op.A;
        elseif strcmp(meas_op.problem_type,'MaxCut CD')
            sps_plc = reshape(meas_op.Phi'*(X_c.res_range).',[n.d1,n.d2])+...
                X_c.Afac.*(2.*meas_op.A+speye(nn)-ones(nn,nn));
        end
    end
else
    if strcmp(meas_op.problem_type,'MaxCut fixedval')
        nn = n.d1;
        sps_plc_fast = meas_op.sps_plc_fast;
        setSval(sps_plc_fast,X_c.res_range(1:nn),nn);
        sps_plc = sps_plc_fast + (X_c.res_range(nn+1)/meas_op.normA).*meas_op.A;
    elseif strcmp(meas_op.problem_type,'MultidimensionalScaling')
%         sps_plc_y = reshape(meas_op.Phi_W'*(X_c.ytilde).',[n.d1,n.d2]);
%         sps_plc = reshape(meas_op.Phi'*(X_c.XcOm).',[n.d1,n.d2]);
        if contains(meas_op.mode_linsolve,'tangspace')
%             sps_plc = meas_op.kappa(meas_op.kappainv(meas_op.kappainv(reshape(meas_op.Samp_opmat'*(X_c.res_range).',n.d1,n.d1))));
            sps_plc = reshape(applyMeasopBackward((X_c.res_range).',meas_op),[n.d1,n.d2]);%(-2).*meas_op.kappa(meas_op.kappainv(reshape(meas_op.Samp_opmat'*(X_c.res_range).',n.d1,n.d1)));
        else
            if isfield(meas_op,'efficiencymode') && strcmp(meas_op.efficiencymode,'fast')
                sps_plc_fast = MeasopBackwardMDS_sparsepart((X_c.res_range).',meas_op);
            else
                sps_plc = reshape(meas_op.Phi'*(X_c.res_range).',[n.d1,n.d2]);
            end
        end
    else
        sps_plc = meas_op.sps_plc;
        setSval(sps_plc,X_c.res_range,length(X_c.res_range));
    end
end
if isfield(meas_op,'efficiencymode') && (strcmp(meas_op.efficiencymode,'fast') || strcmp(meas_op.efficiencymode,'buildfull')) && not(forcefull)
    X_c_full = [];
% elseif strcmp(meas_op.problem_type,'MultidimensionalScaling')
%     X_c_full = (X_c.U*(X_c.Gam1*(V') +(X_c.Gam2))...
%         + X_c.Gam3*(V') +sps_plc); %+sps_plc_y + sps_plc);
elseif strcmp(meas_op.problem_type,'MultidimensionalScaling') && ~contains(meas_op.mode_linsolve,'tangspace') && isfield(meas_op,'efficiencymode') && strcmp(meas_op.efficiencymode,'fast')
    X_c_full = (U*(X_c.Gam1*(V') +(X_c.Gam2))...
         + X_c.Gam3*V' +sps_plc_fast + ...
         ones(n.d1,1)*X_c.res_range(meas_op.m+1:1:end)./sqrt(n.d1)+...
         X_c.res_range(meas_op.m+1:1:end).'*ones(1,n.d2)./sqrt(n.d1)); %+sps_plc_y + sps_plc);
else
    X_c_full = U*(X_c.Gam1*V' +X_c.Gam2)...
        + X_c.Gam3*V' +sps_plc;
%     X_c_full = real(reshape(meas_op.invWWT05*X_c_full(:),n.d1,n.d2));
end
% if strcmp(objective,'pluseps_squared')    
%     X_c_handle = X_c.U*(X_c.Gam1*V'+X_c.Gam2)+X_c.Gam3*V'+...
%         sps_plc;
% else
if linobjective && isfield(meas_op,'efficiencymode') && strcmp(meas_op.efficiencymode,'fast')
    if strcmp(meas_op.problem_type,'CommunityDetection')
        Xc_hdl  = @(x) (U*(X_c.Gam1*(V'*x) +(X_c.Gam2*x))...
            + X_c.Gam3*(V'*x) +yPhiAd.*x+ yAvg.*ones(nn,1).*sum(x)./nn+ X_c.Afac.*(meas_op.A*x));
        Xct_hdl = @(y) (V*(X_c.Gam1'*(U'*y)+X_c.Gam3'*y)...
            + X_c.Gam2'*(U'*y) +conj(yPhiAd).*y+ yAvg.*ones(nn,1).*sum(y)./nn+ X_c.Afac.*(meas_op.A'*y));
    elseif strcmp(meas_op.problem_type,'MaxCut CD')
        Xc_hdl  = @(x) (U*(X_c.Gam1*(V'*x) +(X_c.Gam2*x))...
            + X_c.Gam3*(V'*x) +(yPhiAd+ X_c.Afac).*x-(X_c.Afac.*ones(nn,1)).*sum(x)+ ...
            (2*X_c.Afac).*(meas_op.A*x));
        Xct_hdl = @(y) (V*(X_c.Gam1'*(U'*y)+X_c.Gam3'*y)...
            + X_c.Gam2'*(U'*y) +(conj(yPhiAd)+X_c.Afac).*y-(X_c.Afac.*ones(nn,1)).*sum(y)+ ...
            (2*X_c.Afac).*(meas_op.A'*y));
            %yAvg.*ones(nn,1).*sum(y)./nn+ X_c.Afac.*(meas_op.A'*y));
    elseif strcmp(meas_op.problem_type,'MaxCut')
        Xc_hdl  = @(x) (U*(X_c.Gam1*(V'*x) +(X_c.Gam2*x))...
            + X_c.Gam3*(V'*x) +yPhiAd.*x + X_c.Afac.*(meas_op.A*x));
        Xct_hdl = @(y) (V*(X_c.Gam1'*(U'*y)+X_c.Gam3'*y)...
            + X_c.Gam2'*(U'*y) +conj(yPhiAd).*y + X_c.Afac.*(meas_op.A'*y));
    end
elseif strcmp(meas_op.problem_type,'MultidimensionalScaling') && ~contains(meas_op.mode_linsolve,'tangspace') && isfield(meas_op,'efficiencymode') && strcmp(meas_op.efficiencymode,'fast')
    Xc_hdl  = @(x) U*(X_c.Gam1*(V'*x) +(X_c.Gam2*x))...
        + X_c.Gam3*(V'*x) +sps_plc_fast*x+ ...
        ones(n.d1,1)*(X_c.res_range(meas_op.m+1:end)*x)./sqrt(n.d1)+...
        X_c.res_range(meas_op.m+1:end).'*sum(x)./sqrt(n.d1);
    Xct_hdl = @(y) V*(X_c.Gam1'*(U'*y)+X_c.Gam3'*y)...
        + X_c.Gam2'*(U'*y) +sps_plc_fast'*y+ ...
        ones(n.d1,1)*(X_c.res_range(meas_op.m+1:end)*y)./sqrt(n.d1)+...
        X_c.res_range(meas_op.m+1:end).'*sum(y)./sqrt(n.d1);
elseif forcefull
    Xc_hdl  = @(x) X_c_full*x;
    Xct_hdl = @(y) X_c_full'*y;
% elseif strcmp(meas_op.problem_type,'MultidimensionalScaling')
%         Xc_hdl  = @(x) (U*(X_c.Gam1*(V'*x) +(X_c.Gam2*x))...
%         + X_c.Gam3*(V'*x) +sps_plc_y*x+sps_plc*x);
%     Xct_hdl = @(y) (V*(X_c.Gam1'*(U'*y)+X_c.Gam3'*y)...
%         + X_c.Gam2'*(U'*y) +sps_plc_y'*y+sps_plc'*y);
else
    Xc_hdl  = @(x) U*(X_c.Gam1*(V'*x) +(X_c.Gam2*x))...
        + X_c.Gam3*(V'*x) +sps_plc*x;
    Xct_hdl = @(y) V*(X_c.Gam1'*(U'*y)+X_c.Gam3'*y)...
        + X_c.Gam2'*(U'*y) +sps_plc'*y;
end
if not(isempty(X_c_full))
    X_c_handle={Xc_hdl,Xct_hdl,n.d1,n.d2,X_c_full};
else
    X_c_handle={Xc_hdl,Xct_hdl,n.d1,n.d2};
end
end