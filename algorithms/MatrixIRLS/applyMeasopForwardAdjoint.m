function [PhiPhiT_y,PhiPhiT_handle] = applyMeasopForwardAdjoint(y,meas_op)
%applyMeasopForwardAdjoint This function computes the image
% (Phi *Phi')(y) of the vector 'y' with with respect to the concatenation 
% of the forward operator Phi and the adjoint operator
% Phi', where Phi is a representation a the measurement operator 
% whose information is provided in the struct 'meas_op'.
% =========================================================================
% Parameters
% ----------
% Returns
% ----------
% =========================================================================
% Author: Christian Kuemmerle, Johns Hopkins University, kuemmerle@jhu.edu,
% 2021.
if strcmp(meas_op.problem_type,'MultidimensionalScaling')
    if isfield(meas_op,'efficiencymode') && strcmp(meas_op.efficiencymode,'fast')
        %applyMeasopForward(applyMeasopBackward(y,meas_op),meas_op);
        %meas_op.PhiPhiT*y;%
        PhiPhiT_handle = @(z) MeasopPhiPhiT_MDS(z,meas_op);
        PhiPhiT_y = PhiPhiT_handle(y);
    else
        PhiPhiT_handle = @(z) meas_op.Phi*(meas_op.Phi'*z);
        PhiPhiT_y = meas_op.Phi*(meas_op.Phi'*y);
    end
elseif isfield(meas_op,'PhiPhiT')
    PhiPhiT_y = meas_op.PhiPhiT*y;
    PhiPhiT_handle = @(z) meas_op.PhiPhiT*z;
else
    PhiPhiT_handle = @(z) z;
    PhiPhiT_y = y;%error('This function has not been implemented for the provided measurement operator.')
end
end


function Phi_mat = MeasopPhiPhiT_MDS(z,meas_op)

% sps_plc = MeasopBackwardMDS_sparsepart(z,meas_op);
% n = size(sps_plc);
m = meas_op.m;
% Matdiag = diag(sps_plc);
% Phi_Mat_1 = Matdiag(meas_op.rowind);
% Phi_Mat_1 = Phi_Mat_1 + Matdiag(meas_op.colind); 
% Phi_Mat_1 = Phi_Mat_1 - sps_plc(meas_op.Omega);
% Phi_Mat_1 = Phi_Mat_1 - sps_plc(meas_op.OmegaT);
% Phi_Mat_1 = Phi_Mat_1 + sqrt(2./n(1)).*(z(m+meas_op.rowind)+z(m+meas_op.colind));
Phi_Mat_1 = meas_op.T_handle(z(1:m));% = @(y) apply_T_handle(y,meas_op,n);
avrage = mean(z(m+1:end));
Phi_Mat_2 = z(m+1:end)+avrage;
Phi_mat = [Phi_Mat_1;Phi_Mat_2];
end
