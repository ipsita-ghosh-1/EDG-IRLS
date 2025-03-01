function PhiPhiT_inv_y = applyMeasopForwardAdjointInv(y,meas_op,varargin)
%applyMeasopForwardAdjointInv This function computes the image
% (Phi *Phi')^(-1)(y) of the vector 'y' with with respect to the inverse of
% the concatenation of the forward operator Phi and the adjoint operator
% Phi', where Phi is a representation a the measurement operator 
% whose information is provided in the struct 'meas_op'.
% =========================================================================
% Parameters
% ----------
% Returns
% ----------
% =========================================================================
global PhiPhiT_inv_y_c
if isfield(meas_op,'efficiencymode') && strcmp(meas_op.efficiencymode,'fast')   
    if nargin > 2 && strcmp(varargin{1},'V')
        PhiPhiT_inv_y = meas_op.PhiVPhiVT_handle(y);
    else
        hand_PhiPhiT = @(z) MeasopPhiPhiT_MDS(z,meas_op);%MeasopPhiPhiT_MDS(z,meas_op)
    %         hand_PhiPhiT = @(z) applyMeasopForward(applyMeasopBackward(z,meas_op),meas_op);%meas_op.PhiPhiT*z;%
    %         error('Not implemented yet.')
        prec_V = @(z) meas_op.PhiVPhiVT_handle(z);%meas_op.PhiVPhiVT*z;
        tol_Adj_CG = 1e-14;
        N0_Adj_CG = 10; %prec_V;
        [PhiPhiT_inv_y,flag_AdjInv,relres_AdjInv,N_inner_AdjInv] = pcg_1(hand_PhiPhiT,y,tol_Adj_CG,N0_Adj_CG,prec_V,[],PhiPhiT_inv_y_c);%1e-6,50);
    end
%     hand_PhiPhiT = @(z) MeasopPhiPhiT_MDS(z,meas_op);%MeasopPhiPhiT_MDS(z,meas_op)
% %         hand_PhiPhiT = @(z) applyMeasopForward(applyMeasopBackward(z,meas_op),meas_op);%meas_op.PhiPhiT*z;%
% %         error('Not implemented yet.')
%     prec_V = @(z) meas_op.PhiVPhiVT*z;
%     tol_Adj_CG = 1e-14;
%     N0_Adj_CG = 100; %prec_V;
%     PhiPhiT_inv_y = meas_op.PhiVPhiVT*y;
%     PhiPhiT_inv_y = meas_op.PhiVPhiVT*y;
    % [PhiPhiT_inv_y,flag_AdjInv,relres_AdjInv,N_inner_AdjInv] = pcg_1(hand_PhiPhiT,y,tol_Adj_CG,N0_Adj_CG,prec_V,[],PhiPhiT_inv_y_c);%1e-6,50);
    % fprintf(['N_inner_AdjI: %5d   ', ...
    %         'relres_AdjI: %.3e  flag: %1d\n'], ...
    %         N_inner_AdjInv,relres_AdjInv,flag_AdjInv)
    % PhiPhiT_inv_y_c = PhiPhiT_inv_y;
    % [PhiPhiT_inv_y_t,flag_AdjInv_t,relres_AdjInv_t,N_inner_AdjInv_t] = pcg_1(hand_PhiPhiT,y,tol_Adj_CG,N0_Adj_CG);%1e-6,50);
    % fprintf(['N_inner_AdjI_t: %5d   ', ...
    %         'relres_AdjI_t: %.3e  flag_t: %1d\n'], ...
    %         N_inner_AdjInv_t,relres_AdjInv_t,flag_AdjInv_t)
else
    if nargin > 2
        if strcmp(varargin{1},'V')
            PhiPhiT_inv_y = meas_op.PhiVPhiVT*y;
        % PhiPhiT_inv_y = meas_op.PhiPhiT\y;
        else
            PhiPhiT_inv_y = meas_op.PhiPhiT\y;
        end
    else
        PhiPhiT_inv_y = meas_op.PhiPhiT\y;
    end
    % PhiPhiT_inv_y = y;
end
% hand_PhiPhiT = @(z) MeasopPhiPhiT_MDS(z,meas_op);%MeasopPhiPhiT_MDS(z,meas_op)
% hand_PhiPhiT_2 = @(z) MeasopPhiPhiT_MDS(z,meas_op);
end

