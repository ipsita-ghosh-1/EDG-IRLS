function Phi_Mat = MeasopPhiPhiT_MDS(z,meas_op)
m = meas_op.m;
if isfield(meas_op,'T_handle')
    Phi_Mat_1 = meas_op.T_handle(z(1:m));% = @(y) apply_T_handle(y,meas_op,n);
    avrage = mean(z(m+1:end));
    Phi_Mat_2 = z(m+1:end)+avrage;
    Phi_Mat = [Phi_Mat_1;Phi_Mat_2];
else
    sps_plc = MeasopBackwardMDS_sparsepart(z,meas_op);
    % n = size(sps_plc,1);
    Matdiag = diag(sps_plc);        
    Phi_Mat_1 = Matdiag(meas_op.rowind);
    Phi_Mat_1 = Phi_Mat_1 + Matdiag(meas_op.colind); 
    Phi_Mat_1 = Phi_Mat_1 - sps_plc(meas_op.Omega);
    Phi_Mat_1 = Phi_Mat_1 - sps_plc(meas_op.OmegaT);
    %Phi_Mat_1 = Phi_Mat_1 + (z(m+meas_op.rowind)+z(m+meas_op.colind)); % sqrt(2./n).*
    avrage = mean(z(m+1:end));
    Phi_Mat_2 = z(m+1:end)+avrage;
    Phi_Mat = [Phi_Mat_1;Phi_Mat_2];
end
end