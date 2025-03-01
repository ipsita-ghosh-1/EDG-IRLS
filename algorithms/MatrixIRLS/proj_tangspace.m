function gam_T = proj_tangspace(gam,problem,weight_op_0,weight_op)
%proj_tangspace This function projects the matrices from the tangent space
%T0 into the tangent space T
U0 = weight_op_0.U;
U  = weight_op.U;
if not(isfield(weight_op,'symmetricflag')) || ~weight_op.symmetricflag
    V0 = weight_op_0.V;
    V  = weight_op.V;
end

Utild_U = U'*U0;
if strcmp(problem,'MatrixCompletion')
    [d1,r_T0]=size(U0);
    [~,r_T]=size(U);
    if isfield(weight_op,'symmetricflag') && weight_op.symmetricflag
        X1_upper = zeros(r_T0,r_T0);
        gam1 = gam(1:r_T0*(r_T0+1)/2);
        X1_upper(weight_op_0.symInd) = gam1;
        X1_lower = triu(X1_upper,1);
        X1 = X1_upper + X1_lower';
        X2 = reshape(gam(r_T0*(r_T0+1)/2+1:end),d1,r_T0);
        gam_T=zeros(r_T*((r_T+1)/2+d1),1);
        gam_T_1 = Utild_U*(X1*Utild_U'+X2'*U)+U'*X2*Utild_U';
        gam_T_2 = (U0-U*Utild_U)*(X2'*U)+(U0*X1+X2-U*Utild_U*X1-U*(U'*X2))*Utild_U';
        gam_T(1:r_T*(r_T+1)/2)     = gam_T_1(weight_op.symInd);%reshape(gam_T_1,[r_T^2,1]);
        gam_T(r_T*(r_T+1)/2+1:end) = gam_T_2(:);%reshape(gam_T_2,[r_T*d1,1]);
    else
        d2    =size(V0,1);
        V_Vtild = V0'*V;
        gam_T=zeros(r_T*(d1+d2+r_T),1); %gam;
        X1= reshape(gam(1:r_T0^2),[r_T0,r_T0]);
        X2= reshape(gam((r_T0^2+1):(r_T0*(d1+r_T0))),[d1,r_T0]);
        X3= reshape(gam((r_T0*(d1+r_T0)+1):(r_T0*(d2+d1+r_T0))),[r_T0,d2]);
        gam_T_1=Utild_U*(X1*V_Vtild+X3*V)+U'*X2*V_Vtild;
        gam_T_2=(U0-U*Utild_U)*(X3*V)+(U0*X1+X2-U*Utild_U*X1-U*(U'*X2))*V_Vtild;
        gam_T_3=Utild_U*(X1*V0'+X3-X1*V_Vtild*V'-(X3*V)*V')+...
            U'*X2*(V0'-V_Vtild*V');
        gam_T(1:r_T^2)                = reshape(gam_T_1,[r_T^2,1]);
        gam_T((r_T^2+1):(r_T*(r_T+d1))) = reshape(gam_T_2,[r_T*d1,1]);
        gam_T((r_T*(r_T+d1)+1):end)     = reshape(gam_T_3,[r_T*d2,1]);
    end
    
elseif strcmp(problem,'PhaseRetrieval')
    [n,r_T0]=size(U0);
    [~,r_T]=size(U);
    gam_T=zeros(r_T*(n+r_T),1);%gam;
    X1= reshape(gam(1:r_T0^2),[r_T0,r_T0]);
    X2= reshape(gam((r_T0^2+1):(r_T0*(r_T0+n))),[n,r_T0]);
    X2U=X2'*U;
    gam_T_1=Utild_U*(X1*Utild_U'+X2U./sqrt(2))+(U'*X2*Utild_U')./sqrt(2);
    gam_T_2=(U0-U*Utild_U)*X2U+(X2-U*(U'*X2)+sqrt(2).*(U0*X1-U*Utild_U*X1))*Utild_U';
    gam_T(1:r_T^2)                = reshape(gam_T_1,[r_T^2,1]);
    gam_T((r_T^2+1):(r_T*(r_T+n)))   = reshape(gam_T_2,[r_T*n,1]);
else
    error('proj_tangspace.m not yet implemented for this problem.')
end

end

%(U0*X1-U*Utild_U*X1)*Utild_U'