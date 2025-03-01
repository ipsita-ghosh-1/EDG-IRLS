    function [gam,sps_plc] = proj_tangspace_from_Oprange(y,meas_op,weight_op,varargin)
%proj_tangentspace_from_Omega This function projects a m-sparse matrix y
%with values in Omega
%into the tangent space
% T:= {U*Z1' + Z2*V' ; Z1 \in \R^(d1xr), Z2 \in \R^(d2xr)}. The map can be
% written as
%P_T(P_Omega'(y)) = U*U'*P_Omega'(y) + *P_Omega'(y)*V*V' - U*U'**P_Omega'(y)*V*V'.
switch meas_op.problem_type
case {'MatrixCompletion','CommunityDetection',...
        'MaxCut CD','MaxCut','MaxCut fixedval',...
        'MultidimensionalScaling','StructuredCompletion'}
    U      = weight_op.U;
    [d1,R]=size(U);
    if not(isfield(weight_op,'symmetricflag')) || ~weight_op.symmetricflag
        V      = weight_op.V;
        d2      =size(V,1);
    else
        d2 = d1;
    end
    D = max(d1,d2);
%     U       = varargin{1};
%     V       = varargin{2};
%     meas_op = varargin{3};
    mode    = varargin{1};
    if nargin == 6
       increase_antisymmetricweights = varargin{2}; 
    else
       increase_antisymmetricweights = 0;
    end
    switch meas_op.problem_type
        case {'MatrixCompletion','StructuredCompletion','MaxCut CD','MaxCut'}
            if strcmp(meas_op.problem_type,'StructuredCompletion')
                if strcmp(meas_op.structure_type,'block-Hankel-Hermitian')
                    ytmp = transform_y_ylong_structure(y,meas_op,'Hermitian');
                    yused = ytmp(meas_op.permind_sym);
                    sps_plc_triang = meas_op.sps_plc_struc_sym;
                    m = length(yused);
                    setSval(sps_plc_triang,yused,m);
                    sps_plc_sym = (sps_plc_triang+sps_plc_triang')./2;%max(sps_plc_triang,sps_plc_triang');
                    if size(meas_op.sps_plc_struc,2) > size(meas_op.sps_plc_struc,1)
                        sps_plc = sps_plc_sym(1:size(meas_op.sps_plc_struc,1),:);
                        
                    elseif size(meas_op.sps_plc_struc,2) < size(meas_op.sps_plc_struc,1)
                        sps_plc = sps_plc_sym(:,1:size(meas_op.sps_plc_struc,2));
                    else
                        sps_plc = sps_plc_sym;
                    end
%                     sps_plc = sps_plc_triang();
                else
                    ytmp = transform_y_ylong_structure(y,meas_op);
                    yused = ytmp(meas_op.permind);
                    sps_plc = meas_op.sps_plc_struc;
                    m = length(yused);
                    setSval(sps_plc,yused,m);
                end
            else
                yused = y;
                sps_plc = meas_op.sps_plc;
                m = length(yused);
                setSval(sps_plc,yused,m);
            end
        case {'CommunityDetection','MaxCut fixedval'}
            if d1 == d2
                n = d1;
            else
                error('Matrix not square.')
            end
            if strcmp(meas_op.efficiencymode,'fast')
                sps_plc_fast = meas_op.sps_plc_fast;
                setSval(sps_plc_fast,y(1:n),n);
            else
                sps_plc = reshape(meas_op.Phi'*(y),[n,n]);
            end
        case 'MultidimensionalScaling'
            n = d1;
            if strcmp(meas_op.efficiencymode,'fast')
                sps_plc_fast = MeasopBackwardMDS_sparsepart(y,meas_op);
            elseif strcmp(meas_op.efficiencymode,'buildfull')
%                 tmp0 = tangspace_to_matrixspace(X_Tn_c,U,V);
% %             tmp0 = meas_op.Samp_opmat*(reshape(meas_op.kappainv(meas_op.kappainv(meas_op.kappa(tmp0))),size(U,1)*size(U,1),1));
% %             tmp0 = (-2).*meas_op.Samp_opmat*(reshape(meas_op.kappainv(meas_op.kappa(tmp0)),size(U,1)*size(U,1),1));
%                 s = applyMeasopForward(tmp0,meas_op);
                
                sps_plc = reshape(applyMeasopBackward(y,meas_op),d1,d1);
%                         y = matrixspace_to_tangspace(reshape(tmp2,d1,d1),U,V);
            else
                sps_plc = reshape(meas_op.Phi'*(y),[n,n]);
            end
    end
%     vals=nonzeros(sps_plc);
%     UPhst_Yval_tst = multfullsparse(U', vals, sps_plc);
%     Phst_Yval_tst = multsparsefull(vals, V, sps_plc); (this is not
%     faster)
    if strcmp(mode,'rangespace_smallsys')
        UPhst_Yval = U'*sps_plc;
        Phst_Yval = sps_plc*V;
    %%% Standard version for representations without U_{T_c} etc.:
    % y_T(1:R^2)=reshape(-UPhst_Yval*V,[R^2,1]);
    %%% Version for representation with U_{T_c} etc.:
        gam =zeros(R*(d1+d2+R),1);
        gam(1:R^2)=reshape(UPhst_Yval*V,[R^2,1]);
        gam((R^2+1):(R*(R+d2)))=reshape(UPhst_Yval,[R*d2,1]);
        gam((R*(R+d2)+1):end)=reshape(Phst_Yval,[d1*R,1]);
    elseif strcmp(mode,'tangspace')
        if increase_antisymmetricweights
            UPhst_Yval = U'*sps_plc;
            Phst_Yval = sps_plc*V;
            M1= UPhst_Yval*V;
            gam =zeros(R*(d1+d2+R),1);
            M1S = (M1+M1')./2;
            M1T = (M1-M1')./2;
            M1_upper  = triu(M1S);
            M1_lower = tril(M1T,-1);
            M1  = M1_upper+M1_lower;
            if d1 == D
                M2=(Phst_Yval-U*M1+[UPhst_Yval'-V*M1' ;zeros(d1-d2,d1)])./2;
                tmp=Phst_Yval'-M1'*U';
                M3=(UPhst_Yval-M1*V'-tmp(:,1:d2))./2;
                gam(1:R^2)              = reshape(M1,[R^2,1]);
                gam((R^2+1):(R*(R+d1))) = reshape(M2,[R*d1,1]);
                gam((R*(R+d1)+1):end)   = reshape(M3,[R*d2,1]);
            else
               error('To be implemented.') 
            end
        else
            if strcmp(meas_op.problem_type,'MatrixCompletion') || strcmp(meas_op.problem_type,'MaxCut CD') ...
                    || strcmp(meas_op.problem_type,'MaxCut') || not(strcmp(meas_op.efficiencymode,'fast')) ...
                    || strcmp(meas_op.problem_type,'StructuredCompletion')
                gam = matrixspace_to_tangspace(sps_plc,weight_op);
            else
                gam = matrixspace_to_tangspace(sps_plc_fast,weight_op);
                if strcmp(meas_op.problem_type,'CommunityDetection')
                    gam = gam + matrixspace_to_tangspace({y(n+1).*ones(n,1)./n,ones(n,1)},weight_op);
                elseif strcmp(meas_op.problem_type,'MaxCut fixedval')
                    gam = gam + matrixspace_to_tangspace((y(n+1)/meas_op.normA).*meas_op.A,weight_op);
                elseif strcmp(meas_op.problem_type,'MultidimensionalScaling')
                    gam = gam + matrixspace_to_tangspace({ones(n,1),y(meas_op.m+1:end)},weight_op)./sqrt(2*n);
                    gam = gam + matrixspace_to_tangspace({y(meas_op.m+1:end),ones(n,1)},weight_op)./sqrt(2*n);
                end
            end
%             gam(1:R^2)              = reshape(M1,[R^2,1]);
%             gam((R^2+1):(R*(R+d1))) = reshape(Phst_Yval-U*M1,[R*d1,1]);
%             gam((R*(R+d1)+1):end)   = reshape(UPhst_Yval-M1*V',[R*d2,1]);
        end
    end
    
case 'RobustPCA'
    U       = varargin{1};
    V       = varargin{2};
%     [d1,R]=size(U);
%     d2    =size(V,1);
    [d1,d2] = size(y);
    if isequal(size(y),[d1,d2])
        Y=y;
    else
        Y=reshape(y,[d1,d2]);
    end
    weight_op.U = U;
    weight_op.V = V;
    gam = matrixspace_to_tangspace(Y,weight_op);
%     R = size(U,2);
%     gam =zeros(R*(d1+d2+R),1);
%     UPhst_Yval= U'*Y;
%     Phst_Yval = Y*V;
%     M1= UPhst_Yval*V;
%     gam(1:R^2)              = reshape(M1,[R^2,1]);
%     gam((R^2+1):(R*(R+d1))) = reshape(Phst_Yval-U*M1,[R*d1,1]);
%     gam((R*(R+d1)+1):end)   = reshape(UPhst_Yval-M1*V',[R*d2,1]);
case 'PhaseRetrieval'
    % To Do (Apr 14, 2021): Adapt this to new API with meas_op as second
    % argument.
    At = varargin{1};
    AU = varargin{2};
    [~,R] = size(AU);
	n = size(U,1);
    
    yAU=y.*AU;
    X1=AU'*(yAU);
    handle_At=functions(At);
    if contains(handle_At.function,'transposeOperator')
        At_op= @(X) cell2mat(cellfun(At,num2cell(X,1),'UniformOutput',false));
        AtyAU=conj(At_op(conj(yAU)));
        X2=sqrt(2).*(AtyAU-U*(AU'*yAU));
    else
        X2=sqrt(2).*(conj(At(conj(yAU)))-U*(AU'*yAU));
    end
%     X2=sqrt(2).*(conj(At(conj(yAU)))-U*(AU'*yAU)); %At_op(conj(yAU)
    gam = zeros(R*(n+R),1);
    gam(1:R^2)=reshape(X1,[R^2,1]);
    gam((R^2+1):(R*(R+n)))=reshape(X2,[R*n,1]);
    
otherwise
    error('proj_tangspace_from_Oprange.m not yet implemented for this problem.')
end

end

