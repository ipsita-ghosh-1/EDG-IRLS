function y = proj_Oprange_tangspace(gam,meas_op,weight_op,varargin)
%proj_Omega_MC_smallsys This function projects a rank-2R (spaces (U,V))
%matrix onto its entries with indices in Omega
%       Input:      gamma = (R*(d1+d2-R))x1 vector. (for matrix completion)
U      = weight_op.U;
switch meas_op.problem_type
case {'MatrixCompletion','StructuredCompletion','CommunityDetection','MaxCut CD',...
        'MaxCut','MaxCut fixedval','MultidimensionalScaling'}
    [d1,R]  =size(U);
    if not(isfield(weight_op,'symmetricflag')) || ~weight_op.symmetricflag
        V      = weight_op.V;
        d2      =size(V,1);
    else
        V = U;
        d2 = d1;
    end
    D = max(d1,d2);
    mode   = varargin{1};
    if nargin == 6
       increase_antisymmetricweights = varargin{2}; 
    else
       increase_antisymmetricweights = 0;
    end
    switch meas_op.problem_type
        case {'MatrixCompletion','CommunityDetection','MaxCut CD','MaxCut','MaxCut fixedval',...
                'MultidimensionalScaling'}
            rowind = meas_op.rowind;
            colind = meas_op.colind;
            m=length(rowind);
        case {'StructuredCompletion'}
            if strcmp(meas_op.structure_type,'block-Hankel-Hermitian')
                rowind = meas_op.rowind_struc_sym_fat;
                colind = meas_op.colind_struc_sym_fat;
                accumind = meas_op.accumind_sym_fat;
            else
                rowind = meas_op.rowind_struc;
                colind = meas_op.colind_struc;
                accumind = meas_op.accumind;
            end
            m=length(rowind);
    end

    %%% Standard version for representations without U_{T_c} etc.:
    % Z2 = reshape(gam((R^2+1):(R*(d2+R))),[R,d2]);
    % Z3 = reshape(gam((R*(d2+R)+1):(R*(d2+d1+R))),[d1,R]);
    % if isreal(gam)
    %     gam_Omega = partXY((U*reshape(gam(1:R^2),[R,R])).',V',rowind,colind,m)';
    %     gam_Omega = gam_Omega+partXY(U.',Z2,rowind,colind,m)';
    %     gam_Omega = gam_Omega+partXY((Z3).',V',rowind,colind,m)';
    % else
    %     gam_Omega = partXY_cmplx((U*reshape(gam(1:R^2),[R,R])).',V',rowind,colind,m)';
    %     gam_Omega = gam_Omega+partXY_cmplx(U.',Z2,rowind,colind,m)';
    %     gam_Omega = gam_Omega+partXY_cmplx((Z3).',V',rowind,colind,m)';
    % end
    %%% Version where a (R x R) block is set 0 -> problem later on:
    % Z2 = [zeros(R,R),reshape(gam((R^2+1):(R*(d2))),[R,d2-R])];
    % Z3 = [zeros(R,R);reshape(gam((R*d2+1):(R*(d2+d1-R))),[d1-R,R])];
    %%% Version without 0 setting, but still representation with U_{T_c}:
    switch mode
        case 'rangespace_smallsys'
            M2 = reshape(gam((R^2+1):(R*(d2+R))),[R,d2]);
            M3 = reshape(gam((R*(d2+R)+1):(R*(d2+d1+R))),[d1,R]);
            Z2V = M2*V;
            UZ3 = U'*M3;

            if isreal(gam)
                y = partXY((U*reshape(gam(1:R^2),[R,R])).',V',rowind,colind,m)';
                y = y+partXY(U.',M2-Z2V*V',rowind,colind,m)';
                %   gam_Omega+partXY(U.',reshape(gam((R^2+1):(R*(d2))),[R,d2]),rowind,colind,m)';
                y = y+partXY((M3-U*UZ3).',V',rowind,colind,m)';
                %   gam_Omega = gam_Omega+partXY((reshape(gam((R*(R+d2)+1):end),[d1,R])).',V',rowind,colind,m)';
            else
                y = partXY_cmplx((U*reshape(gam(1:R^2),[R,R])).',V',rowind,colind,m)';
                y = y+partXY_cmplx(U.',M2-Z2V*V',rowind,colind,m)';
                y = y+partXY_cmplx((M3-U*UZ3).',V',rowind,colind,m)';
            end
        case 'tangspace'
            [~,M1,M2,M3] = get_Tk_matrices(gam,weight_op);
%             M1=reshape(gam(1:R^2),[R,R]);
%             M2=reshape(gam((R^2+1):(R*(d1+R))),[d1,R]);
%             M3=reshape(gam((R*(d1+R)+1):(R*(d2+d1+R))),[R,d2]);
            if increase_antisymmetricweights
                M1S_u = triu(M1);
                M1S_l = triu(M1,1)';
                M1S   = M1S_u+M1S_l;
                M1T_l = triu(M1,-1);
                M1T_u = -M1T_l';
                M1T   = M1T_u+M1T_l;
                Z1    = M1S+M1T;
                
                if d1 == D
                    Z2    = M2+[M3',zeros(d1,d1-d2)];
                    tmp=M2;
                    Z3    = M3+tmp(1:d2,:)';
                else
                    error('To be implemented.')
                end
                M1 = Z1;
                M2 = Z2;
                M3 = Z3;
            end
            switch meas_op.problem_type
                case {'MatrixCompletion','MaxCut CD','MaxCut','StructuredCompletion'}
                    if isreal(gam)
            %                 rowindit = uint32(rowind);
            %                 colindit = uint32(colind);
            %                 tmp = spmaskmult((U*M1+M2), V', rowindit, colindit); this
            %                 is not faster...
%                         UM2 = U'*M2;
%                         M3V = M3*V;
%                         y = partXY((U*M1+M2-U*UM2)',V',rowind,colind,m)';
%                         y = y+partXY(U.',M3-M3V*V',rowind,colind,m)';
%                         
                        y = partXY((U*M1+M2)',V',rowind,colind,m)';
                        y = y+partXY(U.',M3,rowind,colind,m)';
            %             y = partXY((U*reshape(gam(1:R^2),[R,R])+reshape(gam((R^2+1):(R*(d1+R))),[d1,R])).',V',rowind,colind,m)';
            %             y = y+partXY(U.',reshape(gam((R*(d1+R)+1):(R*(d2+d1+R))),[R,d2]),rowind,colind,m)';
                    else
                        y = partXY_cmplx((U*M1+M2).',V',rowind,colind,m)';
                        y = y+partXY_cmplx(U.',M3,rowind,colind,m)';
                    end
                    if strcmp(meas_op.problem_type,'StructuredCompletion')
                        ylong = y;
                        y = accumarray(accumind,ylong,[]);
                    end
                case {'CommunityDetection','MaxCut fixedval'}
                    if strcmp(meas_op.efficiencymode,'fast')
                        LeftMat1 = U*M1+M2;
                        y_1 = partXY(LeftMat1.',V',rowind,colind,d1)';
                        y_1 = y_1 + partXY(U.',M3,rowind,colind,d1)';
                        y_2 = sum(LeftMat1,1)*sum(V,1).';
                        y_2 = y_2 + sum(U,1)*sum(M3',1).';
                        y_2 = y_2./d1;
                        y = [y_1;y_2];
                    else
                        Mat1 = (U*M1+M2)*V';
                        Mat1 = Mat1 + U*M3;
                        y = meas_op.Phi*Mat1(:);
                    end
                case 'MultidimensionalScaling'
                    if strcmp(meas_op.efficiencymode,'fast')
                        LeftMat1 = U*M1+M2;
                        y_1 = partXY(LeftMat1.',V',rowind,rowind,m)';
                        y_1 = y_1 + partXY(LeftMat1.',V',colind,colind,m)';
                        y_1 = y_1 - partXY(LeftMat1.',V',rowind,colind,m)';
                        y_1 = y_1 - partXY(LeftMat1.',V',colind,rowind,m)';
                        y_1_b = LeftMat1*(sum(V,1).'./sqrt(2*d1));
                        y_1_b = y_1_b + (sum(LeftMat1,1)*(V.'./sqrt(2*d1))).';
                        
                        y_2 = partXY(U.',M3,rowind,rowind,m)';
                        y_2 = y_2 + partXY(U.',M3,colind,colind,m)';
                        y_2 = y_2 - partXY(U.',M3,rowind,colind,m)';
                        y_2 = y_2 - partXY(U.',M3,colind,rowind,m)';
                        y_2_b = U*(sum(M3',1).'./sqrt(2*d1));
                        y_2_b = y_2_b + (sum(U,1)*(M3./sqrt(2*d1))).';
                        y = [y_1+y_2;y_1_b+y_2_b];
                        %error('Fast efficiency mode not implemented yet.')
                    elseif strcmp(meas_op.efficiencymode,'buildfull')
%                         tmp0 = tangspace_to_matrixspace(gam,U,V);
%                         y = applyMeasopForward(tmp0,meas_op);
                        
                        LeftMat1 = U*M1+M2;
                        y_1 = partXY(LeftMat1.',V',rowind,rowind,m)';
                        y_1 = y_1 + partXY(LeftMat1.',V',colind,colind,m)';
                        y_1 = y_1 - partXY(LeftMat1.',V',rowind,colind,m)';
                        y_1 = y_1 - partXY(LeftMat1.',V',colind,rowind,m)';
                        
                        y_2 = partXY(U.',M3,rowind,rowind,m)';
                        y_2 = y_2 + partXY(U.',M3,colind,colind,m)';
                        y_2 = y_2 - partXY(U.',M3,rowind,colind,m)';
                        y_2 = y_2 - partXY(U.',M3,colind,rowind,m)';
                        y = y_1+y_2;
                    else
                        Mat1 = (U*M1+M2)*V';
                        Mat1 = Mat1 + U*M3;
                        y = meas_op.Phi*Mat1(:);
                    end
            end
    end
    
case 'PhaseRetrieval'
    % To Do (Apr 14, 2021): Adapt this to new API with meas_op as second
    % argument.
    A       = varargin{1};
    AU      = varargin{2};
    [~,R] = size(AU);
    n   = size(U,1);
%   X1= reshape(gam(1:R^2),[R,R]);
%   X2= related to reshape(gam((R^2+1):end),[n,R]);
    handle_A=functions(A);
    if contains(handle_A.function,'fourierMeasurementOperator')
        Aop= @(X) cell2mat(cellfun(A,num2cell(X,1),'UniformOutput',false));
        tmp=Aop(reshape(conj(gam((R^2+1):end)),[n,R]));
    else
        tmp=A(reshape(conj(gam((R^2+1):end)),[n,R]));
    end
    y = sum((AU*reshape(gam(1:R^2),[R,R])).*conj(AU),2);
    y = y + sqrt(2).*real(sum(AU.*tmp,2));
%     y = y + sqrt(2).*real(sum(AU.*(Aop(reshape(conj(gam((R^2+1):end)),[n,R]))),2));
otherwise
    error('proj_Oprange_tangspace.m not yet implemented for this problem.')
end
    
end

