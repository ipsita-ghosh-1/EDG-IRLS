function gam_out = apply_weight_vec(gam,weight_vec,...
    increase_antisymmetricweights,varargin)
%apply_weight_vec This function just applies an elementwise multiplication
% of gam with the elements of weight_vec in the case of most 'type_mean'.
% However, for cases where weights are chosen according to a symmetric-skew-
% symmetric splitting as type_mean = 'Hessianupper_2', the multiplication
% is not diagaonal, but takes this splitting into account.
if increase_antisymmetricweights
    d1 = varargin{1};
    d2 = varargin{2};
    r=round((-(d1+d2)+sqrt((d1+d2)^2+4*length(gam)))./2);
    if not(length(gam) == r*(d1+d2+r))
        error('Error in the dimensionality.')
    end
    
    gam_Tkanti = get_antisymmTkbasis_from_Tk(gam,d1,d2,r);
    gam_Tkanti = weight_vec.*gam_Tkanti;
    gam_weighted = get_Tkbasis_from_antisymmTk(gam_Tkanti,d1,d2,r);
%     gam_out=proj_tangspace(gam_weighted,'MatrixCompletion',U,V,U,V);
    gam_out=gam_weighted;
% %     case {'Hessianupper_2','Hessianupper_Ribeiro'}
%         d=min(d1,d2);
%         D=max(d1,d2);
%         R=round((-(d1+d2)+sqrt((d1+d2)^2+4*length(gam)))./2);
%         gam_weighted=zeros(R*(d1+d2+R),1);
%         Z1=reshape(gam(1:R^2),[R,R]);
%         Z2=reshape(gam((R^2+1):(R*(R+d1))),[d1,R]);
%         Z3=reshape(gam((R*(R+d1)+1):end),[R,d2]);
%         
%         [Z1_S] = symmpart(Z1);
%         [Z1_A] = skewsymmpart(Z1);
%         
%         Z23_S = zeros(D,R);
%         Z23_A = zeros(d,R);
%         Z23_S(1:d,:)=0.5.*(Z2(1:d,:)+Z3(:,1:d)');
%         Z23_A(1:d,:)=0.5.*(Z2(1:d,:)-Z3(:,1:d)');
%         if d1 == D
%             Z23_S((d+1):D,:)=Z2((d+1):D,:);
%             H_symm_23=reshape(weight_vec((R^2+1):(R*(d1+R))),[d1,R]);
%             H_sksy_23=reshape(weight_vec((R*(d1+R)+1):(R*(d2+d1+R))),[R,d2])';
%         else
%             Z23_S((d+1):D,:)=Z3(:,(d+1):D)';
%             H_sksy_23=reshape(weight_vec((R^2+1):(R*(d1+R))),[d1,R]);
%             H_symm_23=reshape(weight_vec((R*(d1+R)+1):(R*(d2+d1+R))),[R,d2])';
%         end
%         H_symm_1=reshape(weight_vec(1:R^2),[R,R]);
%         H_sksy_1=H_symm_1;
%         for i=1:R
%             for j=1:R
%                 if j > i
%                     H_sksy_1(i,j)=H_sksy_1(j,i);
%                     H_symm_1(j,i)=H_symm_1(i,j);
%                 end
%             end
%         end
%         gam_weighted(1:R^2)=reshape(H_symm_1.*Z1_S+H_sksy_1.*Z1_A,[R^2,1]);
%         Z23_S=H_symm_23.*Z23_S;
%         Z23_A=H_sksy_23.*Z23_A;
%         if d1 == D
%            Z2=[Z23_S(1:d,:)+Z23_A;Z23_S((d+1):D,:)];
%            Z3=(Z23_S(1:d,:)-Z23_A)';
%         else
%            Z2=Z23_S(1:d,:)+Z23_A;
%            Z3=[Z23_S(1:d,:)-Z23_A;Z23_S((d+1):D,:)]';
%         end
%     gam_weighted((R^2+1):(R*(d1+R)))=...
%     reshape(Z2,[d1*R,1]);
%     gam_weighted((R*(d1+R)+1):(R*(d2+d1+R)))  =...
%     reshape(Z3,[d2*R,1]);
else
        gam_out = weight_vec.*gam;
end
end

