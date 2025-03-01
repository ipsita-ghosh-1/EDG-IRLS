function weight_vec = get_weight_vec(d1,d2,H,dH,tangent_para,...
    weight_op,increase_antisymmetricweights)
%get_weight_vec This function returns a weight vector "weight_vec" from the
% information of H and dH. In particular, it covers the case of weight operators 
% with symmetric-skewsymmetric splitting if increase_antisymmetricweights
% == 1.
R=size(H{1},1);
d=min(d1,d2);
D=max(d1,d2);
if isfield(weight_op,'symmetricflag')
    symmetricflag = weight_op.symmetricflag;
else
    symmetricflag = false;
end
if strcmp(tangent_para,'intrinsic')
    dim_tangent = R*(d1+d2-R);
    if weight_op.symmetricflag
       error('This needs to be updated.') 
    end
else
    if symmetricflag
        dim_core = R*((R+1)/2);
        dim_tangent = R*((R+1)/2+d1);
    else
        dim_core = R^2;
        dim_tangent = R*(d1+d2+R);
    end
end
weight_vec =zeros(dim_tangent,1);

if symmetricflag
    weight_vec(1:dim_core)        = H{1}(weight_op.symInd);
%     gam(r*(r+1)/2+1:end)    = reshape((UtX'-U*M1),[r*d1,1]);
else
    weight_vec(1:dim_core)=H{1}(:);
%     gam(1:r^2)              = reshape(M1,[r^2,1]);
%     gam((r^2+1):(r*(r+d1))) = reshape((XV-U*M1),[r*d1,1]);
%     gam((r*(r+d1)+1):end)   = reshape((UtX-M1*V'),[r*d2,1]);
end

if increase_antisymmetricweights
    if d1 == D
        if strcmp(tangent_para,'intrinsic')
            weight_vec((dim_core+1):(R*(d1)))=...
            reshape(kron(dH{1},ones(1,d1-R)).',[(d1-R)*R,1]);
            weight_vec((R*d1+1):dim_tangent)  =...
            reshape(kron(dH{2},ones(1,d2-R)),[(d2-R)*R,1]);
        else
            weight_vec((dim_core+1):(dim_core+R*d1))=...
            reshape(kron(dH{1},ones(1,d1)).',[d1*R,1]);
            if ~symmetricflag
                weight_vec((dim_core+R*d1+1):dim_tangent)  =...
                reshape(kron(dH{2},ones(1,d2)),[d2*R,1]);
            end
        end
    else
        if strcmp(tangent_para,'intrinsic')
            weight_vec((dim_core+1):R*d1)=...
            reshape(kron(dH{2},ones(1,d1-R)).',[(d1-R)*R,1]);
            weight_vec((R*d1+1):dim_tangent)  =...
            reshape(kron(dH{1},ones(1,d2-R)),[(d2-R)*R,1]);
        else
            weight_vec((dim_core+1):(dim_core+R*d1))=...
            reshape(kron(dH{2},ones(1,d1)).',[d1*R,1]);
            if ~symmetricflag
                weight_vec((dim_core+R*d1+1):dim_tangent)  =...
                reshape(kron(dH{1},ones(1,d2)),[d2*R,1]);
            end
        end
    end
else
    if strcmp(tangent_para,'intrinsic')
        weight_vec((dim_core+1):R*d1)=...
        reshape(kron(dH{2},ones(1,d1-R)).',[(d1-R)*R,1]);
        weight_vec((R*d1+1):dim_tangent)  =...
        reshape(kron(dH{1},ones(1,d2-R)),[(d2-R)*R,1]);
    else
        weight_vec((dim_core+1):(dim_core+R*d1))=...
        reshape(kron(dH{2},ones(1,d1)).',[d1*R,1]);
        if ~symmetricflag
            weight_vec((dim_core+R*d1+1):dim_tangent)  =...
            reshape(kron(dH{1},ones(1,d2)),[d2*R,1]);
        end
    end
end 
end

