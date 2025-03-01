function [X_c,X_Tk_1,X_Tk_2,X_Tk_3] = get_Tk_matrices(gam,weight_op)%d1,d2,r)
%get_Tk_matrices Summary of this function goes here
%   Detailed explanation goes here
[d1,r] = size(weight_op.U);
if ~isfield(weight_op,'symmetricflag') || ~weight_op.symmetricflag
    d2 = size(weight_op.V,1);
else
    d2 = d1;
end
if isfield(weight_op,'symmetricflag') && weight_op.symmetricflag
    X_Tk_1_upper = zeros(r,r);
    gam1 = gam(1:r*(r+1)/2);
    X_Tk_1_upper(weight_op.symInd) = gam1;
    X_Tk_1_lower = triu(X_Tk_1_upper,1);
    X_Tk_1 = X_Tk_1_upper + X_Tk_1_lower';
    X_Tk_2 = reshape(gam(r*(r+1)/2+1:end),d1,r);
    X_Tk_3 = X_Tk_2';%[];
    X_c.Gam1 = X_Tk_1;
    X_c.Gam2 = X_Tk_3;
    X_c.Gam3 = X_Tk_2;
else
    if not(length(gam) == r*(d1+d2+r))
        error('Error in the dimensionality.')
    end
    X_Tk_1=reshape(gam(1:r^2),[r,r]);
    % X_Tk_2=reshape(gam((r^2+1):(r*(r+d1))),[d1,r]);
    X_Tk_2=reshape(gam((r^2+1):(r*(r+d1))),[d1,r]);
    % X_Tk_2=X_Tk_2;
    X_Tk_3=reshape(gam((r*(r+d1)+1):end),[r,d2]);

    X_c.Gam1 = X_Tk_1;
    X_c.Gam2 = X_Tk_3;
    X_c.Gam3 = X_Tk_2;
end
end

