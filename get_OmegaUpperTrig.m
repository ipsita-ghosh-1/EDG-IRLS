function [OmegaUpperTrig,OmegaMatrix] = get_OmegaUpperTrig(n,OmegaLin)
%OmegaUpperTrig Creates an upper triangluar sparse matrix with 1's at
%indices of a sampling set Omega \subset \{(i,j):j > i\} provided by a
%vector OmegaLin with linear indices of Omega (between 1 and n*(n-1)/2).
Hermitian_flag = true;
sampling_vector_withdiag_all = zeros(n*(n+1)/2,1);
ind_c = -n;
for i = 1:n
    ind_c = ind_c+(n-i+2);
    diag_ind(i) = ind_c;
end
nondiag_ind = setdiff([1:n*(n+1)/2],diag_ind);
sampling_vector_withdiag_all = zeros(n*(n+1)/2,1);
sampling_vector_withdiag_all(nondiag_ind(OmegaLin)) = 1;
Weight_all = 2.*get_matrix_from_vector(sampling_vector_withdiag_all,n,Hermitian_flag);
Weight_all_upper = triu(Weight_all,1);
sampling_vector_withdiag = zeros(n*(n+1)/2,1);
sampling_vector_withdiag(nondiag_ind(OmegaLin)) = 1;
Weight = 2.*get_matrix_from_vector(sampling_vector_withdiag,...
    n,Hermitian_flag);
Weight_upper = triu(Weight,1);%triu(Weight,1);
OmegaMatrix = sparse(Weight);
OmegaUpperTrig = sparse(Weight_upper);

end