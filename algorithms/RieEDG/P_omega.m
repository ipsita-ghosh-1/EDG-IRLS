function X = P_omega(X_in,samples)

%this operator takes in a distance matrix and outputs a masked distance
%matrix
%this is currently written for symmetric samples. plays nicely with
%Construct_Samples, but if you're interested in an unsymmetric sampling
%pattern you'll need to change this!

    n = size(X_in,1);
    m = size(X_in,2);
    X = zeros(n,m);
    IJ = sub2ind([n,m],samples(:,1),samples(:,2));
    X(IJ) = X_in(IJ);        

return