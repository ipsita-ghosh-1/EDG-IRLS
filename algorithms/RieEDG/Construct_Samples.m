function [samples,Weight] = Construct_Samples(n,m,L,rate,sampling_mode)
%n is rows, m is columns. works for fat matrices. not hard to change for
%tall but code is written for fat.

if sampling_mode == 0
    Weight= rand(n,m);
    Weight(Weight>1-rate)=1;
    Weight(Weight<1)=0;
    Weight(Weight>0)=1;
    for l=1:n
        %diagonal entries are known i.e. D_ii = 0
        Weight(l,l)= 1;
        for p=l+1:n
          %make the weight matrix symmetric
          Weight(l,p)= Weight(p,l);
        end
    end
elseif sampling_mode == 1
    samps = ceil(L*rate);
    inds = randsample(L,samps,false);
    inds = sort(inds,'ascend');
    I = zeros(samps,1);
    J = zeros(samps,1);
    row = 1;
    thresh = n-row;
    for q=1:samps
        while inds(q) > thresh
            row = row+1;
            thresh = thresh + n - row;
        end

        I(q) = row;
        J(q) = inds(q) - (thresh - n);
    end
    A = sub2ind([n m],I,J);
    %for speed purposes, this does not need to be formed. I
    %form it for visualization/other testing purposes if
    %necessary, since it isn't that expensive in the grand
    %scheme of the algorithm.
    Weight = zeros(n,m);
    Weight(A) = 1;
    for l=1:n
        %diagonal entries are known i.e. D_ii = 0
        Weight(l,l)= 1;
        for p=l+1:n
          %make the weight matrix symmetric
          Weight(p,l)= Weight(l,p);
        end
    end
else
    disp 'choose 0 or 1 for sampling_mode'
end

[I,J] = find(Weight==1);
samples = [I,J];

return






