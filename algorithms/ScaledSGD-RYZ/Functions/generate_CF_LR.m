function [spdata, n_movie] = generate_CF_LR(data, train_size, test_size, UseParall)
% Generate samples of ground truth item-item matrix from MovieLens dataset
if nargin < 4; UseParall = true; end
worker = 10;
m = train_size + test_size;
million = 1e6;

% Form the user-item matrix
n_user = max(data(:,1)); n_movie = max(data(:,2));
UserItem = sparse(data(:,1),data(:,2),data(:,3),n_user,n_movie);
clear data

% Delete items with no users rating and update the number of movies
ItemKeep = full(sum(logical(UserItem))) > 0;
UserItem = UserItem(:,ItemKeep);
ItemNorm = full(sqrt(sum(UserItem.^2,1)));
n_movie = numel(ItemNorm);
clear ItemKeep

% Generate random item-item content
% Note that we sample 3 times more because some of the samples might have
% Mij = Mik, which we will delete it later.
sample = 3*m; sample = sample - mod(sample,worker); batch = sample/worker;
[i,j,k] = ind2sub(n_movie*[1,1,1], randperm(n_movie^3, sample));
i = reshape(i,batch,worker); j = reshape(j,batch,worker); k = reshape(k,batch,worker);
Mij = zeros(batch,worker); Mik = zeros(batch,worker);
if UseParall
    parfor pdx = 1:worker
        w = 0;
        ipar = i(:,pdx); jpar = j(:,pdx); kpar = k(:,pdx);
        UI = UserItem; IN = ItemNorm;
        Mijpar = zeros(batch,1); Mikpar = zeros(batch,1);
        for idx = 1:batch
            ii = ipar(idx); jj = jpar(idx); kk = kpar(idx);
            Mijpar(idx) = UI(:,ii)'*UI(:,jj)/(IN(ii)*IN(jj));
            Mikpar(idx) = UI(:,ii)'*UI(:,kk)/(IN(ii)*IN(kk));
            if pdx == worker && mod(idx*worker,million)==0
                w = fprintf([repmat('\b',1,w),'sample: %dM/%dM\n'],idx*worker/million,sample/million) - w; 
            end
        end
        Mij(:,pdx) = Mijpar; Mik(:,pdx) = Mikpar;
    end
else
    w = 0;
    for pdx = 1:worker
        ipar = i(:,pdx); jpar = j(:,pdx); kpar = k(:,pdx);
        UI = UserItem; IN = ItemNorm;
        Mijpar = zeros(batch,1); Mikpar = zeros(batch,1);
        for idx = 1:batch
            ii = ipar(idx); jj = jpar(idx); kk = kpar(idx);
            Mijpar(idx) = UI(:,ii)'*UI(:,jj)/(IN(ii)*IN(jj));
            Mikpar(idx) = UI(:,ii)'*UI(:,kk)/(IN(ii)*IN(kk));
            if mod(idx+batch*(pdx-1),million)==0
                w = fprintf([repmat('\b',1,w),'sample: %dM/%dM\n'],(idx+batch*(pdx-1))/million,sample/million) - w; 
            end
        end
        Mij(:,pdx) = Mijpar; Mik(:,pdx) = Mikpar;
    end
end
i = i(:)'; j = j(:)'; k = k(:)'; Mij = Mij(:)'; Mik = Mik(:)';
Yijk = sign(Mij - Mik);
clear UserItem ItemNorm Mij Mik

% Keep only the first m samples that Mij ~= Mik
keep = find(Yijk); keep = keep(1:m);
i = i(keep); j = j(keep); k = k(keep); Yijk = (Yijk(keep)+1)/2;
clear keep

% Split into train and test set
perm = randperm(m);
i_train = i(perm(test_size+1:end)); j_train = j(perm(test_size+1:end));
k_train = k(perm(test_size+1:end)); Yijk_train = Yijk(perm(test_size+1:end));
i_test = i(perm(1:test_size)); j_test = j(perm(1:test_size));
k_test = k(perm(1:test_size)); Yijk_test = Yijk(perm(1:test_size));
clear i j k Yijk perm

% Output sparse data
spdata.train = [i_train(:),j_train(:),k_train(:),Yijk_train(:)];
spdata.test = [i_test(:),j_test(:),k_test(:),Yijk_test(:)];

end