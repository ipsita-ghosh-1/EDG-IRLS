function [U,sing,V,singmin] = update_matrixdec_from_handle(X_handle,...
    r,N_SVD,decomposition,opts,varargin)
%update_matrixdec_from_handle Given a cell 'X_handle' containing the
%information about a matrix, computes the 'r' first eigenvector/ singular
%vector pairs and corresponding eigen/ singular values, where
%'decomposition' indicates with matrix decomposition exactly to use.
% =========================================================================
% Parameters
% ----------
% X_handle{1}: Matrix handle for matrix-vector multiplication X*v
% X_handle{2}: Matrix handle for matrix-vector multiplication X'*w
% X_handle{3}: Number of rows of X
% X_handle{4}: Number of columns of X
% X_handle{5}: If not empty: Full matrix X.
%
% r: int. Rank parameter.
%
% decomposition: character string.
%
% Returns
% ----------
Uinit = [];
Vinit = [];
symmetricflag = false;

if ~isempty(varargin)
    for tt = 1:2:length(varargin)
        switch varargin{tt}
            case 'Uinit' 
                Uinit = varargin{tt+1};
            case 'Vinit'
                Vinit = varargin{tt+1};
                
        end
    end
end
if isfield(opts,'symmetricflag')
   symmetricflag = opts.symmetricflag; 
end

OPTS.isreal = true; % for primme algorithms
OPTS.tol = opts.tol_CG;

if (strcmp(decomposition,'svdexact') && length(X_handle) < 5)
   error("Please provide the full matrix in the fifth cell of the matrix handle of X, since you chose an exact decomposition.")
end
if strcmp(decomposition,'svd')
        [U, sing, V]=bksvd_mod(X_handle(1:4),...
            r, N_SVD);
        if symmetricflag
           % signdifferences = sign(sum(sign(U)./sign(V),1));
           % V = U.*signdifferences;
           V = U;
        end
        %V = U;
        singmin = [];
elseif strcmp(decomposition,'svd_primme')
    X_handle_primme = @(v,flag) get_primme_handle(v,flag,X_handle);
%     if ~isempty(Uinit)
%         OPTS.u0 = Uinit;
%     end
%     if ~isempty(Vinit)
%         OPTS.v0 = Vinit;
%     end
    [U, sing, V]= primme_svds(X_handle_primme,X_handle{3},X_handle{4},r,'L',OPTS);
    singmin = [];
elseif strcmp(decomposition,'svdexact')
    [U, sing, V] = svds(X_handle{5},r);
%         V = U;
    singmin = [];
elseif strcmp(decomposition,'eig_abs')
    d = min(X_handle{3},X_handle{4}); 
%     if ~isempty(Uinit)
%         OPTS.v0 = Uinit;
%     end
    TARGET = 'LM';
    [U, sing]= primme_eigs(X_handle{1},d,r,TARGET,OPTS);
    V = U*sign(sing);
    sing = abs(sing);
    singmin = [];
elseif strcmp(decomposition,'eig')
    d = min(X_handle{3},X_handle{4}); 
%     if ~isempty(Uinit)
%         OPTS.v0 = Uinit;
%     end
    TARGET = 'LA';
    [U, sing]= primme_eigs(X_handle{1},d,r,TARGET,OPTS);
    V = U*sign(sing);
    sing = abs(sing);
    singmin = [];
elseif strcmp(decomposition,'eigexact_abs')
    d = min(X_handle{3},X_handle{4}); 
    [U, sing] = eigs(X_handle{1},d,r,'largestabs');
    V = U*sign(sing);
    sing = abs(sing);
    singmin = [];
elseif strcmp(decomposition,'eigexact')
    d = min(X_handle{3},X_handle{4}); 
    [U, sing] = eigs(X_handle{1},d,r,'largestreal');
    V = U;
%         [Umin, singmin] = eigs(X_c_handle{5},2,'smallestreal');
% %         sing=diag(sing);
%         singmin=diag(singmin);
%         sing=[sing;singmin];
%         sing=diag(sing);
    singmin = [];
else 
    error('Please indiciate whether svd or eig is used in the weight updates.')
end
end

function  Xtimesvec = get_primme_handle(v,flag,X_handle)
if strcmp(flag,'notransp')
    Xtimesvec = X_handle{1}(v);
elseif strcmp(flag,'transp')
    Xtimesvec = X_handle{2}(v);
end
end