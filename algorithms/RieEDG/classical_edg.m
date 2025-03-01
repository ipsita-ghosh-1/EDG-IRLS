function P = classical_edg(D,d,distance_matrix)
% input is a euclidean distance matrix
% d is dimension of said matrix
% distance_matrix is 1 if the matrix given is a distance matrix, 0 if its a gram
% matrix then leave it alone

n = size(D,1);
if distance_matrix == 1
    J = eye(n) - (1/n)*ones(n,1)*ones(1,n);
    G = -(1/2)*J*(D)*J;
else
    G = D;
end

[U,E] = eig(G);
E = diag(E);
E = real(E);
[E,idx] = sort(E,'descend');
E = E(1:d);
U = U(:,idx(1:d));
E = diag(E);
P = U*sqrt(abs(E));

return