function cliques = maximalCliques_bruteforce(A)

n = size(A,1);
cliques = {};

if n == 0
    return;
end

% Bron-Kerbosch 算法
R = [];
P = 1:n;
X = [];
cliques = bron_kerbosch(R, P, X, A);

end

function cliques = bron_kerbosch(R, P, X, A)

cliques = {};

if isempty(P) && isempty(X)
    cliques = {R};
    return;
end

P_local = P;

for v = P_local
    Nv = find(A(v,:) > 0);

    R_new = [R, v];
    P_new = intersect(P, Nv);
    X_new = intersect(X, Nv);

    sub_cliques = bron_kerbosch(R_new, P_new, X_new, A);
    cliques = [cliques, sub_cliques]; %#ok<AGROW>

    P(P == v) = [];
    X = [X, v];
end

end