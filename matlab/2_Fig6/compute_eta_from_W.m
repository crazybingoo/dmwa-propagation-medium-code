function eta_val = compute_eta_from_W(W)

if isempty(W)
    eta_val = 0;
    return;
end

W = double(W);
n = size(W, 1);

if size(W, 2) ~= n
    error('Input W must be a square matrix.');
end

if n < 2
    eta_val = 0;
    return;
end

W(~isfinite(W)) = 0;
W(W < 0) = 0;
W(1:n+1:end) = 0;

sum_W = sum(W(:));
if sum_W <= 0
    eta_val = 0;
    return;
end

lambda_1 = max(abs(eig(W)));
if ~isfinite(lambda_1) || lambda_1 <= eps
    eta_val = 0;
else
    eta_val = (sum_W / n) / lambda_1;
end

end
