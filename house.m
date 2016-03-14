function [q, r] = house(a)
    [m, n] = size(a);
    r = a;
    w = zeros(m, n);
    for k = 1:n
        x = r(k:m, k);
        v = sign(x(1)) * norm(x) * eye(m-k+1, 1) + x;
        v = v ./ norm(v);
        r(k:m, k:n) = r(k:m, k:n) - 2 * v * (v' * r(k:m, k:n));
        w(k:m, k) = v;
    end
    
    q = zeros(m, m);
    for j = 1:m
        e = zeros(m, 1);
        e(j) = 1;
        for k = n:-1:1
            e(k:m) = e(k:m) - 2 * w(k:m, k) * (w(k:m, k)' * e(k:m));
        end
        q(:, j) = e;
    end
end