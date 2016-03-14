function [q, r] = mgs(a)
    [m, n] = size(a);
    q = zeros(m, m); 
    r = zeros(m, n);
    v = a;
    for i = 1:n
        r(i, i) = norm(v(:, i));
        q(:, i) = v(:, i) / r(i, i);
        for j = i+1:n
            r(i, j) = q(:, i)' * v(:, j);
            v(:, j) = v(:, j) - r(i, j) * q(:, i);
        end
    end
end