function [Q,R] = house(A)
[m,n] = size(A);
V = zeros(m,n);
R = A;
for k = 1:n
x = R(k:m,k);
if (x(1) == 0)
s =1;
else
s = sign(x(1));
end
x(1)= s * norm(x) + x(1);
v =x;
v = v/norm(v);
V(k:m, k) = v;
R(k:m, k:n) = R (k:m, k:n) - 2 *v*(v' * R(k:m,k:n));
end
R = R(1:n,1:n);
Q = eye(m,n);
for j = 1:n
for k = n:-1:1
Q(k:m, j) = Q(k:m, j) - 2 * V(k:m, k) * (V(k:m, k)' * Q(k:m,j));
end
end