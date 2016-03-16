function [V,R] = house(A)
[m, n] = size(A);
R=A;
V=zeros(m);
for k=1:n
    x = R(k:m,k);
    e = zeros(m-k+1,1); 
    e(1) = 1;
    v = sign(x(1))*norm(x)*e + x;
    if norm(v) ~= 0
        v = v./norm(v);
    end
    R(k:m, k:n) = R(k:m, k:n) -2*(v)*v'*R(k:m, k:n);    
    V(k:m,k) = v;
end