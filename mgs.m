function [Q, R] = mgs(A)
[m,n] = size(A);
R = zeros(m,n);
Q = zeros(m);
for i=1:min(m,n)
    R(i,i) = norm(A(:,i));
    if R(i,i) ~= 0
        Q(:,i) = A(:,i)/R(i,i);
    end
    for j=i+1:n
        R(i,j) = Q(:,i)'*A(:,j);
        A(:,j) = A(:,j) - R(i,j)*Q(:,i);
    end
end