m=50;
n=12;
t = linspace(0,1,m)';

V = vander(t);
V = fliplr(V);
A = V(1:m,1:n);
b = cos(4*t);
%a)
x1 = (A'*A)\(A'*b);


%b)
[Q R] = mgs(A);
v2 = Q'*b;
x2 = R\v2;

%c)
[Q R] = house(A);
v3 = Q'*b;
x3 = R\v3;

%d)
[Q R] = qr(A);
v4 = Q'*b;
x4 = R\v4;

%e)
x5 = A\b;

%f)
[U S V] = svd(A);
v6 = S\(U'*b);
x6 = V*v6;

%g)
%The normal equations exhibit instability.

