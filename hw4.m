%%
%%% Constants %%%
epsilon = 1e-05;

%
%%% 9.3.(a) %%%
hello = [ 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
          0 1 1 0 0 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
          0 1 1 0 0 1 1 0 0 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
          0 1 1 0 0 1 1 0 0 1 1 1 1 1 1 0 0 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
          0 1 1 1 1 1 1 0 0 1 1 0 0 0 0 0 0 1 1 0 0 0 0 0 0 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0;
          0 1 1 1 1 1 1 0 0 1 1 1 1 1 1 0 0 1 1 0 0 0 0 0 0 1 1 0 0 0 0 0 0 1 1 1 1 1 1 0;
          0 1 1 0 0 1 1 0 0 1 1 1 1 1 1 0 0 1 1 0 0 0 0 0 0 1 1 0 0 0 0 0 0 1 1 1 1 1 1 0;
          0 1 1 0 0 1 1 0 0 1 1 0 0 0 0 0 0 1 1 0 0 0 0 0 0 1 1 0 0 0 0 0 0 1 1 0 0 1 1 0;
          0 1 1 0 0 1 1 0 0 1 1 1 1 1 1 0 0 1 1 0 0 0 0 0 0 1 1 0 0 0 0 0 0 1 1 0 0 1 1 0;
          0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 0 0 1 1 1 1 1 1 0 0 1 1 0 0 0 0 0 0 1 1 0 0 1 1 0;
          0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 0 0 1 1 1 1 1 1 0 0 1 1 0 0 1 1 0;
          0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 0 0 1 1 1 1 1 1 0;
          0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 0;
          0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
          0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
[m, n] = size(hello);
spy(hello);

%%
%%% 9.3.(b) %%%
[U,S,V] = svd(hello);
plot(U);
plot(S);
plot(V);
semilogy(U);
semilogy(S);
semilogy(V);
r=rank(hello); % 10- The number of columns that are not dependent which is in this case the number of columns that are not equivalent
nzc = sum(diag(S)> epsilon); % 10 the number of non zero values on the digonal of S matrix gives the rank.

%%
%%% 9.3.(c) %%%
B=zeros(m,n,r);
S_ = S;
figure
colormap(gray);
for i=r:-1:1
    B(:,:,i) = U*S_*V';
    S_(i,i) = 0;
    subplot(2,5,i);
    pcolor(B(:,:,i));
    title(strcat('B_{', num2str(i),'}'))
end

%%
%%% 10.4.(a) %%%
% F reflects the plane across the vector whose angle from the positive x-axis is
% (pi-theta)/2
% J rotates the plane by the angle theta clockwise.

%%
%%% 10.4.(b)(c) %%%

% Householder
R = hello;
for k = 1:n
    x = R(k:m,k);
    e = zeros(m-k+1,1); e(1) = 1;
    q = sign(x(1))*norm(x)*e + x
    if norm(q) ~= 0
        q = q./norm(q);
    end
    % This is the core part
    R(k:m, k:n) = R(k:m, k:n) -2*(q)*q'*R(k:m, k:n);
    % The number of the operations in each iteration:
    % (n-k)(m-k) * (2 mult + 1 addition + 1 subtraction) = 4(n-k)(m-k)
    % Therefore the total number of operations is:
    % N = Sum_k^n 4(m-k)(n-k) ~ 4n^2(m-n/3)
    
%     Q(k:m,k) = q;
end

% Givens 
Q = eye(m);
R = hello;
for j = 1:n
    for i = m:-1:(j+1)

        x = R(i-1,j);
        y = R(i,j);
        if y==0
            continue;
        end
        % Get the givens
        c = x/norm([x,y]');
        s = -y/norm([x,y]');
        
        % This is the core part
        R(i-1:i,j:end)= [c -s; s c]*R(i-1:i,j:end); 
        % The number of the operations in each iteration:
        % 2(n-j)* (2 mult + 1 addition) = 6(n-j)
        % The total number of operations for each j:
        % (m-j-1)*6(n-j)
        % Therefore the total number of operations is:
        % N = Sum_j^n (m-j-1)*6(n-j) ~ 6n^2(m-n/3)
        % which is 1.5 times of the total flops of householder algorithm.
               
%         G = eye(m);
%         G([i-1, i],[i-1, i]) = [c -s; s c];
%         Q = Q*G';
    end
end

%%
%%% 11.3.(a)(b)(c)(d)(e)(f)(g) %%%
m=50;
n=12;

t=linspace(0,1,m);
B = fliplr(vander(t));
A = B(1:m,1:n);
b = cos(4*t)';

x = zeros(n,6);

format long;
x(:,1) = linsolve(A,b) % (a)

[Q,R]= mgs(A);
opts.UT = true;
d = Q'*b;
x(:,2) = linsolve(R,d,opts) % (b)

[V, R] = house(A);
d = b;
for k=1:n
    d(k:m,1) = d(k:m,1)-2*V(k:m,k)*V(k:m,k)'*d(k:m,1);
end
x(:,3) = linsolve(R,d,opts) % (c)

[Q,R]= qr(A);
d = Q'*b;
x(:,4) = linsolve(R,d,opts) % (d)

x(:,5) = A\b %(e)

[U,S,V]=svd(A,0);
x(:,6) = V*((U'*b)./diag(S)) % (f)

image(log(x)); % (g)





