function hw3_11_3
    m = 50;
    n = 12;

    t = linspace(0,1,m);

    tmp = fliplr(vander(t));
    A = tmp(:,1:n);
    b = cos(4*t);

    %a-------------------
    x1 = (A'*A) \ (A'*b');
    %--------------------

    %b-------------------
    [Q, R] = mgs(A);
    v = Q' * b';
    x2 = R\v; 
    %--------------------

    %c-------------------

    [W, R] = house(A);
    Q = zeros(m, m);
    for i = 1:m
        e = zeros(m,1);
        e(i) = 1;
        Q(:,i) = formQx(W,e);
    end

    vec = b';
    for i = 1:n
        vec(i:m)=vec(i:m)-2*W(i:m,i)*((W(i:m,i))'*vec(i:m));
    end
    x3 = R\vec;

    %--------------------

    %d-------------------
    [Q, R] = qr(A);
    v = Q' * b';
    x4 = R\v; 
    %--------------------

    %e-------------------
    x5 = A\b';
    %--------------------

    %f-------------------
    [U,S,V] = svd(A);
    x6 = V*(S\(U'*b'));
    %--------------------

    [x1 x2 x3 x4 x5 x6]   
    
   function [Q, R] = mgs(A)
       n = size(A,2);
       V = A;
       Q = zeros(size(A,1), size(A,2));
       R = zeros(n,n);

       for k = 1:n
           R(k,k) = norm(V(:,k));
           Q(:,k) = V(:,k)/R(k,k);
           for j = k+1:n
               R(k,j) = Q(:,k)' * V(:,j);
               V(:,j) = V(:,j) - R(k,j) * Q(:,k);
           end
       end
   end
    function [W,R]=house(A) 
        m=size(A,1);
        n=size(A,2);
        W=zeros(m,n);
        for k=1:n
            x=A(k:m,k);
            e1=eye(m-k+1,1);
            if (x(1)==0) coef=1;
            else coef=sign(x(1));
            end;
            W(k:m,k)=coef*norm(x,2)*e1+x;
            W(k:m,k)=W(k:m,k)/norm(W(k:m,k),2);
            A(k:m,k:n)=A(k:m,k:n)-2*W(k:m,k)*((W(k:m,k))'*A(k:m,k:n));
        end;
        R=triu(A,0);
    end

    function a=formQx(W,vect)
        m=size(W,1);
        n=size(W,2);
        for k=n:-1:1
            vect(k:m)=vect(k:m)-2*W(k:m,k)*((W(k:m,k))'*vect(k:m));
        end;
        a=vect;
    end
end