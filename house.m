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