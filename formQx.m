function a=formQx(W,vect)
    m=size(W,1);
    n=size(W,2);
    for k=n:-1:1
        vect(k:m)=vect(k:m)-2*W(k:m,k)*((W(k:m,k))'*vect(k:m));
    end;
    a=vect;
end