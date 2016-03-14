%a)
A=zeros(15,40);
%H
A(2:9,2:3)=1;
A(5:6,4:5)=1;
A(2:9,6:7)=1;
%E
A(3:10,10:11)=1;
A(3:4,12:15)=1;
A(6:7,12:15)=1;
A(9:10,12:15)=1;
%L
A(4:11,18:19)=1;
A(10:11,20:23)=1;
%L
A(5:12,26:27)=1;
A(11:12,28:31)=1;
%O
A(6:13,34:35)=1;
A(6:13,38:39)=1;
A(6:7,36:37)=1;
A(12:13,36:37)=1;

%b)
S=svd(A);
plot(S)
title('Plot S');
figure;
semilogy(S);title('Semilogy S');
R=rank(A)

[U S V]=svd(A);
%c
for i=1:rank(A)
   B = U(:,1:i)*S(1:i,1:i)*V(:,1:i)';
   figure
   colormap(gray);
   pcolor(B);
end