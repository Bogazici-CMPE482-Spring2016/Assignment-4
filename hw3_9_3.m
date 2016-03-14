A = zeros(15,40);
%a-----------
vertical = [2 2; 2 3; 2 6; 2 7; 3 10; 3 11; 4 18; 4 19; 5 26; 5 27; 6 34; 6 35; 6 38; 6 39];
horizontal = [5 2; 6 2;3 10;4 10; 6 10; 7 10;9 10;10 10;10 18; 11 18;11 26;12 26;6 34; 7 34; 12 34; 13 34];


for i = 1:length(vertical)
    A(vertical(i,1):vertical(i,1)+7,vertical(i,2)) = 1;
end

for i = 1:length(horizontal)
    A(horizontal(i,1),horizontal(i,2):horizontal(i,2)+5) = 1;
end
%spy(A);
%a------------

%b------------
S = svd(A)
figure
plot(S)
title('plot');
figure
semilogy(S)
title('semilogy');
rank(A)
%b------------

%c------------
[U, S, V] = svd(A);
figure
for i = 1:rank(A)
    B = U(:,1:i) * S(1:i, 1:i) * V(:,1:i)'
    pcolor(B)
    colormap(gray)
end
%c------------