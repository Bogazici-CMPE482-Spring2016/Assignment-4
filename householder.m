function [Q,R] = householder(A)
    [m,n] = size(A);
    Q = eye(m);     
    R = A;
    for k=1:n        
        level = m-k+1; % size
        
        e_1 = zeros(level,1);
        e_1(1) = 1;
        
        x = R(k:m,k);
        v_k = sign(x(1))*norm(x)*e_1 + x;
        v_k = v_k / norm(v_k);
        F = eye(level) - 2 * (v_k * v_k');
        Q_k = eye(m,m);
        Q_k(k:end,k:end) = F;
        
        R = Q_k * R;
        Q = Q * Q_k';                
    end 
end