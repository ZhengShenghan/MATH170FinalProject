function [obj,x,loop,status]=one_norm(A,b)
    %m>n
        [m,n]    = size(A);
        loop = 0;
        obj = -1;
        %msg = ''
        %construct c
        C_1 = zeros(1,2*n);
        C_2 = ones(1,2*m);
        c = [C_1 C_2];
        c = transpose(c);
        %constrcut new A
        I = eye(m);
        A1 = [A -A -I I];
        [x,status,loop] = LPSol(A1,b,c,loop);
        obj = transpose(c)*x;
        
end