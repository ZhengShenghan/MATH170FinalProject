function [x,y,B,flg,loop] = Simplex(A,b,c,B,loop)
%
% this function solves the Phase II LP 
% with the simplex method
%
% first set up the simplex tableau
%
format short g;
f       = ['Starting Phase II Simplex Iteration... '];
flg = nan;
%
%
[T,x,y] = GetT(A,b,c,B);
obj     = transpose(c)*x;
[m,n]   = size(A);
disp(['Initial Objective = ', num2str(obj)]);
%
ITER    = 0;
simplex = 1;
%
while (simplex == 1)
%
% determine the next s and r values.
%
    loop = loop + 1;
    y        = transpose(T(end,2:end));
    [zmin,s] = min(c-transpose(A)*y); 
%
% check for convergence.
%
   if (abs(zmin) < 1e-14)
       disp('Simplex Method has converged');
       simplex = 0;
%       disp('Displaying Optimal Basis');
%       disp(transpose(B));
       x    = zeros(n,1);
       x(B) = T(1:end-1,1);
       obj  = transpose(c)*x;
       disp(['Optimal Objective = ', num2str(obj),' after ', num2str(ITER), ' iterations']);
%       disp('Displaying Optimal solution x, c-A^T*y and their componentwise product');
%       disp([x c-transpose(A)*y x.*(c-transpose(A)*y)]);
       continue;
   end

   t        = T(1:end-1,2:end)*A(:,s);
   [flg,r]  = Revisedgetr(n,s,B,T,t);
   if (flg == 1)
       disp('LP is degenerate');
       simplex = 0;
       continue;
   end
   if (r < 1)
       disp('LP has no lower bound');
       simplex = 0;
       flg     = - 1;
       continue;
   end
   x    = zeros(n,1);
   x(B) = T(1:end-1,1);
   ITER = ITER + 1;
%   f    = ['Iteration ', num2str(ITER), ' Obj ', num2str(transpose(c)*x), '. Smallest component in c-A^T*y: ', ... 
%         num2str(zmin), ' at s =', num2str(s), '. Component r = ', num2str(r), ' out of basis'];
%
%   disp('Current solution index set')
%   disp(transpose(B))
%   disp(f);
%
% update the revised simplex tableau.
%
   [T,B1,flg]=RevisedSimplexTableau(B,r,s,t,zmin,T);      
   if (flg == 1)
       disp('LP is degenerate');
       simplex = 0;
       continue;
   end
   B   = B1;
end
