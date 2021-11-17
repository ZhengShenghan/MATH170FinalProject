function [obj,x,loop,status]=test(A,b)
    [m,n]    = size(A);
    loop = 0;
        %obj = -1;
        %msg = ''
        %construct c
        C_1 = zeros(1,2*n);
        C_2 = ones(1,2*m);
        c = [C_1 C_2];
        c = transpose(c);
        %constrcut new A
        I = eye(m);
        A = [A -A -I I];
        [x,status] = LPSol(A,b,c);
        obj = transpose(c)*x;
end


function [T,x,y] = GetT(A,b,c,B)
%
% this sets up tablau
%
% 

[m,n] = size(A);
x    = zeros(n,1);
x(B) = A(:,B)\ b;

T    = A(:,B)\[b eye(m)];
y    = T(:,2:end)'*c(B);
T    = [T;[c'*x,y']];
end
function [x,status] = LPSol(A,b,c)

% code minimizes c^T x 
  % subject to: x >=0
%               A*x = b
%
% this function solves an LP
%
[m,n]    = size(A);
msk      = find(b < 0);
b(msk)   = - b(msk);
A(msk,:) = - A(msk,:);
msk      = find(b == 0);
if (length(msk)>0)
    disp(['LP is degenerate']);
    status.feas = 'degenerate';
    x           = NaN*ones(n,1);
    return;
end
% 
% Phase I
%
Ahat     = [A, eye(m)];
chat     = zeros(m+n,1);
chat(n+1:n+m)  = ones(m,1);
Bhat           = transpose((n+1:n+m));
[xhat,yhat,B,simplex]  = Simplex(Ahat,b,chat,Bhat);
if (simplex == 1)
    disp(['LP is degenerate']);
    status.feas = 'degenerate';
    x           = NaN*ones(n,1);
    return;
end
disp('simplex method successful for Phase I')
if (max(B)>n)
    disp(['LP is infeasible, Phase I objective = ', num2str(transpose(chat)*xhat)])
    status.feas = 'infeasible';
    x = [];
    return;
end
% 
% Phase II
%
[x,~,B,simplex]  = Simplex(A,b,c,B);
if (simplex == 1)  
    disp(['LP is feasible but degenerate']);
    status.feas = 'degenerate';
    x           = NaN*ones(n,1);
    return;
end
if (simplex == -1)
    disp(['LP is feasible but without optimal']);
    status.feas = 'unbounded';
    status.B    = B;
    return;
end
disp('simplex method successful for Phase II')
disp(['LP optimal reached at objective = ', num2str(transpose(c)*x)])
status.feas = 'optimal';
status.B    = B;
end
function [Tnew,Bnew,flg] = RevisedSimplexTableau(B,r,s,t,zmin,T)
%
% This function updates a RevisedSimplexTableau
% 
%
% On input: 
% B: current basis index
% T: current RevisedTableau
% r: B(r)-th column is to LEAVE the index.
% t: Pivot column
% s: s-th  column is to JOIN the index.
% zmin: s-th component in c-A'*y
%

% On output: 
%
% flg:  flg == 0 indicates SUCCESS in updating,
%       flg == 1 indicates FAILURE in updating,
% Bnew: New basis index
% Tnew: New Tableau
%

%
% initialize flg.
%
flg     = 0;
%
% find dimensions of T.
%
[mt,nt] = size(T); 
%
% Set up Bnew
%
B     = B(:);
Bnew  = [B(1:r-1);s;B(r+1:mt-1)];
%
% Setup Tnew
%
Tnew         = zeros(mt,nt); 
if (t(r) == 0)
%
% This is indication of degeneracy. Quit.
%
    flg = 1;
    return;
end
%
% This is the normal case. Proceed. 
%
Temp              = T(r,:)/t(r);
Tnew(1:mt-1,:)    = T(1:mt-1,:) - t*Temp;
Tnew(r,:)         = Temp;
theta0            = zmin/t(r);
Tnew(mt,:)        = T(mt,:) + theta0*T(r,:);
return
end
function [flg, r] = Revisedgetr(n,s,B,T,t)
%
% find the index to kick out
%
% On input: 
% B: current basis index
% T: current Revised Tableau
% t: current pivot column
% s: s-th  column is to JOIN the index.
% n: number of unknowns.
%

% On output: 
%
% r: B(r)-th column is to LEAVE the index.
%    r < 0 indicates unbounded LP.
%

flg = 0;
x   = zeros(n,1);
x(B)= T(1:end-1,1);
if (max(t)<n*eps)
    r = -1;
    return;
end
mask        = find(t>0);
[lambda, r] = min(x(B(mask))./t(mask));
r           = mask(r);
if (lambda < 1e-14)
    flg  = 1;
end
return
end
function [x,y,B,flg] = Simplex(A,b,c,B)
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
end


test_A = []
