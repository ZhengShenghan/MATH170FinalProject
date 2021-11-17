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