function [T,x,y] = GetT(A,b,c,B)
%
% this sets up tablau
%
% 

[m,n] = size(A);
x    = zeros(n,1);
x(B) = A(:,B)\ b;
%x(B) = A(B,:)\ b;

T    = A(:,B)\[b eye(m)];
%T    = A(B,:)\[b eye(m)];
y    = T(:,2:end)'*c(B);
T    = [T;[c'*x,y']];