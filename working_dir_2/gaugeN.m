function [X,Y] = gaugeN(Un,Vn)

%
%gauge transformation of the normalization part in the pencil
%
% the solution is such that: X*Un*X=Id, Y.'*Vn*Y.'=Id

X = sqrtm(pinv(Un));
Y = sqrtm(pinv(Vn));
