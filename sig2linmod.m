<<<<<<< HEAD
function [ th, P, lam, epsi ]= sig2linmod(y,Phi)
%SIG2LINMOD Summary of this function goes here
%   Detailed explanation goes here
th=Phi\y;
R=Phi'*Phi;
epsi=y-Phi*th;
lam = epsi'*epsi/length(y);
P=lam*inv(R);

end

=======
function [th, P, lam, epsi] = sig2linmod(y, Phi);
%Estimate parameters in a linear model
th = Phi\y;
epsi = y-Phi*th;
R = Phi'*Phi;
lam = epsi'*epsi/length(y);
P = lam*inv(R);
>>>>>>> 80433e0f409896159df7ee59761492cc1ca40675
