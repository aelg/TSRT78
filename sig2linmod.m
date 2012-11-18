function [ th, P, lam, epsi ]= sig2linmod(y,Phi)
%SIG2LINMOD Summary of this function goes here
%   Detailed explanation goes here
th=Phi\y;
R=Phi'*Phi;
epsi=y-Phi*th;
lam = epsi'*epsi/length(y);
P=lam*inv(R);

end

