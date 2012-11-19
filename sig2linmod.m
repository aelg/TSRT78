function [th, P, lam, epsi] = sig2linmod(y, Phi);
%Estimate parameters in a linear model
th = Phi\y;
epsi = y-Phi*th;
R = Phi'*Phi;
lam = epsi'*epsi/length(y);
P = lam*inv(R);

