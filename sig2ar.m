<<<<<<< HEAD
function [ th, P, lam, epsi ] = sig2ar( y, na )
%SIG2AR Summary of this function goes here
%   Detailed explanation goes here

N=length(y);
ind=(na:N-1)'*ones(1,na)-ones(N-na,1)*(0:na-1);
Phi=-reshape(y(ind),size(ind));
Y=y(na+1:N);

[ th, P, lam, epsi ]= sig2linmod(Y,Phi);
if nargout==0
    disp([' a = ', num2str([th'])])
    disp(['std = ', num2str(diag(P)')])
    disp(['lam = ', num2str(lam)])

end

=======
function  [th, P, lam, epsi] = sig2ar(y, na);

%Estimate paramaters in an auto-regressive model
N=length(y);
ind = (na:N-1)'*ones(1,na) - ones(N-na,1)*(0:na-1);
Phi = -reshape(y(ind),size(ind));
Y=y(na+1:N);
[th, P, lam, epsi] = sig2linmod(Y,Phi);
if nargout == 0
    disp([' a = ', num2str([th'])])
    disp(['std = ', num2str(diag(P)')])
    disp(['lam = ', num2str(lam)])
end

end
>>>>>>> 80433e0f409896159df7ee59761492cc1ca40675
