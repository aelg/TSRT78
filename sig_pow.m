function res = sig_pow(x, Ts)

%sum(abs(x).^2)
%length(x)
s = size(x);
if(s(1) == 1)
    res = Ts*x*x';
elseif s(2) == 1
    res = Ts*x'*x;
else
    res = 0;
end

end