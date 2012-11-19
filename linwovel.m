%linwovel.m

Fs=8000;
load aaa.mat
aa = detrend(y(2, 6000:end));
load ooo.mat
oo = detrend(y(2, 6000:end));
len = length(aa);
x = 0:2/len:1.9999999999;

figure(1);
plot(x, abs(fft(oo)));
figure(2);
%plot(x, abs(spectrum(oo)));
%spectrum(oo);
plot(psd(spectrum.welch, (oo)));

na = 20

lambada = [];
for i = (1:50)
    [th, P, lam, epsi] = sig2ar(oo',i);
    lambada = [lambada lam];
end
figure(4);
plot(lambada);


[th, P, lam, epsi] = sig2ar(oo',na);
%len = len*2;
W0 = 0.035
%W0 = 0.016
root_angle = angle(roots(th))
for i = 1:length(root_angle)
    if root_angle(i) == 0
        root_angle(i) = 100000;
    end
end
pulse = floor(2/W0)
A = 0.001
%pulse = floor(8/(min(abs(root_angle))/pi))
est = rand(na,1);
%e = randn(len,1).*(rem((1:len),pulse)<pulse/2)';
%e = randn(len,1).*(rem((1:len),pulse) == 0)';
e = A*sqrt(len/pulse)*(rem((1:len),pulse) == 0)';
e = e+rand(len,1)*A*1;
%e = e .* randn(size(e));
%e = e + 0.005*randn(size(e));
%e = conv(e, triang(pulse/2), 'same');
for i = na+1:len
    est = [est; est(i-1:-1:i-na)'*th+e(i)];
end
th
est = filter(1,[1 th'], e);


figure(3);
%plot(x,abs(fft(est)));
plot(psd(spectrum.welch, (est)));
plot(x, abs(fft(est)));
est = est(200:end-200);
%sound(est)
%sound(oo)