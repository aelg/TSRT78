%linwovel.m

% Load and detrend data
Fs=8000;
load aaa.mat
aa = detrend(y(2, 6000:end));
load ooo.mat
oo = detrend(y(2, 6000:end));
len = length(aa);
x = 0:Fs/len:Fs-1e-4;

% Divide into estimation and validation data.
aa_est = aa(1:end*2/3);
aa_val = aa(end*2/3+1:end);
oo_est = oo(1:end*2/3);
oo_val = oo(end*2/3+1:end);
len_est = length(aa_est);
len_val = length(aa_val);
x_est = 0:2/len_est:Fs-1e-9;
x_val = 0:2/len_val:Fs-1e-9;

figure(1);
subplot(2, 1, 1);
length(x)
length(abs(fft(aa)))
plot(x, abs(fft(aa)));
xlabel('Hz');
ylabel('Amplitude');
legend('DFT a-sound');
subplot(2, 1, 2);
plot(x, abs(fft(oo)));
legend('DFT o-sound');
xlabel('Hz');
ylabel('Amplitude');
%spectrum(oo);
%plot(psd(spectrum.welch, (oo)));
%%

% Validation of different orders.
varaa = [];
varoo = [];
for na = (1:20)
    m = ar(aa_est, na);
    residual = filter([m.a(1:end)], 1, aa_val);
    varaa = [varaa residual*residual'/length(residual)];
    
    m = ar(oo_est, na);
    residual = filter([m.a(1:end)], 1, oo_val);
    varoo = [varoo residual*residual'/length(residual)];
    %subplot(2, 1, 1);
    %plot(aa_val)
    %subplot(2, 1, 2);
    %plot(residual)
end
figure(2);
subplot(2,1,1);
plot(varaa);
title('Residual Variance a-sound');
ylabel('Variance');
xlabel('Model Order');
subplot(2,1,2);
plot(varoo);
title('Residual Variance o-sound');
ylabel('Variance');
xlabel('Model Order');

%%

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