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
title('DFT a-sound');
subplot(2, 1, 2);
plot(x, abs(fft(oo)));
title('DFT o-sound');
xlabel('Hz');
ylabel('Amplitude');
%spectrum(oo);
%plot(psd(spectrum.welch, (oo)));
%%

% Validation of different orders.
varaa = [];
varoo = [];
resaa = [];
resoo = [];
for na = (1:20)
    m = ar(aa_est, na);
    residual = filter(m.a, 1, aa_val);
    varaa = [varaa residual*residual'/(aa_val*aa_val')];
    resaa = [resaa residual'];
    
    m = ar(oo_est, na);
    residual = filter(m.a, 1, oo_val);
    varoo = [varoo residual*residual'/(oo_val*oo_val')];
    resoo = [resoo residual'];
end

% Plot loss function for different model orders.
figure(2);
subplot(2,1,1);
plot(varaa);
title('Loss function of validation data a-sound');
ylabel('Variance');
ylim([0 0.08]);
xlabel('Model Order');
subplot(2,1,2);
plot(varoo);
[mi, i] = min(varoo)
title('Loss function of validation data o-sound');
ylabel('Variance');
ylim([0 0.08]);
xlabel('Model Order');
%%

% Calculate whiteness for whiteness test for different model orders.
aa_pr_s = [];
oo_pr_s = [];
for k = 1:20
    aa_pr_s = [aa_pr_s sum(abs(diff(sign(resaa(:,k)))))/2/len_val];
    oo_pr_s = [oo_pr_s sum(abs(diff(sign(resoo(:,k)))))/2/len_val];
end
aa_pr_s
oo_pr_s

%%
% Calculate residual covariance for different model orders.
start = 10;
for k = (start:start+3)
    figure(3);
    subplot(4,1,k-start+1);
    [Reyaa, kk] = sig2crosscovfun(resaa(:,k), aa_val', 10);
    plot(kk, Reyaa/sqrt(var(resaa(:,k))*var(aa_val')));
    aa_max = max(xcorr(resaa(:,k), aa_val'))
    title(sprintf('Residual Covariance plot a-sound order %d', k));
    ylabel('Correlation');
    xlabel('k');
    xlim([1 10]);
    if k ~= 1; ylim([-0.2 0.2]); end;
    %plot(resaa(:,k));
    figure(4);
    subplot(4,1,k-start+1);
    [Reyoo, kk] = sig2crosscovfun(resoo(:,k), oo_val', 10);
    plot(kk, Reyoo/sqrt(var(resoo(:,k))*var(oo_val')));
    title(sprintf('Residual Covariance plot o-sound order %d', k));
    ylabel('Correlation');
    xlabel('k');
    xlim([1 10]);
    if k ~= 1; ylim([-0.2 0.2]); end;
end
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
%e = e+rand(len,1)*A*1;
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