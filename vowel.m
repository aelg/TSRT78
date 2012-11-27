%vowel.m

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
start = 1;
for k = (start:start+3)
    figure(3);
    subplot(4,1,k-start+1);
    % Check course book for code for this function
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
% Generate soundfiles.

maa = ar(aa_est, 8);
moo = ar(oo_est, 4);
pulseaa = floor(8000/130);
pulseoo = floor(8000/130);

A = 0.1;
eaa = A*(rem((1:len),pulseaa) == 0)';
eoo = A*(rem((1:len),pulseoo) == 0)';

estaa = filter(1,maa.a, eaa);
estoo = filter(1,moo.a, eoo);
wavwrite(estaa, 'a.wav');
wavwrite(estoo, 'o.wav');
