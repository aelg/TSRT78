load whistle.mat

Fs=8000;
w = y(2, 4300:end);
len = length(w);
x = 0:2/len:1.9999999999;

freq_dom = 0.28;
B = 0.01;

[LPb, LPa] = butter(5, 0.2, 'high');
[BPb,BPa] = butter(5, [freq_dom-B freq_dom+B]);

%w = filtfilt(HPb, HPa, w);
%w = detrend(w, 'constant');
w_filt = filtfilt(BPb, BPa, w);

figure(1);
plot(x,abs(fft(w)));
figure(2);
plot(abs(fft(w_filt)));

E_tot = sig_pow(w)
E_dom = sig_pow(w_filt)

distorsion_time = 1-E_dom/E_tot

w_fft = fft(w);

E_tot = sig_pow(w_fft)
w_fft_dom = [w_fft(floor((freq_dom-B)*len/2):floor((freq_dom+B)*len/2)) w_fft(floor(len-(freq_dom+B)*len/2):floor(len-(freq_dom-B)*len/2))];
E_dom = sum(abs(w_fft_dom).^2)/len;

distorsion_freq = 1- E_dom/E_tot

na = 8
[th, P, lam, epsi] = sig2ar(w',na)
A = 0.1
est = randn(na,1)*A;
for i = na+1:len-na
    if(mod(i, 10) == 0)
        est = [est; est(i-1:-1:i-na)'*-th+randn(1)*A];
    else
        est = [est; est(i-1:-1:i-na)'*-th];
    end
end
plot(abs(fft(detrend(filtfilt(BPb,BPa,est), 'constant'))));
plot(est);
sound(est*0.1);
sound(w);
