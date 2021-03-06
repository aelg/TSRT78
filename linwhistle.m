load whistle.mat

Fs=8000;
w = detrend(y(2, 4300:end));
len = length(w);
x = 0:2/len:1.9999999999;
a = 0.27:2/len:0.29

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

Ts = 1/8000;

E_tot = sig_pow(w, Ts)
E_dom = sig_pow(w_filt, Ts)

distorsion_time = 1-E_dom/E_tot

w_fft = fft(w);

%E_tot = sig_pow(w_fft, Ts)
E_tot = w_fft*w_fft' * Ts/len
w_fft_dom = [w_fft(floor((freq_dom-B)*len/2):floor((freq_dom+B)*len/2)) w_fft(floor(len-(freq_dom+B)*len/2):floor(len-(freq_dom-B)*len/2))];
E_dom = w_fft_dom*w_fft_dom'*Ts/len;

distorsion_freq = 1- E_dom/E_tot


%{
lambada = [];
for i = (1:50)
    [th, P, lam, epsi] = sig2ar(w',i);
    lambada = [lambada lam];
end
figure(4);
plot(lambada);
%}


na = 2
[th, P, lam, epsi] = sig2ar(w',na);
rr = roots([1 th'])
root_angle = angle(roots([1 th']))*2./pi
root_abs = abs(rr)
A = 0.1;
est = rand(na,1)*A;
for i = na+1:len-na
    if(mod(i, 1) == 0)
        est = [est; est(i-1:-1:i-na)'*-th+randn(1)*A];
    else
        est = [est; est(i-1:-1:i-na)'*-th];
    end
end
figure(5);
plot(abs(fft(detrend(filtfilt(BPb,BPa,est), 'constant'))));
plot(est);
%sound(est*0.1);
%sound(w);
figure(3);
plot(abs(fft(filter(1,th,[zeros(200, 1); 1; zeros(200,1)]))));
[b,a] = arma2spec(1, th', 1);
est2 = filter(b,a,randn(10000,1).*(rem((1:10000),10)==0)'*0.1);
figure(6);
plot(psd(spectrum.welch, (w)));
figure(7);
plot(psd(spectrum.welch, (est)));
%plot(abs(fft(est2)));