Fs=8000;
load sentence.mat
sen = detrend(y(2,:)*100);
len = length(sen);
x = 0:2/len:1.9999999999;
%[b,a] = butter(5, [0.2 0.8]);
%sen = filtfilt(b,a,sen);
%sound(sen)
%%
%figure(1);
%plot(x, abs(fft(sen)));
%figure(2);
%plot(x, abs(spectrum(oo)));
%spectrum(oo);
%plot(psd(spectrum.welch, (sen)));

na = 8
%{
lambada = [];
for i = (1:50)
    [th, P, lam, epsi] = sig2ar(oo',i);
    lambada = [lambada lam];
end
plot(lambada);
%}
%%
cov = [];
est = zeros(length(sen),1);
for i = 1:160:length(sen)-160
    seg = i:i+160-1;
    m = ar(detrend(sen(seg)), na);
    e = filter(m.a,1,detrend(sen(seg))'); % m1 <-> AR model of the segment
    r = abs(covf(e,100));
    cut = 19;
    r = r(cut+1:end);
    cov = [cov; r];
    [ma, ind] = max(r);
    %est = [est; rand(na,1)];
    A = ma;
    pulse = (ind+cut);
    disp(pulse)
    %pulse = floor(pulse/2);
    if max(abs(roots(m.a))) >= 1
        disp(['Unstable: abs = ' num2str(max(abs(roots(m.a))))]);
        r = roots(m.a);
        r2 = [];
        for j = r(:)
            if abs(j) > 1
                r2 = [r2; 1/j];
            else
                r = [r2; j];
            end
        end
        m.a = poly(r2);
    end
    train = sqrt(A*(160/pulse))*(rem((1:160),pulse) == 0)';
    %cc = size(cov)
    %size(est(seg))
    %ll = length([zeros(1, cut) cov zeros(1, 60)])
    %est(seg) = sqrt(1*(160/pulse))*[zeros(1, cut) cov zeros(1, 60)]';
    
    %est(seg)
    %est(seg) = hamming(length(seg)).*est(seg);
    %est(seg) = (rem((1:160),pulse) == 0)';
    est(seg) = filter(1,m.a,train);
    %for j = seg+na
    %    est(j) = [-m.A*est(j:-1:j-na)];
    %end
end
plot(psd(spectrum.welch, (est)));
sound(10*est);
wavwrite(10*est, 'out.wav');
%sound(sen);
%%
%{
[th, P, lam, epsi] = sig2ar(sen',na);
%len = len*2;
W0 = 0.035
W0 = 0.016
root_angle = angle(roots(th))
for i = 1:length(root_angle)
    if root_angle(i) == 0
        root_angle(i) = 100000;
    end
end
pulse = floor(2/W0)
A = 0.1
%pulse = floor(8/(min(abs(root_angle))/pi))
est = rand(na,1);
%e = randn(len,1).*(rem((1:len),pulse)<pulse/2)';
%e = randn(len,1).*(rem((1:len),pulse) == 0)';
e = A*sqrt(len/pulse)*(rem((1:len),pulse) == 0)';
%e = e+rand(len,1)*0.01
%e = e .* randn(size(e));
%e = e + 0.005*randn(size(e));
%e = conv(e, triang(pulse/2), 'same');
for i = na+1:len
    est = [est; est(i-1:-1:i-na)'*-th+e(i)];
end

figure(3);
%plot(x,abs(fft(est)));
plot(psd(spectrum.welch, (est)));
est = est(200:end-200);
sound(est)
sound(sen)
%}