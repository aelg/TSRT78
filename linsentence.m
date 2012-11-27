clear all
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

%{
lambada = [];
for i = (1:50)
    [th, P, lam, epsi] = sig2ar(oo',i);
    lambada = [lambada lam];
end
plot(lambada);
%}
%%

na = 20;
cov = [];
est = zeros(length(sen),1);
delay = 0;
Z = [];
alltrain = [];
last_pulse = 20;
for i = 1:160:length(sen)-160
    seg = i:i+159;
    m = ar(detrend(sen(seg)), na);
    e = filter(m.a,1,detrend(sen(seg))'); % m1 <-> AR model of the segment
    r = abs(covf(e,100));
    %plot(r);
    %pause
    cut = 19;
    r = r(cut+1:end);
    cov = [cov; r];
    [ma, ind] = max(r);
    A = ma;
    pulse = ind+cut;
    %pulse = floor(((ind+cut)+last_pulse)/2)
    last_pulse = pulse;
    if max(abs(roots(m.a))) >= 1
        disp(['Unstable: abs = ' num2str(max(abs(roots(m.a))))]);
        r = roots(m.a);
        r2 = [];
        for j = r(:)
            if abs(j) >= 1
                r2 = [r2; 1/j];
            else
                r2 = [r2; j];
            end
        end
        m.a = poly(r2);
    end
    %A = 1;
    %train = sqrt(A*(160/pulse))*(rem((1:160)+delay,pulse) == 1)';
    train = sqrt(A*(rem((1:160)+delay,pulse) == 1))';
    %train = sqrt(A*(rem((1:160)+delay,pulse))) == 1)';
    train = sqrt(A)*(rem((1:160)+delay,pulse) == 1)';
    if A < 0.1
        train = 0.3*sqrt(mean(r))*randn(160, 1)+train;
    end
    alltrain = [alltrain train'];
    %train = train + sqrt(A2*(160/pulse2))*(rem((1:161),pulse2) == 1)';
    delay = rem(160+delay, pulse);
    delay = 0;
    %cc = size(cov)
    %size(est(seg))
    %ll = length([zeros(1, cut) cov zeros(1, 60)])
    %est(seg) = sqrt(1*(160/pulse))*[zeros(1, cut) cov zeros(1, 60)]';
    
    %est(seg)
    %est(seg) = hamming(length(seg)).*est(seg);
    %est(seg) = (rem((1:160),pulse) == 0)';
    [est(seg), Z] = filter(1,m.a,train, Z);
    %disp(Z)
    %for j = seg+na
    %    est(j) = [-m.A*est(j:-1:j-na)];
    %end
end
plot(psd(spectrum.welch, (est)));
%sound(10*est);
[bpb, bpa] = butter(5, [0.12 0.7]);
%est = filtfilt(bpb, bpa, est);
wavwrite(est, 'out.wav');