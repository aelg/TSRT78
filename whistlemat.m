whistle = y;
whistle = whistle(2,4300:end);
x = 2/length(whistle):2/length(whistle):2;
WHISTLE = zeros(1,length(whistle));
WHISTLE= fft(whistle);
%%
soundsc(whistle,8000)
%%
figure;
plot(x,whistle)
figure;
plot(x,abs(fft(whistle)));

%%
wn= [0.27 0.29];
[B A] = butter(5, wn);
filtwhistle = filtfilt(B, A, whistle);

energytot = sum(abs(whistle).^2)/8000;
energydom = sum(abs(filtwhistle).^2)/8000;


puritytime = 1- energydom./energytot



%%
domlimits=(floor(0.27*length(whistle)/2):floor(0.29*length(whistle)/2));

Energytot = sum(abs(WHISTLE).^2)/8000;
Energydom = 2*sum(abs(WHISTLE(domlimits)).^2)/8000;

purityfreq = 1 - Energydom./Energytot

%%
order = 2;
[theta, P1, lam, epsilon] = sig2ar(whistle', order);

e = zeros(1,length(whistle));
e(1) = rand(1,1);
e(2)= rand(1,1);
for k = 3:length(whistle)
if (rem(k, 5)== 0) 
    e(k) = randn(1,1);
end
    for i=1:order
        if (k-i > 0)
            e(k) = e(k) - theta(i)*e(k-i);
        end
    end
end

e= e/max(e);
theta

%%    
soundsc(e,8000,[-1 1])

lam

%%
%armodell enligt matlab
Model = ar(whistle, 2);





%%
%%kn� f�r armodell(2) f�r vissling ;
error = zeros(1,30);
for order = 1:30;
    [theta, P1, lam, epsilon] = sig2ar(whistle', order);
    error(order) = lam;
end

plot(error);

%%

aaa = y;
aaa = aaa(2,4300:end);

%%
%Kn� f�r aaa 

error = zeros(1,30);
for order = 1:30;
    [theta, P1, lam, epsilon] = sig2ar(aaa', order);
    error(order) = lam;
end

plot(error);

%%
order = 8;
[theta, P1, lam, epsilon] = sig2ar(aaa', order);

e = zeros(1,length(whistle));
e(1) = rand(1,1);
e(2)= rand(1,1);
for k = 3:length(whistle)
if (rem(k, 5)== 0) 
    e(k) = randn(1,1);
end
    for i=1:order
        if (k-i > 0)
            e(k) = e(k) - theta(i)*e(k-i);
        end
    end
end

e= e/max(e);


%%
%psd cehck f�r aaa
A=psd(spectrum.welch, aaa);
B = psd(spectrum.welch, e);
figure;
plot(A); 
figure;
plot(B)

%%

ooo = y;
ooo = ooo(2,4300:end);
%%
order = 3;
[theta, P1, lam, epsilon] = sig2ar(ooo', order);

e = zeros(1,length(whistle));
e(1) = rand(1,1);
e(2)= rand(1,1);
for k = 3:length(whistle)
if (rem(k, 5)== 0) 
    e(k) = randn(1,1);
end
    for i=1:order
        if (k-i > 0)
            e(k) = e(k) - theta(i)*e(k-i);
        end
    end
end

e= e/max(e);

%%
%Kn� f�r ooo 

error = zeros(1,30);
for order = 1:30;
    [theta, P1, lam, epsilon] = sig2ar(ooo', order);
    error(order) = lam;
end

plot(error);

%%
%psdCheck f�r ooo
A=psd(spectrum.welch, ooo);
B = psd(spectrum.welch, e);
figure;
plot(A); 
figure;
plot(B)



