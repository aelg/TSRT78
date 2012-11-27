whistle = y;
whistle = whistle(2,4300:end);
x = 1/length(whistle):8000/length(whistle):8000;
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

energytot = sum(abs(whistle).^2)/8000
energydom = sum(abs(filtwhistle).^2)/8000


puritytime = 1- energydom./energytot



%%
domlimits=(floor(0.27*length(whistle)/2):floor(0.29*length(whistle)/2));

Energytot = sum(abs(WHISTLE).^2)/(length(WHISTLE)*8000)
Energydom = 2*sum(abs(WHISTLE(domlimits)).^2)/(length(WHISTLE)*8000)

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
abs(roots(Model.a))



%%
%%knä för armodell(2) för vissling ;
error = zeros(1,30);
for order = 1:30;
    [theta, P1, lam, epsilon] = sig2ar(whistle', order);
    error(order) = lam;
end

plot(error);

%%

%psd cehck för vissling
A=psd(spectrum.welch, whistle);
B = psd(spectrum.welch, e);
figure;
plot(A); 
figure;
plot(B)
%%

aaa = y;
aaa = aaa(2,4300:end);

%%
%Knä för aaa 

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
%psd cehck för aaa
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
%Knä för ooo 

error = zeros(1,30);
for order = 1:30;
    [theta, P1, lam, epsilon] = sig2ar(ooo', order);
    error(order) = lam;
end

plot(error);

%%
%psdCheck för ooo
A=psd(spectrum.welch, ooo);
B = psd(spectrum.welch, e);
figure;
plot(A); 
figure;
plot(B)

%%
sentence = y(2,:);


%%
%sentence
sentence = y(2,:);
interval= 160;
order = 10;
intervals=floor(length(sentence)/160);

modsen= [];
for k= 1:intervals
sentence((k-1)*160+1:k*160+1) = detrend(sentence((k-1)*160+1:k*160+1));
m1 = ar(sentence((k-1)*160+1:k*160+1),8);
m1.a;
e=filter(m1.a,1,sentence((k-1)*160+1:k*160+1)'); % m1 <-> AR model of the segment
r=covf(e,100)
r=r(20:end);
[A,D] = max(r);
Roots=roots(m1.a);
for i = 1:length(Roots)
    if (abs(Roots(i)) > 1)
        Roots(i) = 1/Roots(i);
    end
end
m1a= poly(Roots);
input = 1:160;
ehat = sqrt(160*A/D)*(rem(input,D)== 0);
yhat=filter(1,m1a,ehat);
modsen= [modsen yhat];

end

%%
soundsc(modsen,8000)

%%
figure(1);
plot(abs(fft(sentence)));
figure(2);
plot(abs(fft(modsen)));
  


    
  



