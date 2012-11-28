%load the signal
load('whistle.mat');
whistle = y;
whistle = whistle(2,4300:end);
x = 1/length(whistle):8000/length(whistle):8000;
WHISTLE= fft(whistle);
%%
%sound of the recorded signal 
soundsc(whistle,8000)
%%
%plots of the recorded signal and dft of the recorded signal.
figure;
plot(x,whistle)
figure;
plot(x,abs(fft(whistle)));

%%
%Energy and purity calculations in the timedomain 
wn= [0.27 0.29];
[B A] = butter(5, wn);
filtwhistle = filtfilt(B, A, whistle);

energytot = sum(abs(whistle).^2)/8000
energydom = sum(abs(filtwhistle).^2)/8000

puritytime = 1- energydom./energytot

%%
%Energy and purity calculations in the fourier domain.
domlimits=(floor(0.27*length(whistle)/2):floor(0.29*length(whistle)/2));

Energytot = sum(abs(WHISTLE).^2)/(length(WHISTLE)*8000)
Energydom = 2*sum(abs(WHISTLE(domlimits)).^2)/(length(WHISTLE)*8000)

purityfreq = 1 - Energydom./Energytot

%%
% making an AR-model of order 2, resample it and play it
order = 2;
[Model refl] = AR(whistle,order, 'burg','now' ) ;

e = randn(1,length(whistle));
yhat = filter(1,Model.a,e);
soundsc(yhat,8000)

%%
%counting absolute value of the poles of the AR-filter.
abs(roots(Model.a))

%%
%%knee for AR-models of whistle ;
error = zeros(1,30);
e = randn(1,length(whistle));
whistle2 = whistle/max(whistle);
for order = 1:30;
[Model refl] = AR(whistle,order, 'burg','now' ) 
error(order) = refl(2,order+1)
end

plot(1:30,error);

%%
%psd check for whistle
A=psd(spectrum.welch, whistle);
B = psd(spectrum.welch, yhat);
figure;
plot(A); 
figure;
plot(B)


