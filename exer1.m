%{
N=32;
omega = 1;
[X, w] = dtft(cos(omega*(0:1:N-1)));
plot(w-pi, abs(X))

y= []
for i = 1:3:length(u2)
    y = [y; u2(i)]
end

T = 3.9e-4
figure(1);
plot(abs(dtft(u2, T)), 'r');
hold on;
plot(abs(dtft(y, T)));
hold off;
figure(2);
plot(u2, 'r');
hold on;
plot((1:3:length(u2)), y);
hold off;
figure(3);
[U2, omega2] = dtft(u2, T*3);
plot(omega2, abs(U2), 'r');
hold on;
[Y, omega] = dtft(decimate(u2, 3), T*3);
plot(omega/3, abs(Y));
hold off;
%}

N = 16;
p = 16;
x0 = [sin(2*pi/sqrt(17)*(0:15))];
x1 = [x0 zeros(1, N*(p-1))];

figure(1);
plot((1:p:length(x1)), abs(fft(x0)));
hold on;
plot((1:1:length(x1)), abs(fft(x1)), 'r');
hold off;