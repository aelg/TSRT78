clear all

% Load and detrend data.
Fs=8000;
load sentence.mat
sen = detrend(y(2,:)*100);
len = length(sen);
x = 0:2/len:1.9999999999;

% Create AR models of all the segments
na = 8; % Model order
cov = [];
est = zeros(length(sen),1);
delay = 0;
Z = []; % Save filter coefficents to avoid steps in signal.
for i = 1:160:length(sen)-160
    seg = i:i+159;
    m = ar(detrend(sen(seg)), na);
    e = filter(m.a,1,detrend(sen(seg))'); % m1 <-> AR model of the segment
    r = abs(covf(e,100)); % Get the residual covariance.
    cut = 19; % Used to not make the period to short.
    r = r(cut+1:end);
    cov = [cov; r];
    [ma, ind] = max(r); % Find max.
    A = ma; % Set amplitude.
    pulse = ind+cut; % Set pulse length.
    
    % Detect and fix unstable poles.
    if max(abs(roots(m.a))) >= 1
        % Print warning.
        disp(['Unstable: abs = ' num2str(max(abs(roots(m.a))))]);
        
        r = roots(m.a); % Get roots
        r2 = [];
        for j = r(:)
            if abs(j) >= 1
                r2 = [r2; 1/j]; % Reflect if outside unit circle.
            else
                r2 = [r2; j];
            end
        end
        m.a = poly(r2); % Put the fixed polynom back.
    end
    
    %A = 1;
    
    % Create pulse train.
    train = sqrt(A)*(rem((1:160)+delay,pulse) == 1)';
    train = 0.3*sqrt(mean(r))*randn(160, 1)+train;
    
    % Make delays overlap nicely between segements.
    delay = rem(160+delay, pulse);
    
    % Z Makes the signal smoother between segements.
    [est(seg), Z] = filter(1,m.a,train, Z);
end
% Generate sound file
wavwrite(est, 'out.wav');