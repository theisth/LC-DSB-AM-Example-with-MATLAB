% Ahmet Ali Tilkicioğlu / ELEC361 Project / DSB-LC-AM 210102002163
clc
clear all

t0 = .1;                                        % 0.1 second variable
ts = 0.0001;                                    % sample time variable (second)  
t = 0:ts:t0-ts;                                 % 0 to 0.1 sec. sampled time
fs = 1/ts;                                      % sample frequency
df = 1/t(end);                                  % sampled frequency values

m = 5*cos(2*pi*100.*t) + 10*cos(2*pi*200.*t);   % message signal m(t)
c = 2*cos(2*pi*1000.*t);                        % carrier signal c(t)
x = m + c;                                      % x signal -> x(t) = m(t) + c(t)
z = (x.*x) + 60.*x;                             % z signal -> z(t) = x^2(t) + 60x(t)


M = fft(m);                                     % Fourier Transform m(t) -> M(f)
X = fft(x);                                     % Fourier Transform x(t) -> X(f)
Z = fft(z);                                     % Fourier Transform z(t) -> Z(f)

shifted_M = abs(fftshift(M)/length(M));         % M(f) value shifted real frequency value
shifted_X = abs(fftshift(X)/length(X));         % X(f) value shifted real frequency value
shifted_Z = abs(fftshift(Z)/length(Z));         % Z(f) value shifted real frequency value

L = length(m);                                  % Lenght of Message signal
f = fs/length(M)*(-length(M)/2:length(M)/2-1);  % sampled all frequency values

% Band Pass Filter Settings
gain = 1;                                       % gain value
bandwidth = 400;                                % bandwidth value
center_Frequency = 1000;                        % center frequency value

%BPF frequency limit calculations
freq_Low_BPF = center_Frequency-(bandwidth/2);  %band pass filter low frequency
freq_High_BPF = center_Frequency+(bandwidth/2); %band pass filter high frequency


r_df = round(df);                               % for getting integer value

% create equal Z(f) vector list with zeros
Y = zeros(1,L);
% positive frequency side for BPF
Y(freq_Low_BPF/r_df:freq_High_BPF/r_df+1) = Z(freq_Low_BPF/r_df:freq_High_BPF/r_df+1)* gain;   
% negative frequency side for BPF
Y((fs/r_df)-(freq_High_BPF/r_df):(fs/r_df)-(freq_Low_BPF/r_df)+1) = Z((fs/r_df)-(freq_High_BPF/r_df):(fs/r_df)-(freq_Low_BPF/r_df)+1) * gain;

% WARNING: Vector list is not 0 to 1000. it is 1 to 1000, so we must add 1
% last limit to get spectrum value in 1200 and -1200 hz.


shifted_Y = abs(fftshift(Y)/length(Y));         % Z(f) value shifted real frequency value
y = ifft(Y);                                    % inverse fourier transform


env = envelope(real(y));                        % envelope y(t)
%WARNING: Input must be real, so we used that real function. 

ft_env = fft(env);                              %envelope fourier transform
shifted_ft_env = abs(fftshift(ft_env)/length(ft_env));  %  envelope func. shifted real frequency value



% PLOT m(t), c(t) ,z(t) ,y(t)
figure(1)

subplot(2,2,1)
plot(t,m(1:length(t)))                          % plotting m(t)
xlabel('time (s)')
ylabel('m(t)')
title('m signal')

subplot(2,2,2)
plot(t,x(1:length(t)))                          % plotting x(t)
xlabel('time (s)')
ylabel('y(t)')
title('x signal')

subplot(2,2,3)
plot(t,z(1:length(t)))                          % plotting z(t)
xlabel('time (s)')
ylabel('z(t)')
title('z signal')


subplot(2,2,4)
plot(t,y(1:length(t)))                          % plotting y(t)
xlabel('time (s)')
ylabel('y(t)')
title('y signal')


% PLOT M(f), C(f) ,Z(f) ,Y(f)
figure(2)

subplot(2,2,1)
plot(f,shifted_M)                               %plotting M(f)
xlabel('frequency (hz)')
ylabel('M(f)')
title('M spectral')

subplot(2,2,2)
plot(f,shifted_X)
xlabel('frequency (hz)')                        %plotting X(f)
ylabel('X(f)')
title('X spectral')

subplot(2,2,3)
plot(f,shifted_Z)                               %plotting Z(f)
xlabel('frequency (hz)')
ylabel('Z(f)')
title('Z spectral')

subplot(2,2,4)
plot(f,shifted_Y)                               % plotting Y(f)
xlabel('frequency (hz)')
ylabel('Y(f)')
title('Y spectral')

figure(3)

subplot(2,1,1)

[up,lo] = envelope(real(y),100,'analytic');
plot(t,up,'-',t,lo,'--')                        % plotting envelope
xlabel('time (s)')
ylabel('~m(t)')
title('y(t) envelope signal')


subplot(2,1,2)
plot(f,shifted_ft_env)                               % plotting ~M(f)
xlabel('frequency (hz)')
ylabel('~M(f)')
title('~M spectral')



