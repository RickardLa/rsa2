%%% Lab 2 
%% Task 1
clc
clf
close all
clear all

% 1 
A = [1 -1.5 0.64];
B = 1;
f0 = 0:0.001:0.5;
G = freqz(B,A,2*pi*f0);

plot(f0,20*log10(abs(G)),'LineWidth',3)
xlabel('Normalized Frequency')
ylabel('Magnitude (dB)')
title('True Spectrum')

% 2
N = 1024;
L = 50; 
f = 0:1/N:(N-1)/(2*N);
x = randn([N+L,1]);
xFilt = filter(B,A,x);
x = xFilt(L+1:end); 

X = fft(x,N);
P = X.*conj(X)/N;
P = P(1:N/2); 
figure;
plot(f,10*log10(abs(P)))
xlabel('Normalized Frequency')
ylabel('Magnitude (dB)')
title('Periodogram')

% 3
K = 4;
M = N/K;
xx = reshape(x,M,K);
XX = fft(xx,N);
PP = XX.*conj(XX)/M;
PB = mean(PP');
PB = PB(1:N/2)';
plot(f,10*log10(abs(PB)),'g')
hold on 

K = 16;
M = N/K;
xx = reshape(x,M,K);
XX = fft(xx,N);
PP = XX.*conj(XX)/M;
PB = mean(PP');
PB = PB(1:N/2)';
plot(f,10*log10(abs(PB)),'r')
xlabel('Normalized Frequency')
ylabel('Magnitude (dB)')
title('True Spectrum vs. Periodogram Averaging')
legend('True Spectrum', 'Bartlett큦 method, K = 4', 'Bartlett큦 method, K = 16')


%% Task 2 
clc
clf
close all
clear all


% 1 
A = [1 -1.5 0.64];
B = 1;

N = 1024;
L = 50; 
f = 0:1/N:(N-1)/(2*N);
x = randn([N+L,1]);
xFilt = filter(B,A,x);
x = xFilt(L+1:end);
window = hamming(64);
PW = pwelch(x,window,[],2*pi*f)*2*pi; 
plot(f,10*log10(abs(PW)))
xlabel('Normalized Frequency')
ylabel('Magnitude (dB)')
title('Spectral Estimation using Hamming Window (M = 64)')

%% Task 3 
clc
clf
clear all 
close all


A = [1 -1.5 0.64];
B = 1;
N = 1024;
L = 50; 

x = randn([N+L,1]);
xFilt = filter(B,A,x);
x = xFilt(L+1:end); 

NFFT = 1024; 
M = 20;                         % Choose suitable number of lags

[PBT,fgrid] = btmethod(x,M,NFFT);
plot(fgrid,PBT)
xlabel('Normalized Frequency')
ylabel('Magnitude (dB)')
title('Spectral Estimation using Blackman-Tukey큦 Method')




%% Task 4
clc
clf
clear all
close all


A = [1 -1.5 0.64];
B = 1;

N = 1024;
L = 50; 

x = randn([N+L,1]);
xFilt = filter(B,A,x);
x = xFilt(L+1:end);

NFFT = 1024; 
f = 0:1/NFFT:(NFFT-1)/(2*NFFT);
p = 100;                              % Order

PAR = pyulear(x,p,NFFT)*pi;
PAR = PAR(1:NFFT/2);

plot(f,10*log10(abs(PAR)))
xlabel('Normalized Frequency')
ylabel('Magnitude (dB)')
title('Spectral Estimation using Parametric AR Modeling (p = 100)')



%% Task 5
clc
clf
clear all
close all

% 1 

A = [1 -0.24 0.08 -0.37 0.52];
B = [1 0.56 0.81]; 

f0 = 0:0.001:0.5;
G = freqz(B,A,2*pi*f0);

plot(f0,20*log10(abs(G)))
hold on

% 2

N = 8192;
L = 50; 

x = randn([N+L,1]);
xFilt = filter(B,A,x);
x = xFilt(L+1:end); 

NFFT = 2048; 
M = 30;                         % Choose suitable number of lags

[PBT,fgrid] = btmethod(x,M,NFFT);
plot(fgrid,PBT)
xlabel('Normalized Frequency')
ylabel('Magnitude (dB)')
title('Spectral Estimation using Blackman-Tukey큦 Method (M = 30)')
legend('True Sprectrum', 'Estimated Spectrum')


% 3
p = 8;                              % Order

PAR = pyulear(x,p,NFFT)*pi;
PAR = PAR(1:NFFT/2);
f = 0:1/NFFT:(NFFT-1)/(2*NFFT);

plot(f,10*log10(abs(PAR)))
xlabel('Normalized Frequency')
ylabel('Magnitude (dB)')
title('Spectral Estimation using Parametric AR Modeling (p = 8)')
legend('True Sprectrum', 'Estimated Spectrum')












