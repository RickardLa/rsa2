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
hold on
% xlabel('Frequency')
% ylabel('Magnitude (dB)')
% title('True Spectrum')

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
% plot(f,10*log10(abs(P)))
% xlabel('Frequency')
% ylabel('Magnitude (dB)')
% title('Periodogram')

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
xlabel('Frequency')
ylabel('Magnitude (dB)')
title('True Spectrum vs. Periodogram Averaging')
legend('True Spectrum', 'Bartlett´s method, K = 4', 'Bartlett´s method, K = 16')


%% Task 2 
clc
clf
close all
clear all


% 1 
A = [1 -1.5 0.64];
B = 1;
f0 = 0:0.001:0.5;
G = freqz(B,A,2*pi*f0);

N = 1024;
L = 50; 
f = 0:1/N:(N-1)/(2*N);
x = randn([N+L,1]);
xFilt = filter(B,A,x);
x = xFilt(L+1:end);
window = hamming(64);
PW = pwelch(x,window,[],2*pi*f)*2*pi; 
plot(f,10*log10(abs(PW)))
xlabel('Frequency')
ylabel('Magnitude (dB)')
title('Spectral Estimation with Hamming Window (M = 64)')
































