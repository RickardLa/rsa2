

 
function [PBT,fgrid] = btmethod(x,M,NFFT)

window = hamming(2*M+1);                % Generate Hamming window 
rx = xcorr(x,M,'biased');               % Compute sample autocorrelation function
fgrid = 0:1/NFFT:(NFFT-1)/(2*NFFT);     % The frequency grid it should be computed at

BT = fft(window.*rx,NFFT);              % Compute the FFT of the product 
BT = BT(1:NFFT/2);                      % Symmetric -> Only need half
PBT = 10*log10(abs(BT));                % Convert to dB 

end