%% Determining qualities of input signal    BP TESTING
%determine frequency, amplitude, and mean (report all as mean \pm std)
% e.g. IIIA5mV.lvm

[frequency, amplitude, mean] = inChar('../data/IIIA5V.lvm') %input name of file of interest


%freq is principal frequency, A is a tuple with [mean amplitude, std
%amplitude] and mu is the midline value. frequency in Hz, amplitude +
%midline in V. 
function [freq, A, mu] = inChar(filename)
    bpArray = load(filename);   %load .lvm recorded by LabVIEW
    x = bpArray(:,1);           %input (blood pressure) vector
    y = bpArray(:,2);           %Calculated MAP vector

    mu = mean(x);               %calculate mean value
    x = x - mu;                 %recenter x at midline for amplitude calculations
    
    N = length(x);              %signal length
    fs = 1000;                  %sampling frequency (Hz)
    dt = 1 / fs;                 %sampling period
    df = fs/N;                  %frequency step size (Hz)
    f = (0:N-1) * df;           %frequency array (Hz)
    t = (0:N-1) * dt;           %Time array (s)
    X = abs(1/N * fft(x));      %find fft of x to determine amplitude and frequency

    figure(1)                   %plotting input waveform
    plot(t(t>0.5& t<0.6),x(t>0.5 & t<0.6), 'LineWidth',2)
    xlabel('Time (s)')
    ylabel('Amplitude (V)')
    title('5V Square Wave, 50Hz')

    figure(2)                   %plotting amplitdue spectral density
    plot(f,X)
    xlabel('Frequency (Hz)')
    ylabel('Amplitude (V)')
    title('5V Square Wave, 50Hz, ASD')

    [~, maxIndex] = max(X(1:round(N/2)));   %Find index of primary frequency
    freq = f(maxIndex);                     %Find primary frequency

    pks = findpeaks(x);                 %find local maxima in signal
    A = [mean(pks), std(pks)];          %Find mean + std of local maxima
end
