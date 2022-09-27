%% Determining qualities of input signal    BP TESTING
%determine frequency, amplitude, and mean (report all as mean \pm std)
% e.g. IIIA5mV.lvm

[frequency, amplitude, mean] = inChar('../data/IIIA5V.lvm') %input name of file of interest


%freq is principal frequency, A is a tuple with [mean amplitude, std
%amplitude] and mu is the midline value. frequency in Hz, amplitude +
%midline in V. 


%% Checking that filter did what we intended + Characterizing the Transfer Function MAP TESTING
%Reproduce results in MATLAB
% e.g. IIIA10.lvm

for N = [10, 100, 1000, 10000]  %iterate over the number of points used
    filename = "../data/IIIA" + num2str(N) + ".lvm";    %directory to data file
    mapArray = load(filename);  %load data from data file
    x = mapArray(:,3);  %extract x vector (not ordered according to lab instructions)
    y = mapArray(:,2);  %extract y vector
    newY = zeros(size(y));  %initialize MATLAB replication of calculated MAP
    for i = 1:length(y)     %populate new calculated MAP vector
        window = y(max(i-N+1,1):i); 
        newY(i) = sum(window)/length(window);
    end

    M = length(x);              %signal length
    fs = 1000;                  %sampling frequency (Hz)
    dt = 1 / fs;                 %sampling period
    t = (0:M-1) * dt;           %Time array (s)

    figure(N)                       %Plot original and recreated MAP
    plot(t,x,t,newY, "LineWidth", 2)
    legend("LabVIEW", "MATLAB")
    xlabel("Time (s)")
    ylabel("Arterial Pressure (V)")
    title("N="+num2str(N)+" Filter Verification")

    figure(N + 1)                   %Plot measured bp against calculated MAP
    plot(t,y,t,x,"LineWidth", 2)
    legend("Measured BP", "Calculated MAP")
    xlabel("Time (s)")
    ylabel("Arterial Pressure (V)")
    title("N="+num2str(N)+" Performance")
end





%% Finding delay and confirming 5mV resolution  OUTPUT TESTING
%Find delay as argmax xcov, take delayed difference, form CI
% e.g. IIIAZ5mV.lvm
clear;clc;
outputArray = load("../data/IIIAZ5V.lvm");      %load recorded data
y = outputArray(:,1);                           %extract columns
x = outputArray(:,2);
z = outputArray(:,3);

N = length(x);              %signal length
fs = 1000;                  %sampling frequency (Hz)
dt = 1 / fs;                 %sampling period (s)
df = fs/N;                  %frequency step size (Hz)
f = (0:N-1) * df;           %frequency array (Hz)
t = (0:N-1) * dt;           %Time array (s)

mu_y = mean(y);             %reset y to have midline at zero
y = y - mu_y;
mu_z = mean(z);             %reset z to have midline at zero
z = z - mu_z;

offset = mu_z - mu_y

delay = finddelay(y,z);         %use cross correlation to find most likely delay (indices)
delayTime = delay * dt          %most likely delay (s)

gain = mean(z(delay+1:N) ./ y(1:N-delay))



%% 90 bpm VALIDATION
clear, clc

data = sprintf('../data/IIIB90bpm10000.lvm'); % reading in 90bpm N = 10,000 data
[x y z] = textread(data,'%n %n %n','headerlines',23);
% data1 = sprintf('IIIB90bpm1000.lvm'); % reading in 90bpm N = 1,000 data
% [x1 y1 z1] = textread(data1,'%n %n %n','headerlines',23);
% data2 = sprintf('IIIB90bpm100.lvm'); % reading in 90bpm N = 100 data
% [x2 y2 z2] = textread(data2,'%n %n %n','headerlines',23);

t = 1:2000; % time in milliseconds


xa = x(9001:11000); % N = 10,000 x data
ya = y(9001:11000); % N = 10,000 y data

x1 = x(8991:10990); % shift data back in time by 10 ms
y1 = y(8991:10990); 

x2 = x(9011:11010); % shift data forward in time by 10 ms
y2 = y(9011:11010); 


length(t)
length(x2)
figure(1)
subplot(2,1,1);
pa = plot(t/1000,ya*100,'-k',t/1000,xa*100,'-b','LineWidth',2);
hold on 
pa1 = plot(t/1000,x1*100,'-g','LineWidth',2);
pa2 = plot(t/1000,x2*100,'-c','LineWidth',2);
grid on;
legend('Arterial Pressure','T_{MAP} = 667 ms','T_{MAP} = 657 ms','T_{MAP} = 677 ms');
xlabel('Time (seconds)')
ylabel('Pressure (mmHg)')
hold off



subplot(2,1,2)
pb = plot(t/1000,ya*100,'-b','Linewidth',2);
hold on
pb1 = plot(t/1000,y1*100,'-g','Linewidth',2);
pb2 = plot(t/1000,y2*100,'-c','Linewidth',2);
axis([0 2 72 75]);
grid on;
xlabel('Time (seconds)')
ylabel('MAP (mmHg)')
hold off
%% 180 bpm 
clear, clc
figure(2)


t = 1:2000; % time

data3 = sprintf('../data/IIIB180bpm10000.lvm'); % reading in 180bpm N = 10,000 data
[X Y Z] = textread(data3,'%n %n %n','headerlines',23);


x3 = X(10001:12000); % N = 10,000 x data
y3 = Y(10001:12000); % N = 10,000 y data

figure(2)
subplot(2,1,1)
pc = plot(t/1000, y3*100,'-b', t/1000, x3*100, '-r','Linewidth', 2)
legend('Arterial Pressure','T MAP = 667 sec')
grid on;
title('90 Beats per Minute')
xlabel('Time (seconds)')
ylabel('Pressure (mmHg)')


subplot(2,1,2)
pd = plot(t/1000,y3*100,'-b','LineWidth',2);
axis([0 2 72 75]); grid on;
ylabel('MAP, mmHg');
xlabel('Time, s');


%% 60 bpm 
clear, clc

data3 = sprintf('../data/IIIB60bpm10000.lvm'); % reading in 180bpm N = 10,000 data
[X Y Z] = textread(data3,'%n %n %n','headerlines',23);

t = 1:2000; % time (miliseconds)
x = X(10001:12000); % N = 10,000 x data
y = Y(10001:12000); % N = 10,000 y data

subplot(2,1,1)
pe = plot(t/1000, y*100, '-r', t/1000, x*100, '-b', 'Linewidth', 2)
grid on;
legend('Arterial Pressure','T_{MAP} = 667 ms');
xlabel('Time (seconds)')
ylabel('Pressure (mmHg)')

subplot(2,1,2)
pd = plot(t/1000,y*100,'-b','LineWidth',2);
axis([0 2 72 75]); grid on;
ylabel('MAP, mmHg');
xlabel('Time (seconds)');



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





