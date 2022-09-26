%% Finding delay and confirming 5mV resolution
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