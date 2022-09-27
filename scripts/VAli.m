%% 90 bpm 
clear, clc

data = sprintf('IIIB90bpm10000.lvm'); % reading in 90bpm N = 10,000 data
[x y z] = textread(data,'%n %n %n','headerlines',23);

t = 1:2000; % time in milliseconds


xa = x(10001:12000); % N = 10,000 x data
ya = y(10001:12000); % N = 10,000 y data

x1 = x(9991:11990); % shift data back in time by 10 ms
y1 = y(9991:11990); 

x2 = x(10011:12010); % shift data forward in time by 10 ms
y2 = y(10011:12010); 


figure(1)
subplot(2,1,1);
set1 = plot(t/1000,ya*100,'-k',t/1000,xa*100,'-b','LineWidth',2);
hold on 
set2 = plot(t/1000,x1*100,'-g','LineWidth',2);
set3 = plot(t/1000,x2*100,'-c','LineWidth',2);
grid on;
legend('Arterial Pressure','T_{MAP} = 667 ms','T_{MAP} = 657 ms','T_{MAP} = 677 ms');
xlabel('Time (seconds)')
title('90 Beats per Minute')
ylabel('Pressure (mmHg)')
hold off



subplot(2,1,2)
pa = plot(t/1000,ya*100,'-k',t/1000,xa*100,'-b','LineWidth',2);
hold on 
set2 = plot(t/1000,x1*100,'-g','LineWidth',2);
set3 = plot(t/1000,x2*100,'-c','LineWidth',2);
axis([0 2 72 78]); % zooming axes of above graph
grid on;
xlabel('Time (seconds)')
ylabel('MAP (mmHg)')
hold off
%% 180 bpm 
clear, clc
figure(2)

t = 1:2000; % time

data3 = sprintf('IIIB180bpm10000.lvm'); % reading in 180bpm N = 10,000 data
[X Y Z] = textread(data3,'%n %n %n','headerlines',23);


x3 = X(10001:12000); % N = 10,000 x data
y3 = Y(10001:12000); % N = 10,000 y data

figure(2)
subplot(2,1,1)
set1 = plot(t/1000, y3*100,'-b', t/1000, x3*100, '-r','Linewidth', 2)
legend('Arterial Pressure','T_{MAP} = 667 sec')
grid on;
title('180 Beats per Minute')
xlabel('Time (seconds)')
ylabel('Pressure (mmHg)')


subplot(2,1,2)
set1 = plot(t/1000, y3*100,'-b', t/1000, x3*100, '-r','Linewidth', 2)
% pd = plot(t/1000,y3*100,'-b','LineWidth',2);
axis([0 2 72 78]); grid on;
ylabel('MAP, mmHg');
xlabel('Time, s');


%% 60 bpm 
clear, clc

data3 = sprintf('IIIB60bpm10000.lvm'); % reading in 180bpm N = 10,000 data
[X Y Z] = textread(data3,'%n %n %n','headerlines',23);

t = 1:2000; % time (miliseconds)
x = X(10001:12000); % N = 10,000 x data
y = Y(10001:12000); % N = 10,000 y data

subplot(2,1,1)
set1 = plot(t/1000, y*100, '-r', t/1000, x*100, '-b', 'Linewidth', 2)
grid on;
legend('Arterial Pressure','T_{MAP} = 667 ms');
xlabel('Time (seconds)')
ylabel('Pressure (mmHg)')
title('60 Beats per Minute')

subplot(2,1,2)
set1 = plot(t/1000, y*100, '-r', t/1000, x*100, '-b', 'Linewidth', 2)
% pd = plot(t/1000,y*100,'-b','LineWidth',2);
axis([0 2 72 78]); grid on;
ylabel('MAP, mmHg');
xlabel('Time (seconds)');







