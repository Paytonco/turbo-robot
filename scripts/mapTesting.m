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
