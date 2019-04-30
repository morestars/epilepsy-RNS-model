%% NII Final Project
% Epilepsy
%% Load in data
%[hdr, record] = edfread('chb0x_yy.edf');
[hdr, chb03_34] = edfread('chb03_34.edf');

%% Analyze data
record = chb03_34;
f = 256;
dt = 1/f;
x = 0:1/f:(1/f)*(length(record)-1);

%% Line Length
step = 100;
channel = 10;
ecog = record(channel,1:step);
dLen = sqrt( (diff(ecog).^2 + dt^2 )); % pythagorean distance
totLen = sum(dLen);
prevLen = totLen;
hitLineLength = [];

for i = step+1:step:length(record)
    ecog = record(channel,i:i+step-1);
    [hit, prevLen] = lineLength(ecog,1/f,prevLen);
    hitLineLength = [hitLineLength hit];
end

seizureTimesLineLength = find(hitLineLength)*step/f;

%% Area
step = 100;
ecog = record(channel,1:step);
prevAbsArea = trapz(abs(ecog))*dt;
hitAbsArea = [];


for i = step+1:step:length(record)
    ecog = record(channel,i:i+step-1);
    [hit, prevAbsArea] = absArea(ecog,1/f,prevAbsArea);
    hitAbsArea = [hitAbsArea hit];
end

seizureTimesAbsArea = find(hitAbsArea)*step/f;

%% Visualize detection
figure()
plot(x,record')
xlabel('Time (s)')
ylabel('23-Channel EEG')
title('Seizure at 1982-2029s')

figure()
histogram((seizureTimesLineLength))
hold on
histogram((seizureTimesAbsArea))
xlabel('Time (s)')
ylabel('Frequency')
legend('Line Length Detection','Area Detection')
hold off

x2 = 0:1/f*step:(1/f)*step*(length(hitLineLength)-1);
figure()
plot(x2,hitLineLength,'LineWidth',3)
hold on
plot(x2,hitAbsArea,'LineWidth',3)
xlabel('Time (s)')
ylabel('Detection (Binary)')
legend('Line Length Detection','Area Detection')
hold off

%% Combine detection
dtctTime = [];
for i = seizureTimesAbsArea
    for j = seizureTimesLineLength
        if round(i) == round(j)
            dtctTime = [dtctTime round(i)];
        end
    end
end

%% Stimulation
lead = dtctTime*f;
%lead = find(x==dtctTime);
lag = length(record)-lead-520;
pulse = stimPattern(lead,lag,5,0.0005,0.005,100);

figure()
subplot(3,1,1)
plot(x,record')
xlabel('Time (s)')
ylabel('23-Channel EEG')
title('Seizure at 1982-2029s')

subplot(3,1,2)
plot(x2,hitLineLength,'LineWidth',3)
hold on
plot(x2,hitAbsArea,'LineWidth',3)
xlabel('Time (s)')
ylabel('Detection (Binary)')
legend('Line Length Detection','Area Detection')
hold off

subplot(3,1,3)
plot(x,pulse,'LineWidth',3)
xlabel('Time (s)')
ylabel('Pulse (mA)')

            


