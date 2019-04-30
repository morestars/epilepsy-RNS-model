
%% Load in data
%[hdr, record] = edfread('chb01_01.edf');
[hdr, chb01_03] = edfread('chb01_03.edf');
[hdr, chb01_04] = edfread('chb01_04.edf');

[hdr, chb03_01] = edfread('chb03_01.edf');
[hdr, chb03_02] = edfread('chb03_02.edf');
[hdr, chb03_03] = edfread('chb03_03.edf');
[hdr, chb03_04] = edfread('chb03_04.edf');
[hdr, chb03_34] = edfread('chb03_34.edf');
[hdr, chb03_35] = edfread('chb03_35.edf');
[hdr, chb03_36] = edfread('chb03_36.edf');

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

%%
histogram((seizureTimesLineLength),20)
hold on
histogram((seizureTimesAbsArea),20)
hold off

%%
plot(hitLineLength)
hold on
plot(hitAbsArea)
hold off

