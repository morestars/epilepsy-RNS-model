
%% Load in data
%[hdr, record] = edfread('chb01_01.edf');
[hdr, chb03_01] = edfread('chb03_01.edf');
[hdr, chb03_02] = edfread('chb03_02.edf');
[hdr, chb03_03] = edfread('chb03_03.edf');
[hdr, chb03_04] = edfread('chb03_04.edf');
[hdr, chb03_34] = edfread('chb03_34.edf');
[hdr, chb03_35] = edfread('chb03_35.edf');
[hdr, chb03_36] = edfread('chb03_36.edf');

%% Analyze data
f = 256;
dt = 1/f;
x = 0:1/f:(1/f)*(length(record)-1);

ecog = record03(1,1:30);
dLen = sqrt( (diff(ecog).^2 + dt^2 )); % pythagorean distance
totLen = sum(dLen);
prevLen = totLen;
hitLineLength = [];

for i = 31:30:length(record03)
    ecog = record03(1,i:i+29);
    [hit, prevLen] = lineLength(ecog,1/f,prevLen);
    hitLineLength = [hitLineLength hit];
end

seizureTimes = find(hitLineLength)*30/f;