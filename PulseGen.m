function pulse = PulseGen(filter,stimamp,stimdur)
pulse = zeros(length(ecog));
intperiod = 10;
threshold = 5;
for i =1:intperiod: ength(filter);
    det = trapz(i:i+intperiod);
    if det > threshold
        pulse(i:i+stimdur)=ones(stimdur)*stimamp;
    end
end
end
