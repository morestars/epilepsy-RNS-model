function filter = freqdetect(ecog,fs,lc,hc)
%ecog data
%fs sampling rate
%lowpass cut off
%high pass cut off

filter = bandstop(ecog, [lc hc],fs);
Thresh = 15;
filter(filter>=Thresh) = 1;
filter(filter<Thresh)=0;
end
