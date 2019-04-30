
T = 200;
dT=.01;
t=0:dT:T;


clear Vend temp
currentLevels = 10; 

I(1:numel(t)) = currentLevels;




gbar_K=36; gbar_Na=120; g_L=.3;
E_K = -12; E_Na=115; E_L=10.6;
C=1;

for k = 1:20;
        Iimp(1:numel(t)) = 0;
        Iimp(9200:9300)=-100*(1/sqrt((14-k)^2+1))-100*(1/sqrt((7-k)^2+1));
       Iimp(9500:9600)=100*(1/sqrt((14-k)^2+1))+100*(1/sqrt((7-k)^2+1));

   


V=0;
alpha_n = .01 * ( (10-V) / (exp((10-V)/10)-1) );
beta_n = .125*exp(-V/80);
alpha_m = .1*( (25-V) / (exp((25-V)/10)-1) );
beta_m = 4*exp(-V/18); 
alpha_h = .07*exp(-V/20); 
beta_h = 1/(exp((30-V)/10)+1);

n(1) = alpha_n/(alpha_n+beta_n);
m(1) = alpha_m/(alpha_m+beta_m); 
h(1) = alpha_h/(alpha_h+beta_h); 


for i=1:numel(t)-1 
    alpha_n(i) = .01 * ( (10-V(i)) / (exp((10-V(i))/10)-1) );
    beta_n(i) = .125*exp(-V(i)/80);
    alpha_m(i) = .1*( (25-V(i)) / (exp((25-V(i))/10)-1) );
    beta_m(i) = 4*exp(-V(i)/18);
    alpha_h(i) = .07*exp(-V(i)/20);
    beta_h(i) = 1/(exp((30-V(i))/10)+1);
    

    I_Na = (m(i)^3) * gbar_Na * h(i) * (V(i)-E_Na); 
    I_K = (n(i)^4) * gbar_K * (V(i)-E_K); 
    I_L = g_L *(V(i)-E_L); 
    I_ion = I(i)+ Iimp(i)- I_K - I_Na - I_L; 
    
    V(i+1) = V(i) + dT*I_ion/C;
    n(i+1) = n(i) + dT*(alpha_n(i) *(1-n(i)) - beta_n(i) * n(i)); %Equation 7
    m(i+1) = m(i) + dT*(alpha_m(i) *(1-m(i)) - beta_m(i) * m(i)); %Equation 15
    h(i+1) = h(i) + dT*(alpha_h(i) *(1-h(i)) - beta_h(i) * h(i)); %Equation 16

end


V = V-70;
Vend(k,:)=V;
end
Thresh = 20;
Vp = Vend;
Vend(Vend<Thresh)=0;
Vend(Vend>=Thresh) = 1;
for o = 1:k
    Vend(o,:)=Vend(o,:)*o;
end
Vend(Vend==0)=NaN;

figure()
plot(t,Vp)
ylabel('Membrane Potential (mV)')
xlabel('Time (ms)')
title('Action Potential Desynchronization')

figure()
plot(t,Vend,'.')
hold on
Iimp = (Iimp/10)+7;
plot(t,Iimp)
Iimp = Iimp+7;
ylabel('Neuron')
xlabel('Time (ms)')
title('Raster Plot')
plot(t,Iimp)


