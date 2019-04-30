% Hodgkin Huxley Model

Kidx = 1;
VK = zeros(500001,5);

for K = [32, 34, 35.7, 35.9, 36]

gbar_Na = 120; % mS/cm^3
gbar_K = K; % mS/cm^3 = 36 % modify for refractory period
gbar_L = 0.3; % mS/cm^3

C_m = 1; %uF/cm^3

E_Na = 50; % mV
E_K = -77; % mV
E_L = -54.4; % mV

V = -70; % mV
Vm = -70; % mV

m = 0.01;
n = 0.01;
h = 0.01;

T = 20; % deg C
Tb = 6.3; % deg C
Q = 3; % giant squid axon
phi = Q^((T-Tb)/10);

%{
a_n = 0.01*(V+55)/(1-exp(-(V+55)/10));
B_n = 0.125*exp(-(V+65)/80);
a_m = 0.1*(V+40)/(1-exp(-(V+40)/10));
B_m = 4*exp(-(V+65)/18);
a_h = 0.07*exp(-(V+65)/20);
B_h = 1/(1+exp(-(V+35)/10));

n_inf = a_n/(a_n+B_n);
tau_n = 1/(a_n+B_n);
m_inf = a_m/(a_m+B_m);
tau_m = 1/(a_m+B_m);
h_inf = a_h/(a_h+B_h);
tau_h = 1/(a_h+B_h);

dn = dt*(n_inf-n)/tau_n;
dm = dt*(m_inf-m)/tau_m;
dh = dt*(h_inf-h)/tau_h;

dV = dt*(1/C_m)*(-gbar_Na*m^3*h*(V-E_Na)-gbar_K*n^4*(V-E_K)-gbar_L*(V-E_L));
%}

dt = 0.001; % 100 us
dur = 500; % 30 s = 30e6 us
j = 1;
I = zeros(dur/dt+1,1);
I(1:5) = 10;
I(50000:50010) = 600;
I(100000:100010) = 600;
I(150000:150010) = 600;
I(200000:200010) = 600;
I(250000:250010) = 600;
I(300000:300010) = 600;
I(350000:350010) = 600;
I(400000:400010) = 600;

I(370000:370010) = 600;
Vm = zeros(dur/dt,1);

for i = 0:dt:dur
    
    a_n = 0.01*(V+55)/(1-exp(-(V+55)/10));
    B_n = 0.125*exp(-(V+65)/80);
    a_m = 0.1*(V+40)/(1-exp(-(V+40)/10));
    B_m = 4*exp(-(V+65)/18);
    a_h = 0.07*exp(-(V+65)/20);
    B_h = 1/(1+exp(-(V+35)/10));

    n_inf = a_n/(a_n+B_n);
    tau_n = 1/(a_n+B_n);
    m_inf = a_m/(a_m+B_m);
    tau_m = 1/(a_m+B_m);
    h_inf = a_h/(a_h+B_h);
    tau_h = 1/(a_h+B_h);

    dn = dt*(n_inf-n)/tau_n;
    dm = dt*(m_inf-m)/tau_m;
    dh = dt*(h_inf-h)/tau_h;
    
    n = n + dn;
    m = m + dm;
    h = h + dh;
    
    I_L = gbar_L * (V - E_L);
    I_Na = gbar_Na * m^3 * h * (V - E_Na);
    I_K = gbar_K * n^4 * (V - E_K);
    I_Cm = I(j) - (I_L + I_Na + I_K);
    
    %dV = dt*(1/C_m)*(-gbar_Na*m^3*h*(V-E_Na)-gbar_K*n^4*(V-E_K)-gbar_L*(V-E_L));
    
    dVm = dt * I_Cm / C_m;
    V = V + dVm;
    %Vm = [Vm; V];
    Vm(j) = V;
    
    j = j + 1;
    
end

VK(:,Kidx) = Vm;
Kidx = Kidx + 1;    
end

plot(VK)
