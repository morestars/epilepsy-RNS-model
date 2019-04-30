function [V] = LPF(R,C,Vi,dt)
    % R, C, dt values:
        %{ 
        R = 1000;
        C = 300e-6;
        dt = 0.01;
        %}

    T = [];
    V = [];
    I = [];
    Vo = 0;
    Io = 0;
    i = 1;
    dur = length(Vi)*dt-dt;

    for t=0:dt:dur
        Io = (Vi(i)-Vo)/R;
        dVo=dt*((Vi(i)-Vo)/(R*C));
        Vo=Vo+dVo;
        I = [I,Io];
        V=[V,Vo];
        T=[T,t];
        i = i+1;
    end
end