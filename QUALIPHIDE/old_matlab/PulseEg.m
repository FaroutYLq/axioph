function [t, pulsetrain] = PulseEg(Ns, fs, t0, tau )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
dt = 1/fs;
t = linspace(0,(Ns-1)*dt,Ns);
pulsetrain0 = exp(-(t-t0)/tau);
pulsetrain0(t<t0)=0;
pulsetrain = smoothdata(pulsetrain0,"gaussian",7);
% figure(9)
% plot(t, pulsetrain, 'Marker','+')
end