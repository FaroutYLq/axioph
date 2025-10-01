function [Ec,Ej] = computeEcEj(Csigma, deltaL, Rn)
    planck = 4.135E-15; % eV.s
    ec = 1.6E-19;  % coulomb, also eVtoJ
    Ec=ec./(2.*(Csigma)); % eV
    Ej= planck.*deltaL./8./(ec.*Rn); % eV