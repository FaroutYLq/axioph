clear, clc, close all

%% Universal Constants

planck = 4.135E-15; % eV.s
kb = 8.617E-5; % eV/K
eVtoJ = 1.6E-19;
ec = 1.6E-19;  % coulomb

%% Aluminum specific
DOS = 1.72E10; % /um3.eV

fermi = 11.6; % eV Fermi level
Tc = 1.2; % AlMn 0.1, Al 1.2 

% Superconducting gap for OCS and also for "left" and "right" sides of QCD
delta = 1.89E-4; % eV

%% Experiment Fixed

T = 0.02; %K Temperature
Rn = 27E3; % ohms normal state resistance
cpl = 6.5e-17;% F/um gate capacitance per micron
clen = 20; % um gate capacitor length in micron
Cg = 1E-15;%clen*cpl; % F gate capacitance

% Computed quantities

e2R = ec*ec*Rn/eVtoJ; % I^2*R*t*t = e2R % E*t --> eV*s

%% Transmon if you want to specify capacitances
Cj = 4.43E-16; % F junction capacitance
Cshunt = 6.5E-14;
Csigma = Cj+Cg+Cshunt; % F total capacitance
[Ec, Ej] = computeEcEj(Csigma, delta, Rn);
curlyN = DOS*sqrt(2*pi*delta*kb*T); % /um3

%% Transmon if you want to directly specify Ej, Ec

ejecratio = 12;
Ej = 0.4*kb; % WashU
Ec = Ej/ejecratio;
%Ec = 0.015*kb; % 

%Ej = 0.295*kb; % Serniak
%Ec = 0.017*kb;

%Ej = 0.250*kb; % WIS target
%Ec = 0.022*kb;

%Ej = 0.27*kb; % WIS target 2
%Ec = 0.015*kb;

Rn = planck.*delta./8./(ec.*Ej);

disp(['Rn kOhm: ' num2str(Rn/1000)])
disp(['Ej GHz: ' num2str(Ej/planck/1E9)])
disp(['Ec GHz: ' num2str(Ec/planck/1E9)])
disp(['Transmon Ej/Ec ~ ' num2str(Ej/Ec)])

%% Energy diagram
% Diagonalizes Hamiltonian to generate qubit energy levels

u=linspace(0,1,500); % Offset Charge
[EE, EO, DE] = solvesystem(Ec,Ej,u, delta, delta);

EE_f = (EE - EE(:,1))./planck./1e9;
EO_f = (EO - EO(:,1))./planck./1e9;

disp(['E01 GHz: ' num2str(EO_f(1,2))])
disp(['E02 GHz: ' num2str(EO_f(1,3))])
disp(['E03 GHz: ' num2str(EO_f(1,4))])
disp(['Transmon dE/Ec ratio: ' num2str(DE(1)/Ec)]);
fr1 = min(EO_f(:,2));
fr3 = min(EO_f(:,4));

FigHandleA = figure;
set(FigHandleA, 'Position', [100, 100, 800, 600]);
plot(u,EE_f,'LineWidth',2);
ax = gca;
ax.ColorOrderIndex = 1;
hold on
plot(u,EO_f,'LineWidth',2,'LineStyle','--');
hold off
xlabel('Offset Charge [CgVg/2e]','Interpreter','latex','FontSize',25);
ylabel('f$_{0j}$ [GHz]','Interpreter','latex','FontSize',25);
title(['$E_J/E_C=$ ' num2str(Ej/Ec)],'Interpreter','latex','FontSize',25);
set(gca,'TickLabelInterpreter','latex','FontSize',25);
set(gca,'YMinorTick','on');
set(gca,'XMinorTick','on');

%% Dispersive Shift

% User specified
g = 150E6; % 150 MHz example, from HFSS
fr = 7.0E9; % Hz, where you place resonator
nlevels = 6;
matrixelems = zeros(length(u),nlevels);
chi_ip = zeros(length(u),nlevels);
for i=1:length(u)
    [fullmat,chi_ip(i,:)]=dispermatrix(Ec,Ej,u(i)+0.5,g,fr,nlevels);
    matrixelems(i,:) = fullmat(1,:);
end

disp(['fr: ' num2str(fr/1E9) ' GHz'])

fr_range = fr:1E6:((fr3*1E9)+2E8);
chi_scan = zeros(1,length(fr_range));
ind = 1;
for fr_span = fr_range
    [~,chi_scan_1]=dispermatrix(Ec,Ej,u(1)+0.5,g,fr_span,nlevels);
    [~,chi_scan_2]=dispermatrix(Ec,Ej,u(250)+0.5,g,fr_span,nlevels);
    chi_scan(ind) = chi_scan_1(1) - chi_scan_2(1);
    ind = ind+1;
end

FigHandleB = figure;
set(FigHandleB, 'Position', [100, 100, 800, 600]);
semilogy(u(1:250),matrixelems(1:250,2:end),'LineWidth',2,'LineStyle','-');
ylim([1E-5 2E0]);
xlabel('Offset Charge [CgVg/2e]','Interpreter','latex','FontSize',25);
ylabel('$\langle j,o|\hat{n}|0,o\rangle$','Interpreter','latex','FontSize',25);
set(gca,'TickLabelInterpreter','latex','FontSize',25);
set(gca,'YMinorTick','on');
set(gca,'XMinorTick','on');
legendCell = cellstr(num2str([1:nlevels-1]','j=%-d'));
leg = legend(legendCell,'location','best','Interpreter','latex','FontSize',25);
legend box off

FigHandleC = figure;
set(FigHandleC, 'Position', [100, 100, 800, 600]);
plot(u(1:125),chi_ip(1:125,1:2)./1e6,'LineWidth',2,'LineStyle','-');
hold on
plot(u(250:375)-u(250),chi_ip(250:375,1:2)./1e6,'LineWidth',2,'LineStyle','-');
hold off
xlim([0 0.25]);
ylim([-20 20]);
xlabel('Offset Charge [CgVg/2e]','Interpreter','latex','FontSize',25);
ylabel('$\chi_{i,p} [\frac{g^2}{200^2\,\rm{MHz}}]$ [MHz]','Interpreter','latex','FontSize',25);
set(gca,'TickLabelInterpreter','latex','FontSize',25);
set(gca,'YMinorTick','on');
set(gca,'XMinorTick','on');
leg = legend({'$|0,o\rangle$','$|1,o\rangle$','$|0,e\rangle$','$|1,e\rangle$'},'location','best','Interpreter','latex','FontSize',25);
legend box off

FigHandleD = figure;
set(FigHandleD, 'Position', [100, 100, 800, 600]);
semilogy(fr_range./1e9,abs(chi_scan./1e6),'LineWidth',2,'LineStyle','-');
xlabel('Frequency [GHz]','Interpreter','latex','FontSize',25);
ylabel('$\Delta\chi_{0} [\frac{g^2}{150^2\,\rm{MHz}}]$ [MHz]','Interpreter','latex','FontSize',25);
set(gca,'TickLabelInterpreter','latex','FontSize',25);
set(gca,'YMinorTick','on');
set(gca,'XMinorTick','on');

anharm = (EO_f(1,3)-2.*EO_f(1,2))*1E9;
chi_res = g.^2.*anharm./abs(fr-fr1*1E9)./(abs(fr-fr1*1E9)+anharm)./1E6;

disp(['fr3: ' num2str(EO_f(1,4)) ' GHz'])
disp(['Chi resonator state shift: ' num2str(chi_res) ' MHz'])
disp(['Chi parity shift: ' num2str(chi_scan(1)./1e6) ' MHz'])
disp(['Resonator FWHM: ' num2str(fr/10000/1E6) ' MHz']);
