function [EE,EO,DE] = solvesystem(Ec,Ej,u, deltaL, deltaR)

    nlevels = 4; % Only care about first n levels
    EE = zeros(length(u),nlevels); % Even
    EO = EE; % Odd
    for i=1:length(u)
        [~,eiva]=eigensystem(Ec,Ej,u(i)+1E-4);
        EE(i,:) = diag(eiva(1:nlevels,1:nlevels));

        [~,eiva]=eigensystem(Ec,Ej,u(i)+0.5+1E-4);
        EO(i,:) = diag(eiva(1:nlevels,1:nlevels));
    end
    
    DE=(EO(:,1)-EE(:,1))+deltaL-deltaR; %difference between even and odd state energies