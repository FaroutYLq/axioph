function [melem, chi]=dispermatrix(ec,ej,u,g,wr,nlevels)

planck = 4.135E-15; % eV.s

% same function as 
n=30; % cutoff
H = zeros(2*n+1,2*n+1);  
for l=1:(2*n+1)
    H(l,l) = 4.*ec.*((l-n-1)-u).^2;
    if (l+1 <= (2*n+1))
        H(l,l+1)=-ej/2;
        H(l+1,l)=-ej/2;
    end
end

[Eivec,Eiva]=eig(H);
[U,K]=sort(diag(Eiva));
V = zeros(2*n+1,n); % eigenvectors
for j=1:2*n+1
    V(:,j)=Eivec(:,K(j)); % sorted from smallest to largest
end
Eivec=V; % column vectors
Eiva=U; 

% Need to solve for matrix elements of dispersion system
% Serniak notation <j,p|n|i,p>, i!=j. Let's compute j=0..nlevels-1
% what is the number operator??
% Serniak uses n ~ sum_j[(j-ng)|j><j|] 
% while shouldn't it be n=i*(Ej/32Ec)^(1/4)*(bdagger-b)?
    % with b equiv sum_j[sqrt(j+1)|j><j+1|] 

prefactor = 1;
nop = zeros(2*n+1);
for l=1:2*n+1
    nop(l,l)=prefactor.*((l-n-1)-u);
end

% now compute dispersive shift
% chi_i = g^2 sum_j!=1[2*wij*|<j,p|n|i,p>|^2 / (wij^2 - wr^2)
chi = zeros(1,nlevels);
melem = zeros(nlevels,nlevels);
for i=1:nlevels
    for j=1:2*n+1
        if (i~=j)
            wij = (Eiva(i)-Eiva(j))/planck; % Hz
            mtemp = abs(Eivec(:,j)'*nop*Eivec(:,i)).^2 ;
            if (j<=nlevels)
                melem(i,j) = mtemp;
            end
            chi_temp = 2.*wij.*mtemp./(wij.^2-wr.^2);
            chi(i)=chi(i)+chi_temp;
        end
    end
end
chi = g.^2 .* chi;