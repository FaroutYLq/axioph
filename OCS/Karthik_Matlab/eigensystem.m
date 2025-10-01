% Reproduced from P. Echternach 
function [Eivec,Eiva]=eigensystem(ec,ej,u)
n=18; % cutoff
H = zeros(2*n+1,2*n+1);  
for l=1:(2*n+1)
    H(l,l) = 4.*ec.*((l-n-1)-u).^2;
    if (l+1 <= (2*n+1))
        H(l,l+1)=-ej./2;
        H(l+1,l)=-ej./2;
    end
end

[Eivec,Eiva]=eig(H);
[U,K]=sort(diag(Eiva));
V = zeros(2*n+1); % eigenvectors
for j=1:2*n+1
    V(:,j)=Eivec(:,K(j));
end
Eivec=V;
Eiva=diag(U);