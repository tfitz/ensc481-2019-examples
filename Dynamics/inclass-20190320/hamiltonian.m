function H = hamiltonian(Z, m, g, l, M, J_G)

ntime = size(Z,2);
H = zeros(ntime,1);

for i = 1:ntime
    z = Z(:,i);
    
    H(i) = (J_G.*M.*g.*l.*m.*cos(z(2)) + J_G.*g.*l.*m.^2.*cos(z(2)) + J_G.*z(3).^2/2 + M.*g.*l.^3.*m.^2.*cos(z(2)) + M.*z(4).^2/2 - g.*l.^3.*m.^3.*cos(z(2)).^3 + g.*l.^3.*m.^3.*cos(z(2)) + l.^2.*m.*z(3).^2/2 - l.*m.*z(3).*z(4).*cos(z(2)) + m.*z(4).^2/2)./(J_G.*M + J_G.*m + M.*l.^2.*m + l.^2.*m.^2.*sin(z(2)).^2);
    
end
