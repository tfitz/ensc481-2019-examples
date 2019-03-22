function dz = cart(t,z, m, g, l, M, J_G, f, torque)


dz = [(J_G.*z(3) + l.^2.*m.*z(3) - l.*m.*z(4).*cos(z(2)))./(J_G.*M + J_G.*m + M.*l.^2.*m + l.^2.*m.^2.*sin(z(2)).^2); (M.*z(4) - l.*m.*z(3).*cos(z(2)) + m.*z(4))./(J_G.*M + J_G.*m + M.*l.^2.*m + l.^2.*m.^2.*sin(z(2)).^2); f(t,z); 2*l.^2.*m.^2.*(J_G.*M.*g.*l.*m.*cos(z(2)) + J_G.*g.*l.*m.^2.*cos(z(2)) + J_G.*z(3).^2/2 + M.*g.*l.^3.*m.^2.*cos(z(2)) + M.*z(4).^2/2 - g.*l.^3.*m.^3.*cos(z(2)).^3 + g.*l.^3.*m.^3.*cos(z(2)) + l.^2.*m.*z(3).^2/2 - l.*m.*z(3).*z(4).*cos(z(2)) + m.*z(4).^2/2).*sin(z(2)).*cos(z(2))./(J_G.*M + J_G.*m + M.*l.^2.*m + l.^2.*m.^2.*sin(z(2)).^2).^2 + torque(t,z) - (-J_G.*M.*g.*l.*m.*sin(z(2)) - J_G.*g.*l.*m.^2.*sin(z(2)) - M.*g.*l.^3.*m.^2.*sin(z(2)) + 3*g.*l.^3.*m.^3.*sin(z(2)).*cos(z(2)).^2 - g.*l.^3.*m.^3.*sin(z(2)) + l.*m.*z(3).*z(4).*sin(z(2)))./(J_G.*M + J_G.*m + M.*l.^2.*m + l.^2.*m.^2.*sin(z(2)).^2)];