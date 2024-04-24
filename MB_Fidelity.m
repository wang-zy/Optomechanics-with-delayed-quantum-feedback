function [F, phi] = MB_Fidelity(r, N)
% calculate the multi bounce fidelity
% parameters: 4g0/kappa, N

alpha = 0;
phi = 0;

for i = 1:4*N
    D = Displace(r, alpha);
    phi = phi + imag(D * alpha');
    alpha = alpha + D;
    alpha = - 1j * alpha;
end

F = exp(-abs(alpha)^2 / 2);
phi = phi;%mod(phi, pi);

end

function p = Displace(r, alpha)
% return the state dependent displacement operator

p = -1j * r / (1 + r^2 * real(alpha)^2);

end