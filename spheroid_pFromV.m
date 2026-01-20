function out = spheroid_pFromV (m,nu,mu, type)

% compute the equivalent pressure for given DeltaV: all equations follow
% from spheroid.m. OR vice versa;  if type is pressure then compute dV if
% type is volume then compute P

a = m(1);
b = m(2);
strength = m(8);

c = sqrt(a^2 - b^2);
b2 = b^2;
C0 = 4*(1 - nu);
%C1 = C0*(1 - 2*nu);
%C2 = C0 - 1;
C3 = C0 + 1;
%C4 = C0 - 3;

    L0 = log((a - c)/(a + c));
    d = 1/(b2*c*(2*b2 - C0*a^2)*L0 - a*c^2*(4*b2 + C0*a^2) + a*b^4*(1 + nu)*L0^2);
    a1 = d*b2*(2*c*(a^2 + 2*b2) + 3*a*b2*L0);
    b1 = d*(a^2*c*C0 + b2*(c + C3*(c + a*L0)));
    
if strcmpi(type,'volume')
    % P in units of mu
    out = 3*strength*mu/(2*pi*b2)/( a1*(a*L0*(-1 + 2*nu) + c*(-5 + 4*nu)) - 2*b1*c^3);
else
    % dV in km^3
    out = strength*(  3*mu/(2*pi*b2)/( a1*(a*L0*(-1 + 2*nu) + c*(-5 + 4*nu)) - 2*b1*c^3))^(-1);
end
    out = real(out);
    if out == 0;
        keyboard
    end

return

% C5 = strength/(4*mu);
% Note     u = C5*(u1 - u2); multiplies displacement%     if nargin > 4
%         if strcmpi(flag,'volume')
%             C5 = strength*3/(8*b2*pi*(-2*b1*c^3 - a1*(a*L0*(1 - 2*nu) + c*C3)));
%         end
%     end

%to invert this



DelV = DelV*1e9  % to get m^3 from km^3




