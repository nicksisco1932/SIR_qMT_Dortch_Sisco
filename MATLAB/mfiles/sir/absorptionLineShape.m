function g = absorptionLineShape(T2,delta,lineShape)
% T2: T2 of pool (sec)
% delta: freq offset of RF pulse (Hz)
% lineShape: 'super-lorentzian','lorentzian','gaussian'
% g(2*pi*delta): absorption lineshape

% Make sure delta is a row vector
if size(delta,1) ~= 1
    delta = delta';
end

% Calculate g for specified lineshape
switch lineShape
    case 'super-lorentzian'
        du = 5e-4; u = 0:du:1;
        f = ones(length(delta),1)*(sqrt(2/pi)*(T2./abs(3*u.^2-1))).*...
            exp(-2*((2*pi*delta*T2)'*(1./abs(3*u.^2-1))).^2);
        g = sum(f,2)*du;
    case 'lorentzian'
        g = (T2/pi)*(1./(1+(2*pi*delta*T2).^2));
    case 'gaussian'
        g = (T2/sqrt(2*pi)).*exp(-(2*pi*delta*T2).^2/2);
end
