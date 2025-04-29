function u = Initial4(epsilon, k0, r, x)
    s=1;
    for k = 1:64
        psik = r(k);
        Ek = (k/k0)^4*exp(-2*(k/k0)^2);
        s = s + epsilon*sqrt(Ek)*sin(2*pi*k*(x+psik));
    end
    u = s;

end