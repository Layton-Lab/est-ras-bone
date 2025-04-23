function h = HillFunc(Conc, n, Km)
    h = (Conc.^n)./(Conc.^n + Km.^n);
end