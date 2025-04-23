function pct_ACEi = get_pct_ACEi(do_ACEi, t, pct_ACEi_inhib, age_ACEi)
   
if do_ACEi
    t_yrs = t/(365*24);
    if t_yrs > age_ACEi
        if t_yrs < age_ACEi + 1
            % gradual decrease for numerics
            pct_ACEi = (t_yrs - age_ACEi)*pct_ACEi_inhib;
        else
            pct_ACEi = pct_ACEi_inhib;
        end
    else
        pct_ACEi = 0;
    end
else
    pct_ACEi = 0;
end

end