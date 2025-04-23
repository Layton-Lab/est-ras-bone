function pct_ARB = get_pct_ARB(do_ARB, t, pct_ARB_inhib, age_ARB)
    t_yrs = t/(365*24); % time in years
    if do_ARB
        if t_yrs > age_ARB
            if t_yrs < age_ARB + 1
                % gradual decrease for numerics
                pct_ARB = (t_yrs - age_ARB)*pct_ARB_inhib;
            else
                pct_ARB = pct_ARB_inhib;
            end
        else
            pct_ARB = 0;
        end
    else
        pct_ARB = 0;
    end
end