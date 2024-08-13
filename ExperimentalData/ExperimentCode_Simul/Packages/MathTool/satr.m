function dqsat = satr(dqi)
delta = 0.018;
c = 1/delta;
    if dqi<=-delta
        dqsat = -1;
    elseif dqi>=delta
        dqsat = 1;
    else 
        dqsat = c*dqi;
    end
end