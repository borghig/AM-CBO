function [igd] = compute_IGD(exact, computed)

Ne = length(exact);

igd = 0;
for i = 1:Ne
    igd = igd + min(vecnorm(exact(i,:) - computed,2,2))/Ne;
end

    
end