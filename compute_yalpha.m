%auxiliary functions

function [yalpha] = compute_yalpha(QoI,E, alpha)

% weigthed average  of QoI according to Gibbs distribution associated to E
% alpha = inverse of temperature param
% E = function to use for Gibbs distribution

weights = exp(-alpha.*(E - min(E)));
yalpha = weights'*QoI/sum(weights);
yalpha = reshape(yalpha,size(QoI(1,:)));

end