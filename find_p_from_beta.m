function [ p ] = find_p_from_beta( beta )
% Given a beta value and the other input included in data. the function
% finds the value of p corresponding to this beta using the polynomial
% approxiamting the relation beta(p)

global data

    for i=1:length(beta)
        p(i)=polyval(data.coeff,beta(i));
    end
end

