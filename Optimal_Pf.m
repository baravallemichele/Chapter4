function [ output ] = Optimal_Pf( C0,CI,H,A,omega,gamma,lambda,distrR,distrS,VR,VS )
%Input:
%     C0,CI,H = fixed construction, marginal construction and failure costs
%     omega = obsolescence rate
%     gamma = interest rate
%     lambda = poisson process rate
%     distrS,R = S and R distriubtion type (1=Normla;2=Lognormal;15=Gumbel)
%     VR,VS = coefficients of variation
%Output:
%     output(1)=Pf_opt
%     output(2)=p_opt

ECtot=@(p) (C0+CI.*p)+((C0+CI.*p+A).*omega./gamma)+((C0+CI.*p+H).*lambda.*form_ferum( p,distrR,distrS,VR,VS )./gamma);

% Simplex minimization method
    %Optimization parameters
    p_start=12.36*VR+11.25*VS; %approximative optimal p
    options = optimset('Display','notify','TolFun',1e-7,'TolX',1e-6);
    [p_opt,fval,exitflag,output]=fminsearch(ECtot,p_start,options);
% Gradient based method
%     options = optimoptions(@fminunc,'Algorithm','quasi-newton','TolFun',1e-6,'TolX',1e-6);
%     [p_opt,fval,exitflag,output]=fminunc(ECtot,2,options);
Pf_opt=form_ferum( p_opt,distrR,distrS,VR,VS );

output=[Pf_opt,p_opt];
end

