function [ output_args ] = EC_tot( mean_beta )
%EC_tot evaluates the total expected costs as a function of the mean beta
%for beta normal, lognormal or weibull distributed with given cov 

%Find p corresponding to mean_beta
   p=find_p_from_beta( mean_beta ); 

global data
    for i=1:length(mean_beta)
        if data.cov_beta==0 %Case for beta known (COV beta == 0)
            ECtot=( (data.C0+data.CI.*p(i)) + (data.C0+data.CI.*p(i)+data.AA).*(data.omega/data.gamma) + (data.C0+data.CI.*p(i)+data.H ).*normcdf(-mean_beta(i)).*data.lambda./data.gamma );
        else
            %Define the integration domain boundaries
                int_low_lim=mean_beta(i)*(1-6*data.cov_beta);
                int_up_lim=mean_beta(i)*(1+6*data.cov_beta);
            switch data.distr_type_beta
                case 1 %Normal distributed beta
                    ECtot=integral(@(b) normpdf(b,mean_beta(i),mean_beta(i).*data.cov_beta).*( (data.C0+data.CI.*find_p_from_beta(b)) + (data.C0+data.CI.*find_p_from_beta(b)+data.AA).*(data.omega/data.gamma) + (data.C0+data.CI.*find_p_from_beta(b)+data.H ).*normcdf(-b).*data.lambda./data.gamma ),int_low_lim,int_up_lim);
                case 2 %Lognormal distributed beta
                    std_beta(i)=mean_beta(i).*data.cov_beta;
                    %Find the parameters of the Lognormal dsitribution
                        mu_beta(i)=log(mean_beta(i)^2/sqrt(mean_beta(i)^2+ std_beta(i)^2));
                        sigma_beta(i)=sqrt(log(data.cov_beta^2+1));
                    ECtot=integral(@(b) lognpdf(b,mu_beta(i),sigma_beta(i)).*( (data.C0+data.CI.*find_p_from_beta(b)) + (data.C0+data.CI.*find_p_from_beta(b)+data.AA).*(data.omega/data.gamma) + (data.C0+data.CI.*find_p_from_beta(b)+data.H ).*normcdf(-b).*data.lambda./data.gamma ),max(0,int_low_lim),int_up_lim);
                case 3 %Weibull
                    %Approximate the Weibull parameters
                        syms alfa bet
                        bet_start=max(0.1,(data.cov_beta)^(-1.086)); %approx starting value
                        alfa_start=max(0.1,mean_beta(i)/gamma(1+1/bet_start));
                    %Exact Weibull parameters
                        S=vpasolve([mean_beta(i)==alfa*gamma(1+1/bet),mean_beta(i).*data.cov_beta==alfa*sqrt((gamma(1+2/bet))-(gamma(1+1/bet))^2)],[alfa,bet],[alfa_start,bet_start]);
                    ECtot=integral(@(b) wblpdf(b,eval(S.alfa),eval(S.bet)).*( (data.C0+data.CI.*find_p_from_beta(b)) + (data.C0+data.CI.*find_p_from_beta(b)+data.AA).*(data.omega/data.gamma) + (data.C0+data.CI.*find_p_from_beta(b)+data.H ).*normcdf(-b).*data.lambda./data.gamma ),max(0,int_low_lim),int_up_lim);
            end %End switch
        end %End if
        if mean_beta(i)<1 
            ECtot=1e+999
        end
        output_args(i)=ECtot;
    end %End cycle on i
end %End function