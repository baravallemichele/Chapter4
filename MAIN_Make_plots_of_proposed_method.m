%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Scripts behind calculations presented in the article:
% A framework for estimating the implicit safety level of existing design codes
% Proceedings of the 12th International Conference on Structural Safety & Reliability (ICOSSAR2017) 6-10 August, 2017 Vienna, Austria.
%     Michele Baravalle a,*
%     Jochen Köhler a
% a)Dept. of Structural Engineering, Norwegian University of Science & Technology, Richard Birkelands vei 1A, 7491 Trondheim, Norway
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
clc
clear global
%% INPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%define global variables
global data
%Random variables
    %Resistance
    data.VR=0.15;   %Coefficient of variation   
    data.distrR=2;  %Distribution type (1 = Normal;2 = Lognormal;3 = Gamma;4 = Shifted Exponential marginal distribution;  5 = Shifted Rayleigh marginal distribution; 6 = Uniform distribution; 7 = Beta; 8 = Chi-square; 11 = Type I Largest Value marginal distribution;12 = Type I Smallest Value marginal distribution                                                             %
                                      %13 = Type II Largest Value marginal distribution;14 = Type III Smallest Value marginal distribution ;15 = Gumbel (same as type I largest value) ;16 = Weibull marginal distribution (same as Type III Smallest Value marginal distribution with epsilon = 0 ) %
    %Load
    data.VS=0.30;   %Coefficient of variation 
    data.distrS=15; %Distribution type
%Costs
    data.C0=1;%     %Fixed constr costs -- OBS:DO NOT CHANGE THIS VALUE --
    data.CI=.01;    %Slope of the lienarized construction costs (C(p)=C0+CI*p)
    data.H =3.2;    %Costs associated to failure
    data.AA=0;      %Costs of demolition
%Parameters in the objective function    
    data.gamma=0.03;%Societal interest rate
    data.omega=0.02;%Obsolescence rate
    data.lambda=1;  %Parameter of the Poisson process representing failure events
    data.distr_type_beta=2; %Type of distribution representing the variation of the code relaibility index -- ONLY Normal, Lognormal and Weibull implemented so far --
%Options for the fminsearch function    
    options = optimset('Display','notify','TolFun',1e-6,'TolX',1e-6);
%% SCRIPTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Approximate the relation between p and beta with a 5th order polynomial
    figure          %figure created for visual inspection of the goodness of the approxiamtion
    dum_p=0.1:0.1:20; %values of p where fitting the polynomial
    dum_beta=-norminv(form_ferum( dum_p,data.distrR,data.distrS,data.VR,data.VS ));
    dum_p=dum_p(dum_beta<Inf);dum_beta=dum_beta(dum_beta<Inf);
    data.coeff = polyfit(dum_beta,dum_p,5);
    plot(dum_beta,dum_p,'.'); hold on
    plot(dum_beta,polyval(data.coeff,dum_beta));
    xlabel('p'); ylabel('\beta(p)');legend('\beta(p) exact','\beta(p) approximated')
    title('Check \beta(p) polynomial approximation')
%% Create plots for COV_beta and CI/H
    figure
    H = 0:0.1:10;   %Define H values 
    Vbeta = [0 0.02:0.02:0.3]; %Defined values of beta coeff. of variation
    [H,Vbeta] = meshgrid(H,Vbeta); 
    CI=[.1 .01 .001]; %Define three CI values for making the three subplots 
    num_cycles=3*size(H,1)*size(H,2); %Number of cycles
    dum=0; %Initialize dummy variable
    for j=1:3
        for i_h=1:size(H,2)
            for i_vbeta=1:size(Vbeta,1)
                dum=dum+1;
                data.CI=CI(j);
                data.H =H(i_vbeta,i_h);
                data.cov_beta=Vbeta(i_vbeta,i_h);
                %Find mean beta minimizing the total Expected costs
                    mean_beta_start_value=3;
                    [mean_beta_opt(i_vbeta,i_h),ECtot_opt(i_vbeta,i_h),exitflag,output]=fminsearch(@EC_tot,mean_beta_start_value,options); 
            end %End cycle on i_vbeta
            disp(strcat(num2str(round(dum/num_cycles*100)),'%')) %percentage of points calculated
        end %End cycle on i_h

        subplot(1,3,j) %Subplot for the j-th calue of CI
        hold on
        set(subplot(1,3,j), 'Position', [0.1+(0.3*(j-1)) 0.17 0.27 0.67 ]);
        C=contour(Vbeta,H,mean_beta_opt,[2.5:0.1:6],'Color',[.8 .8 .8],'LineStyle','-');
        C=contour(Vbeta,H,mean_beta_opt,[2.5:0.5:6],'LineColor','k'); hold on
        clabel(C,'FontSize',15,'Color','black')
        if j==2
            title({'Mean reliability index (\mu_{\beta} ) at optimum';strcat('C_I / C_0 = ',num2str(data.CI/data.C0))})
        else
            title(strcat('C_I / C_0 = ',num2str(data.CI/data.C0)))
        end
        xlabel('COV_{\beta}'); set(gca,'YTick',0:10)
        if j==1
            ylabel('H / C_0');grid on;
        else
            set(gca,'YTickLabel',{})
        end
        axis([0 max(Vbeta(:,1)) 0 max(H(1,:))])
        grid on
        grid minor
        box on
    end %End cycle on j






