function [Pf] = form_ferum( p,distrR,distrS,VR,VS )
%form_ferum execute FORM by FERUM script
%   inputs: 
%       p=central safety factor (p=meanR/meanS)
%   output:
%       distrR,distrS= distribution function type
%       VR,VS=cov for R and S

global probdata gfundata
    for i=1:length(p)
        clear probdata femodel analysisopt gfundata randomfield systems results output_filename
            probdata.marg = [distrR, 1, VR,1-VR*1.6, nan, nan, nan, nan, 0;...
                             distrS, 1, VS,1+VS*1.6, nan, nan, nan, nan, 0];
            % Correlation matrix (square matrix with dimension equal to number of r.v.'s)
            probdata.correlation = eye(2); 
            % Determine the parameters,the mean and standard deviation associated with the distribution of each random variable
            probdata.parameter = distribution_parameter(probdata.marg);
            % define the limit state function
            gfundata(1).parameter = 'yes';
            gfundata(1).thetag = [p(i)];
            gfundata(1).expression = 'gfundata(1).thetag(1).*x(1)-x(2)'; %g=pR-S
            % Give explicit gradient expressions with respect to the involved quantities (in the order x(1), x(2), ...) if DDM is used:
            gfundata(1).dgdq = { 'gfundata(1).thetag(1)' ;
                                 '-1'};
            % Give explicit gradient expressions with respect to the limit-state function parameters 
            %(in the order thetag(1), thetag(2), ...) if DDM is used:
            gfundata(1).dgthetag = {'x(1)'};
            LoadFerumOptions;
            [formresults] = form(1,probdata,analysisopt,gfundata,femodel,randomfield);
            Pf(i)=formresults.pf1;
    end 
end

