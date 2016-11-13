function [slope slopeStd]=linearFitOnlySlope(x,y,yerr,useErrorbars,intercept,showPlots)
% Fits a straight line intercept + slope * x to the data (x,y), only  fitting the slope (the intercept is given)
% yerr is the error of the y data, and is used only if useErrorBars=1

assert(length(x)==length(y));
assert(length(x)==length(yerr));

f = @(pars,x) pars(1)*x + intercept;
parguess = y(end)-y(1); % nonlinear fits need initial guesses
w=1./yerr;

if(useErrorbars==1)
    yw=sqrt(w).*y;
    fw = @(pars,x) sqrt(w).*f(pars,x);
else
    yw=y;
    fw=f;
end

[pars,resid,J] = nlinfit(x,yw,fw,parguess);


alpha = 1-0.68;  % this is for 1-alpha confidence intervals
pars_ci = nlparci(pars,resid,'jacobian',J,'alpha',alpha);

slope=pars(1);
slopeStd    =abs(slope    -max(pars_ci(1,:)));

xfit = linspace(x(1), x(end));   yfit = f(pars,xfit);

if(showPlots==1)
     disp(['slope = ' num2str(slope) ' pm ' num2str(slopeStd) ' at confidence level ' num2str(1-alpha); ]);
     if(useErrorbars==1)
        figure(100),
            errorbare(x,y,yerr,'ro'); hold on;                   
            plot(xfit,yfit,'b-')
            xlabel('x');ylabel('y')
     else
        figure(100),
            plot(x,y,'ro','DisplayName','data'); hold on;                   
            plot(xfit,yfit,'DisplayName','fitted model')
            xlabel('x');ylabel('y')
           legend show
           legend('Location','SouthEast') 
     end

end






