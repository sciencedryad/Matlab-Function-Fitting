function [a b siga sigb r2 q y_fit delta]=linearFit(x,y,dy,useErrorBars,output)
% Fits a straight line a+b*x to the data (x,y), with y errors sigy.
%   If there are no error bars, set useErrorBars=0, otherwise useErrorBars=1
% Returns the fit parameters a, b and their respective errors siga, sigb.
% Also returns the goodness-of-fit-measure q, which should be as close to
%   1 as possible. Basically, q is the probability that a value of chi2 as
%   poor as the obtained one should occur by chance. 
% Also returns the vector y_fit, which are the points of the fit, 
%   and the vector delta. y +/- delta is at least a 50% prediction interval.  
%   For large degrees of freedom, the confidence level approaches approximately 68%.

% More details on q: 
% Null-hypothesis: Deviations from the linear fit are randomly distributed (Gaussian)
% q = probability to observe the given data, assuming that the null hypothesis is true
%   = 'probability that the null hypothesis is true'
% If q is very small, say q<0.1, then the null hypothesis can be rejected,
%      and the linear fit is probably not very good
% If q is large, the data are at least consistent with the null hypothesis.
% If no error bars are available (or used), then q=1, but this is not meaningful.

if(useErrorBars==0)
    sigy=ones(1,length(dy));
else
    sigy=dy;
end

assert(length(x)==length(y));
assert(length(x)==length(dy));
%ndata=length(y);

%============================
% use weighted polynomial fit
%============================
    weights=1./sqrt(sigy);
    [p s] = weightedpolyfit(x, y, weights,1);
        % p contains the coefficients of the polynomial, 
        % Structure s is for use with polyval to get error estimates on the coefficients.
        %   It contains fields R, df, and normr:
        %   R = triangular factor from a QR decomposition of the Vandermonde matrix of x, 
        %   df = the degrees of freedom, 
        %   normr = the norm of the residuals. 
        % If the data y are random, an estimate of the covariance matrix of p 
        %   is (Rinv*Rinv')*normr^2/df, where Rinv is the inverse of R. 
        %   If the errors in the data y are independent normal with constant variance, 
        %   polyval produces error bounds that contain at least 50% of the predictions.
    [y_fit delta] = polyval(p, x, s);
        % y_fit gives the y points of the fit 
        % delta gives errors in y data: for any
        %   degrees of freedom, Y +/- DELTA is at least a 50% prediction interval
        %   in all cases.  For large degrees of freedom, the confidence level
        %   approaches approximately 68%.

    a=p(2);
    b=p(1);
    

%============================
% calculate errors on fit parameters
%============================
    chi2=sum(((y-a-b*x)./sigy).^2);
    Rinv=inv(s.R);
    covMat=(Rinv*Rinv')*s.normr^2/s.df;
    siga=sqrt(covMat(2,2));
    sigb=sqrt(covMat(1,1));
%     if(useErrorBars==0) % no available errors
%         sigfac=sqrt(chi2/s.df);
%         siga=siga*sigfac;
%         sigb=sigb*sigfac;
%     end

%============================
% assess goodness of fit    
%============================
    if(useErrorBars==0) % no error bars
        q=1.0;
    else % with error bars
        q=1-gammainc( 0.5*chi2, 0.5*s.df);
    end
    SSerr=sum((y-y_fit).^2); % residual sum of squares
    SStot=sum((y-mean(y)).^2); % residual sum of squares
    r2=1-SSerr/SStot;

%============================
% output if desired
%============================
if(output==1)
%     if(useErrorBars==0)
%        plot(x,y,'x');
%     else
%         errorbar(x,y,sigy,'x');    
%     end
    hold on
    plot(x,polyval(p,x),'k-')
    %plot(x, y_fit + 2* delta, ':', x, y_fit - 2* delta, ':')
    %plot(x, y_fit + 2* delta, ':', x, y_fit - 2* delta, ':')
%     textsetoff=50;
%     text(min(x)+10,max(y)-1,['slope = ',num2str(b),' pm ',num2str(sigb)]);
%     text(min(x)+10,max(y)-textsetoff,['intercept = ',num2str(a),' pm ',num2str(siga)]);
%     if(useErrorBars==0)
%         text(min(x)+10,max(y)-2*textsetoff,['r2 = ',num2str(r2),...
%             ', fit without errors']);
%     else
%         text(min(x)+10,max(y)-2*textsetoff,['r2 = ',num2str(r2),', q=',num2str(q)]);
%     end

    if(useErrorBars==0)
        disp(['  Fit result: b=',num2str(b),' pm ',num2str(sigb),', a=',num2str(a),' pm ',num2str(siga),...
          ', r2 = ',num2str(r2),]);%', (no estimate q of goodness of fit without errors)']);
    else
        disp(['  Fit result: b=',num2str(b),' pm ',num2str(sigb),', a=',num2str(a),' pm ',num2str(siga),...
           ', r2 = ',num2str(r2),', q=',num2str(q)]);
    end
end




