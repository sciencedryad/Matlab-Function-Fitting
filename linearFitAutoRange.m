function [slope sigSlope intercept sigIntercept r2 xCutoffFound xCutoffUsed]...
            =linearFitAutoRange(x,y,sigy,r2Limit, xCutoffLarge,useAutoDeterminedFitRange,useErrorbars)      
                
% Linear fit to the data (x,y) with errors sigy (only used if useErrorbars=1).
% If useAutoDeterminedFitRange=1, tries to determine a fit region
%       x(1:xCutoffUsed)  for which r2>=r2Limit.
%       Method: fits in the region x(1:xCutoffTest), 
%       varying xCutoffTest from xCutoffLarge:-1:3.
%       If a satisfying fit with (r2>r2Limit) is found (for the first time), the program returns
%            xCutoffFound=1 and xCutoffUsed=xCutoffTest.
%       If no satisfying fit is found (r2 is always too small), then
%           xCutoffFound=0, and the data is fitted for x(1:xCutoffLarge).
% If useAutoDeterminedFitRange=0, 
%    then xCutoffFound=0, and the data is fitted for x(1:xCutoffLarge).

xCutoffFound=0;    
        
%=======================================
% use this for automatically determined fit line
%=======================================
if(useAutoDeterminedFitRange==1)    
       indexMax=find(x<xCutoffLarge,1,'last');
       xValuesToTest=x(indexMax:-1:5);
       %disp(['testing ' num2str(length(xValuesToTest)) ' cutoffs' ])
       for  xCutoffTest=xValuesToTest;
           %disp(['       testing xCutoffTest = ', num2str(xCutoffTest)]);
           xSelect=logical((x<=xCutoffTest));
           xToFit=x(xSelect);   yToFit=y(xSelect);  sigyToFit=sigy(xSelect);
           [intercept slope sigIntercept sigSlope r2 q phiFit delta]  =linearFit(xToFit,yToFit,sigyToFit,useErrorbars,0) ;
            if(r2>r2Limit) 
                %disp('FOUND A GOOD FIT!')
                xCutoffUsed=xCutoffTest;
                xCutoffFound=1;
                break; 
            end;
       end        
end
%=======================================
%use this fpr pre-determined cutoff, or if no reasonable cutoff could be found
%=======================================
if(xCutoffFound==0) % this is also true for useAutoDeterminedFitRange!=1 since then xCutoffFound is never changed from 0
           %disp('DID NOT FIND A GOOD FIT!')  
           xCutoffUsed=xCutoffLarge;
           xSelect=logical((x<xCutoffUsed));
           xToFit=x(xSelect);   yToFit=y(xSelect);  sigy=sigy(xSelect);
           [intercept slope sigIntercept sigSlope r2 q phiFit delta]=linearFit(xToFit,yToFit,sigy,useErrorbars,0) ;           
end

