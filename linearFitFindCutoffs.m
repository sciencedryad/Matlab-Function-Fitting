function [a b siga sigb r2 q  xLowLimit]=linearFitFindCutoffs(x,y,output,r2Limit)
% Fits a straight line a+b*x to the data (x,y), with y errors sigy.
% Searches for the smallest xLimit so that the fit within [xLimit,x(end)] has r2>r2Limit.
% If no such xLimit is found, then xLimit=-1;
% More details on the fit and the other parameters: see linearFit.m
assert(length(x)==length(y));

xLowLimitTestVec=linspace(x(1),x(end-3),2*length(x));
useErrorBars=0;
xLowLimit=-1;

for xLowLimitTest=xLowLimitTestVec
   % disp(['xLowLimitTest=' num2str(xLowLimitTest) ]);
    xSelect=logical((x>xLowLimitTest));
    xs=x(xSelect);  ys=y(xSelect);  dy=ys; 
    %dof =length(xs)-2;
    [a b siga sigb r2 q y_fit delta]=linearFit(xs,ys,dy,useErrorBars,0);
    %r2=1-(1-r2)/dof;         
    if(r2>r2Limit) 
          xLowLimit=xLowLimitTest;
         [a b siga sigb r2 q y_fit delta]=linearFit(xs,ys,dy,useErrorBars,1);
          break; 
    end;
end


end




