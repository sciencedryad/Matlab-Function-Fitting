function plotHandle=errorbareLabel(sty,x,y,xbar,ybar,symbol)

% ERRORBARE Enhanced Errorbar Function.
%   ERRORBARE(STY,X,Y,Xbar,Ybar,symbol) 
%   It can draw errorbar along X/Y/Dual axis 
%   in normal,semilog,loglog coordinate system,
%   and adjust length of top line automatically,
%   can also control dotstyle and color in the same way with errorbar.
%
%   If the lower and upper error range of x/y is different, they should be
%   input as [lower,upper] if x/y is a column vector; 
%   for a row vector, they should be [lower;uper].
%
%   parameter STY include 12 types: 
%   v,h,d,vlogx,hlogx,dlogx,vlogy,hlogy,
%       dlogy,vlogd,hlogd,dlogd 
%   where
%   v stands for vertical errorbar,
%   h draws horizontal errorbar,
%   d means dual direction,
%   logx corresponding to semilogx,can use preffix v/h/d
%   logy corresponding to semilogy,can use preffix v/h/d
%   logd corresponding to loglog,can use preffix v/h/d
%
%   For example,
%	x = 1:10;
%	y = sin(x)+2;
%	e = std(y)*ones(size(x));
%	errorbare(x,y,e)	% use function "errorbar" directly
%	errorbare('v',x,y,e)	% "e" is error of "y" in vertical mode
%	errorbare('v',x,y,[e;2*e])  % different error limits
%	errorbare('hlogx',x,y,e)    % "e" is error of "x" in horizontal mode
%	errorbare('d',x,y,e,2*e,'or') % accept user defined linestyle
%	errorbare('dlogd',x,y,[e;e+1],[e;2*e],'ob') % all included
%
%   by Henry Sting  Email: henrysting@hotmail.com
%   $Revision: 1.2 $  $Date: 2010-3-5 $

lx=[];ux=[];ly=[];uy=[]; % initialize bottom and up limit of error bar
xl=[];xr=[];yl=[];yr=[]; % initialize short lines at the end of error bars

if ~isstr(sty)
	if nargin == 3
		errorbar(sty,x,y)
		return
	elseif nargin == 4
		errorbar(sty,x,y,xbar)
		return
	elseif nargin > 4
		error('Please assign adopted coordinate system with symbol parameters.')
	end
elseif isstr(sty)
    if size(x)~=size(y)
		error('Coordinate array should be equal.')
	end
	if nargin == 4
		symbol='ob';
		if length(x)~=length(xbar)
			error('Format of Xbar is illegal.')
		end
		if sty(1) == 'v'
			ybar=xbar;xbar=[];
		elseif sty(1) == 'h'
			ybar=[];
		elseif sty(1) == 'd'
			error('Parameters are not enough.')
		else
			error('Symbol parameter is illegal.')
		end
	elseif nargin == 5 & ~isstr(ybar)
		symbol='ob';
		if length(x)~=length(xbar)
			error('Format of Xbar is illegal.')
		elseif length(y)~=length(ybar)
			error('Format of Ybar is illegal.')
		end
	elseif nargin == 5 & isstr(ybar)
		symbol=ybar;ybar=[];
		if length(x)~=length(xbar)
			error('Format of Xbar is illegal.')
		end
		if sty(1) == 'v'
			ybar=xbar;xbar=[];
		elseif sty(1) == 'h'
			ybar=[];
		elseif sty(1) == 'd'
			error('Parameters are not enough.')
		else
			error('Symbol parameter is illegal.')
		end
	elseif nargin == 6
		if length(x)~=length(xbar)
			error('Format of Xbar is illegal.')
		elseif length(y)~=length(ybar)
			error('Format of Ybar is illegal.')
		end
		if ~isstr(symbol)
			error('Symbol should be string')
		end
	end
end


[ls,col,mark,msg] = colstyle(symbol); if ~isempty(msg), error(msg); end
symbol = [ls mark col]; % Use marker only on data part
esymbol = ['-' col]; % Make sure bars are solid


% convert input into column vector

[a,b]=size(x);
if a < b
	x=x';y=y';xbar=xbar';ybar=ybar';
	c=a;a=b;b=c;
end

%% when up-limit not equal to bottom-limit
[xa,xb]=size(xbar);
if xb==1
	ux=xbar;lx=xbar; 
elseif xb==2
	lx=xbar(:,1);ux=xbar(:,2);
end

[ya,yb]=size(ybar);
if yb==1
	uy=ybar;ly=ybar; 
elseif yb==2
	ly=ybar(:,1);uy=ybar(:,2);
end

%% plotting

dx=(max(x(:))-min(x(:)))/100;
dy=(max(y(:))-min(y(:)))/100;
logn=10;
if length(sty) == 1
	xl = x-dx; xr = x+dx; yl = y-dy; yr = y+dy; 
	plotHandle=plot(x,y,symbol);hold on
elseif length(sty) == 5 & sty(2:5) == 'logx' 
	dx=(log(max(x(:)))-log(min(x(:))))/100;
	xl = x/logn^dx;xr = x*logn^dx;yl = y-dy; yr = y+dy; 
	semilogx(x,y,symbol);hold on
elseif length(sty) == 5 & sty(2:5) == 'logy' 
	dy=(log(max(y(:)))-log(min(y(:))))/100;
	yl = y/logn^dy;yr = y*logn^dy;xl = x-dx; xr = x+dx; 
	semilogy(x,y,symbol);hold on
elseif length(sty) == 5 & sty(2:5) == 'logd' 
	dx=(log(max(x(:)))-log(min(x(:))))/100;
	dy=(log(max(y(:)))-log(min(y(:))))/100;
	xl = x/logn^dx;xr = x*logn^dx; yl = y/logn^dy;yr = y*logn^dy;
	loglog(x,y,symbol);hold on
end

%% generating vetical error
if sty(1) == 'v' || sty(1) == 'd'
vx = zeros(a*9,b);
vx(1:9:end,:) = x;
vx(2:9:end,:) = x;
vx(3:9:end,:) = NaN;
vx(4:9:end,:) = xl;
vx(5:9:end,:) = xr;
vx(6:9:end,:) = NaN;
vx(7:9:end,:) = xl;
vx(8:9:end,:) = xr;
vx(9:9:end,:) = NaN;

vy = zeros(a*9,b);
vy(1:9:end,:) = y-ly;
vy(2:9:end,:) = y+uy;
vy(3:9:end,:) = NaN;
vy(4:9:end,:) = y-ly;
vy(5:9:end,:) = y-ly;
vy(6:9:end,:) = NaN;
vy(7:9:end,:) = y+uy;
vy(8:9:end,:) = y+uy;
vy(9:9:end,:) = NaN;

plot(vx,vy,esymbol,'markersize',20)
end

%% generating horizontal error
if sty(1) == 'h' | sty(1) == 'd'
hx = zeros(a*9,b);
hx(1:9:end,:) = x-lx;
hx(2:9:end,:) = x+ux;
hx(3:9:end,:) = NaN;
hx(4:9:end,:) = x-lx;
hx(5:9:end,:) = x-lx;
hx(6:9:end,:) = NaN;
hx(7:9:end,:) = x+ux;
hx(8:9:end,:) = x+ux;
hx(9:9:end,:) = NaN;

hy = zeros(a*9,b);
hy(1:9:end,:) = y;
hy(2:9:end,:) = y;
hy(3:9:end,:) = NaN;
hy(4:9:end,:) = yl;
hy(5:9:end,:) = yr;
hy(6:9:end,:) = NaN;
hy(7:9:end,:) = yl;
hy(8:9:end,:) = yr;
hy(9:9:end,:) = NaN;

% p=plot(hx,hy,esymbol,'markersize',20);


end
