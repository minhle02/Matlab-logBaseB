function result = loga(n,A,digit)
%result = loga(n,A,digit) take logarithm base n of each element in matrix A
%This function can take log base n of negative number, but base n must 
%greater than 0
%Input:
%   n: the base of logarithm
%   A: input matrix/array/digit
%   digit: number of significant figures. Default = 3
%Output:
%   result: matrix A with each element is taken logarithm of base n

if nargin<2
    error('Minimum argument is 2');
end
if nargin<3
    digit = 3;
end
if n<0
    error('the base cannot be negative');
end
if digit<1
    error('Minimum number of digits is 1');
end

err = (0.5*10^(2-digit))/100;
a = size(A,1);
b = size(A,2);
result = zeros(a,b);

for i=1:a
    for j=1:b
        y = A(i,j);
        result(i,j) = iloga(n,y,err);
    end
end
end


function result = iloga(n,y,error)
%Take logarithm of a single number
flag = 0;
if n == 1
    result = 0;
    return;
else
    if y>1
        a = n;
        n = y;
        y = a;
        flag = 1;
    end
    if abs(y)<1
        flag = 2;
        y = 1/y;
    end
end


if y>0
    result = posloga(n,y,error);
elseif y==0
    result = Inf;
else
    result = negloga(n,y,error);
end

if flag == 1
    result = 1/result;
elseif flag == 2
    result = -result;
end
end

function result = posloga(n,y,error)
%Take logarithm base n of positive number only
sign = 1;
if y<1
    y = 1/y;
    sign = -1;
end

result = whole(n,y);
remain = y/n^result;
if remain > 1
    func = @(x) remain - n^x;
    root = bisect(func,0,1,error);
    result = result+root;
end
result = result*sign;
end

function result = negloga(n,y,error)
%Take logarithm base n of negative number only
result = posloga(n,-y); %real part of result

%logn(i) = ln(i)/ln(n)
a = posloga(exp(1),n,error);%calculate ln(n)

%Euler identity: e^(i*pi) = -1=i^2, thus ln(i) = (pi/2)*i
%When y<0, logn(y) = logn(-y*-1) = logn(-y*i^2) = logn(-y)+2logn(i)
result = result + (1/a)*pi*1i;
end


function w = whole(n,y)
%Calculate the whole part of final result
w = 0;
while y>=n
    w = w+1;
    y = y/n;
end
end

function [root,fx,ea,iter]=bisect(func,xl,xu,es,maxit,varargin) % bisect: root location zeroes
% [root,fx,ea,iter]=bisect(func,xl,xu,es,maxit,p1,p2,...):
% uses bisection method to find the root of func
% input:
%   func = name of function
%   xl, xu = lower and upper guesses
%   es = desired relative error (default = 0.0001%)
%   maxit = maximum allowable iterations (default = 50)
%   p1,p2,... = additional parameters used by func
% output:
% root = real root
% fx = function value at root
% ea = approximate relative error (%)
% iter = number of iterations
if nargin<3,error('at least 3 input arguments required'),end
test = func(xl,varargin{:})*func(xu,varargin{:});
if test>0,error('no sign change'),end
if nargin<4||isempty(es), es=0.0001;end
if nargin<5||isempty(maxit), maxit=50;end
iter = 0; xr = xl; ea = 100;
while (1)
    xrold = xr;
    xr = (xl + xu)/2;
    iter = iter + 1;
    if xr ~= 0,ea = abs((xr - xrold)/xr) * 100;end 
    test = func(xl,varargin{:})*func(xr,varargin{:}); 
    if test < 0
        xu = xr;
    elseif test > 0
        xl = xr; 
    else
        ea = 0; 
    end
  if ea <= es || iter >= maxit,break,end
end
root = xr; fx = func(xr, varargin{:});
end