function [RealRandomNumbers] = RealRand(interval,r,c)
%Real Random Numbers: 
%The numbers will never be repeated each time the function is run.
%Random numbers will be between -interval and interval, and the matrix will
%be the size of [r,c]
if nargin < 3
    c = r;
end

rng('shuffle');
RealRandomNumbers = -interval+(interval+interval)*rand(r,c);
end