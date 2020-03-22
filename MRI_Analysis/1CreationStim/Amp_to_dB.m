function [amplitude] = Amp_to_dB(decibels)
%Amplitud to decibels conversion: Calculates the equivalence of amplitud to
%decibels through the formula
%Power ratio = x
%Amplitude = sqrt(x)
%Decibels = 10*log10(x)
%   Detailed explanation goes here
amplitude = sqrt(10^(decibels/10));
end