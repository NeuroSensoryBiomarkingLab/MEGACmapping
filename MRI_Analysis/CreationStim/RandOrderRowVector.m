function [RandVector] = RandOrderRowVector(max,rep)
%Growing Integers Row Vector
%Creates a row vector of growing numbers with X number of repetitions for
%each number until the maximum number is reached, and then, creates a list
%of them arranged in a random way.
%MAX is the number of frequencies and REP is the number of repetitions each
%one of them should have.

vector = [];
for ii=0:(max)
    vector(((ii*rep)+1):((ii*rep)+rep))=ii;
end

RandVector=vector(randperm(length(vector)));

end

