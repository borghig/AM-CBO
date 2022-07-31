%Computes the pareto points for a given sample
%Assuming fi(X), i = 1...n, are objective functions
% for a minimization problem.
% A point X* is said to be Pareto optimal
% if there is no X such that fi(X) <= fi(X*) for
% all i = 1...n, with at least one strict inequality.
%
% [dominants, dompos, ndompos] = faster_pareto2(sample),
% sample - m x n input matrix:
% sample = [f1(X1) f2(X1) ... fn(X1);
%           f1(X2) f2(X2) ... fn(X2);
%           .......................
%           f1(Xm) f2(Xm) ... fn(Xm)]
% dominants - an output matrix with rows which are Pareto points of sample.
% [dominants, dompos, ndompos] = fast_pareto(sample) dompos is a row vector
% containing the position of the dominant solution and ndompos is a row
% vector containing the remaining ones
% in the sample.
% Example:
% input
%sample = [14    16     8    12;
%          15    16    15    13;
%          3     4     13    13;
%          15    16    16    7 ;
%          11    16    11    11;
%          3     9     2     4 ;
%          6     14    14    12;
%          10    4     16    2]
% outputs
% dominants = [3     4    13    13;
%              3     9     2    4;
%              10    4    16    2]
%  dompos = [3    6     8]
%  ndompos = [1     2     4     5     7]
% the positions are related to the sample
% by Ahmed  HASSAN (ahah432@yahoo.com) and
%Maurel AZA-GNANDJI(my.lrichy@gmail.com)
%30-07-2015
function [dominants, dompos, ndompos] = faster_pareto2(sample)
dominants = sample;
dompos = 1:size(dominants,1);
ndompos = [];
% Rule : no other row can be better in all elements
i=1;
j=1;
while i <= size(dominants,1)
    r = size(dominants,1);
    bb = ones(r,1)*dominants(i,:)-dominants;
    bb(i,:) = '';
    if any(all(bb'>=0))
        dominants(i,:)=[];
        ndompos(j) = dompos(i);
        dompos(i) = '';      
        
        j=j+1;
        i=i-1;   % continue
    end
    
    
    i = i+1;
end
if isempty(dompos)
    warning('No pareto points found, The result is an empty matrix.')
end
