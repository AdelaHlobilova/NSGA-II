%% This function assigns ranks to all solutions in F in NSGA-II.
% First of all, all solutions are sorted to fronts according to objective
% functions. Not all solutions have to be sorted to save computational
% time. Only at least a number of the new parental population (numNeeded)
% has to be sorted.
%  
% input:    F ......... Objective function values (matrix [numInd,numObj])
%           numNeeded . minimum number of individuals to be sorted in fronts
%                       it is not necessary to sort whole population (parents + offsprings)
% output:   rank ...... ranks for solutions the first front has a rank equal to 1, the second 
%                       front equal to 2, etc. 0 is for solutions that are not necessary to sort
%                       to fronts (they are too bad to be sorted)
%
% tested on Octave 6.3.0 (2021-07-11)
% author:  Adela Hlobilova, adela.hlobilova@gmail.com
% version: 3/7/2014 (originally written), 23/2/2022 (last version)

function [rank] = find_ranks(F,numNeeded)
[numInd,numObj] = size(F);

%% find np (domination count) and Sp (set of solutions that the p solution dominates)
np = zeros(numInd,1);   % domination count
Sp = cell(numInd,1);    % set of solutions that the p solution dominates
index = 1:numInd;       % auxiliary value 1:numInd
rank = zeros(numInd,1); % ranks of solutions
numFront = 1;           % the total number of fronts (changed later)

for i=1:numInd
    %% find solutions that dominate ith solution
    % first condition in dominance
    cond1 = bsxfun(@ge,F(i,:),F);
    % second condition in dominance
    cond2 = bsxfun(@gt,F(i,:),F);
    % solutions that dominates ith solution
    dominates = (sum(cond1,2) == numObj)&(sum(cond2,2) > 0);
    % a number of solutions that dominates ith solution
    np(i) = sum(dominates);
    
    %% find solutions that are dominated by ith solution
    % first condition in dominance (but in opposite way)
    cond1 = bsxfun(@le,F(i,:),F);
    % second condition in dominance (but in opposite way)
    cond2 = bsxfun(@lt,F(i,:),F);
    % solutions that are dominated by ith solution
    dominated = (sum(cond1,2) == numObj)&(sum(cond2,2) > 0);
    Sp{i} = index(dominated);
    
end

rank(np==0) = numFront;     % solutions that have np=0 are in the first front
ActualFront = index(np==0); % indexes of individuals in the first front
numAlreadyAdded = length(ActualFront); % number of already assigned individuals

np(np==0) = Inf;  % those solutions are already used and not needed anymore

while numAlreadyAdded<numNeeded
    for i=ActualFront
        np(Sp{i}) = np(Sp{i})-1;
    end
    
    % solutions that have np=0 are in the actual front now
    numFront = numFront+1;  % a number of the actual front
    rank(np==0) = numFront; % assign actual ranks to solutions in actual front
    ActualFront = index(np==0); % indexes of individuals in the actual front
    numAlreadyAdded = numAlreadyAdded + length(ActualFront); % number of already assigned individuals
    np(np==0) = Inf; % and not needed anymore
end


end