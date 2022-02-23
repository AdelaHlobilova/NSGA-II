%% Crowding distance in NSGA-II
% A crowding distance is another criterion how to efficiently find better
% solutions in the same front. The higher the distance between
% solutions is, the better diversity is in the population. This criterion
% therefore prefers solutions that are at the edges of the fronts and then
% also prefers solutions that have a big distance to other solutions. On
% the other hand, it tries to eliminate solutions that are very close to
% each other in the objective space. An elimination is the sense of the
% very small crowding distance.
%
%  WARNING: It computes crowding distance according to each front! If you want 
%  to compute crowding_distance in overall population for other purposes than 
%  NSGA-II, use ones(numInd,1) instead of ranks in calling this function.
%  
% input:    F ............ objective function values (matrix [numInd,numObj])
%           ranks ........ evaluation of the dominance (the smaller, the
%                          better), for more details see find_ranks.m
% output:   crowd_d_ALL .. crowding distances for each front
%
% tested on Octave 6.3.0 (2021-07-11)
% author:  Adela Hlobilova, adela.hlobilova@gmail.com
% version: 6/3/2015 (originally written), 23/2/2022 (last version)

function [crowd_d_ALL] = crowding_distance(F,ranks)

numFronts = max(ranks);
numAllInd = size(F,1);
crowd_d_ALL = zeros(numAllInd,1);

for iFront = 1:numFronts
    
    % choose solutions in the iFront-th front and add it to an independent vector
    F_iFront = F(ranks==iFront,:);
    [numInd,numObj] = size(F_iFront);
    % sort individuals according to all objectives (ranking is in ind_CDS)
    [~,ind_CDS]  = sort(F_iFront,1);

    crowd_d = zeros(numInd,1);  % sorted as the population in individuals
    for i=1:numObj
        % solutions with the lowest and the highest values of the i-th objective
        % function are prefered with Inf value for their high diversity
        crowd_d(ind_CDS(1,i),1) = Inf;
        crowd_d(ind_CDS(end,i),1) = Inf;
        % maxval and minval serve for the scaling purposes
        maxval = max(F_iFront(:,i));
        minval = min(F_iFront(:,i));
        for j=2:numInd-1
            % scaled distance for the j-th solution from the (j-1)-th to
            % (j+1)-th solution in objective space (cumulated across all
            % objective functions)
            crowd_d(ind_CDS(j,i),1) = crowd_d(ind_CDS(j,i),1)+ ...
                abs(((F_iFront(ind_CDS(j+1,i),i)-F_iFront(ind_CDS(j-1,i),i)))/ ...
                (maxval - minval));
        end
    end
    % localization of crowding distances according to global sorting in F
    crowd_d_ALL(ranks==iFront) = crowd_d;
end
end