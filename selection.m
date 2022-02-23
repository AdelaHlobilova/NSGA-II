%% Tournament selection in NSGA-II
% A tournament selection serves for keeping best individuals in the
% population. There are *three laps* for the tournament fight for every
% randomly chosen pair of individuals. The first competition lap is in
% ranking. The solution with better rank (lower rank) is placed in the next
% offspring population. In case that ranks are the same for both solutions,
% the pair competes in the crowding distances (the higher, the better). In
% case that the crowding distances are the same, randomly chosen individual
% wins. The next chosen pair can contain the worse parent from the previous
% figth and it is possible to choose it to the next offspring population,
% whether it is better in the next fight. It means that we do not kill
% worse after the imediate fight but we return it to the mating pool
% keeping a chance to it in next fight. For more information see book
% Burke: Search methodologies (chap. GA).
%
% Input:    parents ......... parental population
%           ranks ........... ranking for the parental population (for more
%                             details see file find_ranks.m)
%           crowd_dist ...... crowding distances for each front (for more
%                             details see file crowding_distance.m)
%           prob_sel ........ probability of selection (real number in [0,1])
% Output:   offsprings ...... new offspring population
%
% tested on Octave 6.3.0 (2021-07-11)
% author:  Adela Hlobilova, adela.hlobilova@gmail.com
% version: 3/7/2014 (originally written), 23/2/2022 (last version)
function offsprings = selection(parents,ranks,crowd_dist,prob_sel)

[numInd,~] = size(parents);

% tournament selection (TS)
TS_MP_index = ceil(rand(numInd,2)*numInd);
REP_TS_MP_index = (TS_MP_index(:,1) == TS_MP_index(:,2));

% not choosing the same solution for the tournament
if sum(REP_TS_MP_index) > 0
    pos_REP_TS_MP_index = find(REP_TS_MP_index);
    for i=1:length(pos_REP_TS_MP_index)
        while (TS_MP_index(pos_REP_TS_MP_index(i),1) == TS_MP_index(pos_REP_TS_MP_index(i),2))
            TS_MP_index(pos_REP_TS_MP_index(i),:) =  ceil(rand(1,2)*numInd);
        end
    end
end

% tournament selection mating pool
% TS_MP_F = [F_parents(TS_MP_index(:,1)) F_parents(TS_MP_index(:,2))];
TS_MP_Ranks = [ranks(TS_MP_index(:,1)) ranks(TS_MP_index(:,2))];
TS_MP_CrDist = [crowd_dist(TS_MP_index(:,1)) crowd_dist(TS_MP_index(:,2))];

%% who will win the tournament selection?
% --- first lap of the fight
% find solutions with better ranks
[~,TS_winners] = min(TS_MP_Ranks,[],2); % Tournament selection winners!!!

% zero results with equal ranks (they will fight in another lap)
next_fighters = (TS_MP_Ranks(:,1) == TS_MP_Ranks(:,2));
TS_winners(next_fighters) = 0;

% --- second lap of the tournament
[~,TS_winners(next_fighters)] = max(TS_MP_CrDist(next_fighters,:),[],2);

% zero results with equal ranks and equal crowding distances
next2_fighters = false(numInd,1);
next2_fighters(next_fighters) = TS_MP_CrDist(next_fighters,1) == TS_MP_CrDist(next_fighters,2);
TS_winners(next2_fighters) = 0;

% --- third lap of the tournament
% if ranks and crowding distances are equal just choose one solution randomly
TS_winners(next2_fighters) = round(rand(1))+1;

% --- end of all fights

%% choose the other individual and not the best one with (1-prob_sel)
change_ind = bsxfun(@ge,rand(numInd,1),prob_sel);
temp = 1:numInd;
ind_num = temp(change_ind);
one = ind_num(TS_winners(ind_num)==1);
two = ind_num(TS_winners(ind_num)==2);
TS_winners(one) = 2;
TS_winners(two) = 1;

%% assign TS_winners indexes to correct individuals
temp = zeros(numInd,1);
for i=1:numInd
    temp(i) =  TS_MP_index(i,TS_winners(i));
end

%% create offspring population
offsprings = parents(temp,:);

end