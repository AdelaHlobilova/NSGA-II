%% Non-Dominated Sorting Genetic Algoritgm II (NSGA-II) according to
%  Deb, Kalyanmoy, et al. "A fast and elitist multiobjective
%  genetic algorithm: NSGA-II." Evolutionary Computation,
%  IEEE Transactions on 6.2 (2002): 182-197.
%
%  A tournament selection, simulated binary cross-over recombination
%  and Gaussian mutation is used to create offsprings. The algorithm
%  is mu+lambda, a new population is chosen from a mixture of parental
%  and offspring generation by ranking. Solutions from the last front
%  are chosen by a crowding distance operator.
%  A whole generation is stored in the matrix with size [numInd,numVar],
%  objective values are stored in the matrix with size [numInd,numObj],
%  constraint values are stores in the matrix with size [numInd,numCon].
%
% input:    numInd ...... number of individuals in each generation (even positive integer)
%           numGen ...... number of generations (positive integer)
%           numVar ...... number of design variables (positive integer)
%           numObj ...... number of objective functions (numObj > 1, integer)
%           lbDesVar .... row vector that contains all minimum allowable values of design variables
%           ubDesVar .... row vector that contains all maximum allowable values of design variables
%           prob_sel .... probability of selection (real value from an interval [0 1])
%           prob_mut .... probability of mutation (real value from an interval [0 1])
%           sig_mut ..... standard deviation of the Gaussian mutation (real positive value)
%           prob_cross .. probability of recombination (real value from an interval [0 1])
%           graphics .... 'on' for graphical outputs, 'off' for no graphical outputs
%           parallel .... 'on' for parallel evaluation of objective functions, 'off' for serial computations
% output:   parents ..... final population containing design variables
%           F_parents ... final population containing objective functions for design variables in parents
%
% tested on Octave 6.3.0 (2021-07-11)
% author:  Adela Hlobilova, adela.hlobilova@gmail.com
% version: 3/7/2014 (originally written), 23/2/2022 (last version)

function [parents,F_parents] = ...
    NSGA_II(numInd,numGen,numVar,numObj,lbDesVar,ubDesVar,prob_sel,prob_mut,sig_mut,prob_cross,param,graphics,parallel,Xover_operator)


%% initial parental population
parents = rand(numInd,numVar);
parents = transform(parents,zeros(1,numVar),ones(1,numVar),lbDesVar,ubDesVar);

% and corresponding objective function values
F_parents = evaluate(parents,parallel);

%% find ranks of the parental population
ranks = find_ranks(F_parents,numInd);

for id_Gen = 1:numGen
    %% offspring population
    crowd_dist = crowding_distance(F_parents,ranks);

    % selection (based on ranks and crowding distances)
    offsprings = selection(parents,ranks,crowd_dist,prob_sel);

    % cross-over (recombination)
    if strcmp(Xover_operator,'SBX')
        nu = param;
        offsprings = recombination_SBX(offsprings,prob_cross,nu,lbDesVar,ubDesVar);
    elseif strcmp(Xover_operator,'boundedSBX')
        nu = param;   % from 2 to 5 for SBX, 
        offsprings = recombination_bounded_SBX_novec(offsprings,prob_cross,nu,lbDesVar,ubDesVar);
    elseif strcmp(Xover_operator,'BLX')
        alpha = param;  % from 0 to 1 for BLX, 0.5 or 0.333 recommended values
        offsprings = recombination_BLX_novec(offsprings,prob_cross,alpha,lbDesVar,ubDesVar);
    else
        error('Wrong definition of cross-over operator!');
    end
    
    % mutation
    offsprings = mutation(offsprings,prob_mut,sig_mut,lbDesVar,ubDesVar);
    
    %% evaluation of objective functions values for offsprings
    F_offsprings = evaluate(offsprings,parallel);
    
    % pool with both generations
    pool = [parents; offsprings];
    F_pool = [F_parents; F_offsprings];
    
    % find ranks for at least numInd solutions from the pool
    ranks = find_ranks(F_pool,numInd);

    %% next parental population
    if length(nonzeros(ranks)) ~= numInd
        % reset all parental staffs (they are already added in pool)
        F_parents = zeros(numInd,numObj);
        parents = zeros(numInd,numVar);
        new_ranks = zeros(numInd,1);
        
        % number of fronts that will be added to next parental population
        % including the last front
        numFronts = max(ranks);
        
        % add all fronts except for the last one (that is bigger)
        indexes = (ranks~=0) & (ranks~= numFronts);
        numAdded = sum(indexes);  % number of added fronts
        parents(1:numAdded,:) = pool(indexes,:);
        F_parents(1:numAdded,:) = F_pool(indexes,:);
        new_ranks(1:numAdded,:) = ranks(indexes,:);
        
        % last front
        LastFront = pool(ranks==numFronts,:);
        F_LastFront = F_pool(ranks==numFronts,:);
        
        % compute crowding distance for the last front
        crowd_dist_LF = crowding_distance(F_LastFront,ones(size(LastFront,1),1));
        % sort the last front according to the crowding distance in descend order
        [~,ind_CDS] = sort(crowd_dist_LF,'descend');
        % number of free slots in the next parental population
        numNeeded = numInd-numAdded;
        
        % add necessary numbers of individuals from the last front
        parents(numAdded+1:numInd,:) = LastFront(ind_CDS(1:numNeeded),:);
        F_parents(numAdded+1:numInd,:) = F_LastFront(ind_CDS(1:numNeeded),:);
        new_ranks(numAdded+1:numInd) = numFronts*ones(numNeeded,1);
        ranks = new_ranks;
        
    else  % crowding distance is useless to compute (numInd == nonzeros(ranks))
        F_parents = F_pool(ranks~=0,:);
        parents = pool(ranks~=0,:);
        ranks = ranks(ranks~=0);
    end
    
if strcmp(graphics,'on')
    
    pause(0.1);
    clf
    figure(1)
    fprintf('Iteration %d\n',id_Gen);
    if numVar==1
        subplot(1,2,1)
        scatter(parents(:,1),zeros(numInd,1),"b",".")
        xlabel('x'), title('Design space');
%         xlim([lbDesVar ubDesVar]); ylim([-1 1]);
        xlim([min(parents(:,1)) max(parents(:,1))]);  ylim([-1 1]);
        set(gca,'ytick',[],'yticklabel',[]);

    elseif numVar == 2
        subplot(1,2,1)
        scatter(parents(:,1),parents(:,2),"b",".")
        box on
        axis equal
        xlabel('x_1'), ylabel('x_2'), title('Design space');
%         xlim([lbDesVar(1) ubDesVar(1)]); ylim([lbDesVar(2) ubDesVar(2)]);
        xlim([min(parents(:,1)) max(parents(:,1))]); ylim([min(parents(:,2)) max(parents(:,2))]); 
    elseif numVar==3
        subplot(1,2,1)
        scatter3(parents(:,1),parents(:,2),parents(:,3),"b",".")
%         hold on
        xlabel('x_1'), ylabel('x_2'), zlabel('x_3'); title('Design space');
        box on
%         axis equal
        xlim([lbDesVar(1) ubDesVar(1)]); ylim([lbDesVar(2) ubDesVar(2)]);
        zlim([lbDesVar(3) ubDesVar(3)]);
%         xlim([min(parents(:,1)) max(parents(:,1))]); ylim([min(parents(:,2)) max(parents(:,2))]); 
%         zlim([min(parents(:,3)) max(parents(:,3))]);
    end
    
    if numObj==2
        if numVar<=3
            subplot(1,2,2)
        end
        scatter(F_parents(:,1),F_parents(:,2),"b",".")
        xlabel('F_1'), ylabel('F_2'), title('Objective space');
        box on
        axis equal
        xlim([floor(min(F_parents(:,1))), ceil(max(F_parents(:,1))) ]); ylim([floor(min(F_parents(:,2))), ceil(max(F_parents(:,2)))]);
    elseif numObj==3
        subplot(1,2,2)
%         hold on
        scatter3(F_parents(:,1),F_parents(:,2),F_parents(:,3),"b",".")
        xlabel('F_1'), ylabel('F_2'), title('Objective space');
        box on
        axis equal
        xlim([floor(min(F_parents(:,1))), ceil(max(F_parents(:,1))) ]); ylim([floor(min(F_parents(:,2))), ceil(max(F_parents(:,2)))]);
        zlim([floor(min(F_parents(:,3))), ceil(max(F_parents(:,3)))]);
    end
    
% hold on
% problem_no = 2;
% analytical_ZDT   
 
end
end
end




