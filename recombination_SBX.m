%% simulated binary cross-over in NSGA-II
% for more information see reference: 
% Beyer, H.-G.; Deb, K., "On self-adaptive features in real-parameter ...
% evolutionary algorithms," Evolutionary Computation, IEEE Transactions on ,
% vol.5, no.3, pp.250,270, Jun 2001
%
% inputs:   offsprings ..... new offspring population ([numInd,numVar] size)
%           prob_cross ..... probability of cross-over (recombination)
%                            (real value in [0,1])
%           nu ............. distance parameter of cross-over distribution
%           lbDesVar ....... lower bounds of design variables (a row vector
%                            with numVar size (minX1 minX2 minX3 ... minXnumVar)
%           ubDesVar ....... upper bounds of design variables (similar to lbDesVar)
% outputs:  new_off ........ new offspring generation after recombination
%
% tested on Octave 6.3.0 (2021-07-11)
% author:  Adela Hlobilova, adela.hlobilova@gmail.com
% version: 24/9/2015 (originally written), 23/2/2022 (last version)

function new_off = recombination_SBX(offsprings,prob_cross,nu,lbDesVar,ubDesVar)

[numInd,numVar] = size(offsprings);
if(mod(numInd,2)~=0), error('Number of individuals has to be even!'); end
numPairs = numInd/2;  

% choosing pairs of individuals that are going to be recombinated
CO_MP_index = reshape(randperm(numInd),numPairs,2);

% simulated binary crossover (SBX)
u_for_beta = rand(numPairs,1);
beta = zeros(numPairs,1);

% nu is distance parameter of cross-over distribution
Case1 = (u_for_beta <=0.5);
beta(Case1) = (2*u_for_beta(Case1)).^(1/(nu+1));
Case2 = (u_for_beta > 0.5);
beta(Case2) = (2*(1-u_for_beta(Case2))).^(-1/(nu+1));

beta = kron(ones(1,numVar),beta);

new_off = zeros(2*numPairs,numVar);
new_off(1:numPairs,:) = 1/2*((1-beta).*offsprings(CO_MP_index(:,1),:)+(1+beta).*offsprings(CO_MP_index(:,2),:));
new_off(numPairs+1:end,:) = 1/2*((1+beta).*offsprings(CO_MP_index(:,1),:)+(1-beta).*offsprings(CO_MP_index(:,2),:));

% new individuals after recombination are accepted with desired probability
prob_cross_ind = rand(numPairs,1);
% the rest remains the same
cross_id=(prob_cross_ind > prob_cross);
cross_id = repmat(cross_id,2,1);
CO_MP_index = reshape(CO_MP_index,[],1);

new_off(cross_id) = offsprings(CO_MP_index(cross_id));

% returning to the boundary (to fulfill the lower and the upper bound)
new_off = bsxfun(@min,new_off,ubDesVar);
new_off = bsxfun(@max,new_off,lbDesVar);

end


