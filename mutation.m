%% Standard Gaussian mutation in NSGA-II
% This function mutates solutions with a probability prob_mut by adding a
% random value that is standard normally distributed with mean 0 and
% standard deviation sig_mut. To keep the solution in the feasible design
% space, lbDesVar and ubDesVar are used to return the solution on the boundary, 
% that was mutated to the unfeasible design space.
%  
% input:    offsprings ... all new offspring population (in design space)
%                          matrix with size [numInd,numVar]
%           prob_mut ..... probability of mutation (positive real value [0,1]) 
%           sig_mut ...... standard deviation of Gaussian mutation
%           lbDesVar ..... lower bounds of design variables (a row vector
%                          with numVar size (minX1 minX2 minX3 ... minXnumVar)
%           ubDesVar ..... upper bounds of design variables (similar to lbDesVar)
% output:   newOffsprings ...... mutated offsprings
%
% tested on Octave 6.3.0 (2021-07-11)
% author:  Adela Hlobilova, adela.hlobilova@gmail.com
% version: 7/4/2015 (originally written), 23/2/2022 (last version)

function newOffsprings = mutation(offsprings,prob_mut,sig_mut,lbDesVar,ubDesVar)

[numInd,numVar] = size(offsprings);
newOffsprings = offsprings;

% choosing which individuals are mutated
prob_mut_ind = rand(numInd,1);
[mut_id,~] =find(prob_mut_ind < prob_mut);

%% this version mutates all variables in an appropriate individual 
%  a probability of mutation is defined by user in prob_mut variable
%  EXAMPLE: A problem has a dimension equal to 10 (10 variables, numVar =
%  10). The k-th individual (from total numInd individuals) is decided by
%  an algorithm to be mutated with a probability prob_mut. Therefore all 
%  X_k_1, X_k_2, ..., X_k_10 variables of k_th individual are changed 
%  with a probability equal to 1.

if ~isempty(mut_id)
    for i=1:length(mut_id)
%         mut_var = floor(rand(1)*numVar)+1;  % choosing which variable is mutated
        newOffsprings(mut_id(i),:) = offsprings(mut_id(i),:)+normrnd(zeros(1,numVar),sig_mut,1,numVar);
    end
end

%% this version mutates only selected variables in an appropriate individual
%  a probability of mutation of an individual is defined by user in prob_mut variable 
%  a probability of mutation of the individual variables is equal to 1/numVar (previous version 0.5)
%  EXAMPLEs: A problem has a dimension equal to 10 (10 variables, numVar =
%  10). The k-th individual (from total numInd individuals) is decided by
%  an algorithm to be mutated with a probability prob_mut. The probability 
%  of mutation of each variable is not equal to 1 as in version above but 
%  to 1/numVar. X_k_1, X_k_2, ..., X_k_10 variables of k_th individual are 
%  changed with a probability 1/numVar. 
%  NOTE: This version is less violent and in some cases can help to find
%  better optima.

% if ~isempty(mut_id)
%     for i=1:length(mut_id)
%         for j=1:numVar
%             if rand(1) <= 1/numVar    % previous version rand(1) <= 0.5 
%                 newOffsprings(mut_id(i),j) = offsprings(mut_id(i),j)+normrnd(0,sig_mut(j));
%             end
%         end
%     end
% end

% returning to the boundary
newOffsprings = bsxfun(@min,newOffsprings,ubDesVar);
newOffsprings = bsxfun(@max,newOffsprings,lbDesVar);

end