%% simulated binary cross-over
% for more information see reference: 
% Beyer, H.-G.; Deb, K., "On self-adaptive features in real-parameter ...
% evolutionary algorithms," Evolutionary Computation, IEEE Transactions on ,
% vol.5, no.3, pp.250,270, Jun 2001

% tested on Octave 6.3.0 (2021-07-11)
% author:  Adela Hlobilova, adela.hlobilova@gmail.com
% version: 23/2/2022 (last version)

function offsprings = recombination_bounded_SBX_novec(parents,prob_cross,nu,lbDesVar,ubDesVar)

[numInd,numVar] = size(parents);
if(mod(numInd,2)~=0), error('Number of individuals has to be even!'); end
numPairs = numInd/2;  

% choosing pairs of individuals that are going to be recombinated
CO_MP_index = reshape(randperm(numInd),numPairs,2);

% sorting: p1 < p2
par1 = parents(CO_MP_index(:,1),:) ;
par2 = parents(CO_MP_index(:,2),:) ;

for i=1:numVar
    for j=1:numPairs
        if par1(j,i) > par2(j,i)
            temp = par1(j,i);
            par1(j,i) = par2(j,i);
            par2(j,i) = temp;
        end
    end
end

offsprings = zeros(numInd,numVar);

for i=1:numPairs
    
    if rand(1,1)<prob_cross
        
        u = rand(1,1);
        
        for j=1:numVar
            % offspring 1
            beta_a = 1+(par1(i,j)-lbDesVar(1,j))/(par2(i,j)-par1(i,j));
            gamma_a = 1/(2*beta_a^(nu+1));
            
            if u<=(0.5/(1-gamma_a))
                beta_1 = (2*u*(1-gamma_a))^(1/(nu+1));
            else
                beta_1 = (1/(2*(1-u*(1-gamma_a))))^(1/(nu+1));
            end
            
            offsprings(i,j) = 0.5*(1+beta_1)*par1(i,j) + 0.5*(1-beta_1)*par2(i,j);
            
            % offspring 2
            beta_b = 1+(ubDesVar(1,j)-par2(i,j))/(par2(i,j)-par1(i,j));
            gamma_b = 1/(2*beta_b^(nu+1));
            
            if u<=(0.5/(1-gamma_b))
                beta_2 = (2*u*(1-gamma_b))^(1/(nu+1));
            else
                beta_2 = (1/(2*(1-u*(1-gamma_b))))^(1/(nu+1));
            end
            
            offsprings(i+numPairs,j) = 0.5*(1-beta_2)*par1(i,j)+0.5*(1+beta_2)*par2(i,j);
            
            % returning to boundary: first offspring
            if offsprings(i,j)<lbDesVar(1,j)
               offsprings(i,j) = lbDesVar(1,j);
            end
            if offsprings(i,j)>ubDesVar(1,j)
                offsprings(i,j) = ubDesVar(1,j);
            end
            
            % returning to boundary: second offspring
            if offsprings(i+numPairs,j)<lbDesVar(1,j)
               offsprings(i+numPairs,j) = lbDesVar(1,j);
            end
            if offsprings(i+numPairs,j)>ubDesVar(1,j)
                offsprings(i+numPairs,j) = ubDesVar(1,j);
            end
            
            
        end
    else
        for j=1:numVar
           offsprings(i,j) = par1(i,j);
           offsprings(i+numPairs,j) = par2(i,j);
        end
    end
end

% % new individuals after recombination are accepted with desired probability
% prob_cross_ind = rand(numPairs,1);
% % the rest remains the same
% cross_id=(prob_cross_ind > prob_cross);
% cross_id = repmat(cross_id,2,1);
% CO_MP_index = reshape(CO_MP_index,[],1);
% 
% offsprings(cross_id,:) = parents_backup(CO_MP_index(cross_id),:);

% % returning to the boundary (to fulfill the lower and the upper bound)
% offsprings = bsxfun(@min,offsprings,ubDesVar);
% offsprings = bsxfun(@max,offsprings,lbDesVar);

if any(any(bsxfun(@lt,offsprings,lbDesVar)+ bsxfun(@gt,offsprings,ubDesVar) >= 1))
    fprintf('CHYBA!!!\n');
    
    
    
    
end

end


