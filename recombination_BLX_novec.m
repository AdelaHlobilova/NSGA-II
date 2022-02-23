%% simulated binary cross-over
% for more information see reference: 
% Beyer, H.-G.; Deb, K., "On self-adaptive features in real-parameter ...
% evolutionary algorithms," Evolutionary Computation, IEEE Transactions on ,
% vol.5, no.3, pp.250,270, Jun 2001

% tested on Octave 6.3.0 (2021-07-11)
% author:  Adela Hlobilova, adela.hlobilova@gmail.com
% version: 23/2/2022 (last version)

function offsprings = recombination_BLX_novec(parents,prob_cross,alpha,lbDesVar,ubDesVar)

[numInd,numVar] = size(parents);
if(mod(numInd,2)~=0), error('Number of individuals has to be even!'); end
numPairs = numInd/2;  

% choosing pairs of individuals that are going to be recombinated
CO_MP_index = reshape(randperm(numInd),numPairs,2);

% sorting: p1 < p2
par1 = parents(CO_MP_index(:,1),:) ;
par2 = parents(CO_MP_index(:,2),:) ;

offsprings = zeros(numInd,numVar);

for i=1:numPairs
    
    if rand(1,1)<prob_cross
        
        for j=1:numVar
            d = abs(par1(i,j)-par2(i,j));
            xL = min([par1(i,j), par2(i,j)])-alpha*d;
            xU = max([par1(i,j), par2(i,j)])+alpha*d;
            
            offsprings(i,j) = rand(1,1)*(xU-xL)+xL;
            offsprings(i+numPairs,j) = rand(1,1)*(xU-xL)+xL;
            
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

end


