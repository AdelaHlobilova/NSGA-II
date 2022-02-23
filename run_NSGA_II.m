%% run NSGA_II
% this is a sample code for running the NSGA-II optimization algorithm
% tested on Octave 6.3.0 (2021-07-11)
% author:  Adela Hlobilova, adela.hlobilova@gmail.com
% version: 3/7/2014 (originally written), 23/2/2022 (last version)

pkg load statistics

%% initialization
numInd = 100;  % number of individuals in one generation
numGen = 20;  % number of generations (the initial population is not included in this number)
numObj = 2;    % number of objective (cost/loss) functions, this is a multi-objective algorithm, 
               % the recommendation is to have maximum 3 objectives, for more objectives select
               % a different algorithm (many-objectives would be better, e.g. NSGA-III)
numVar = 1;    % number of design variables

%% bounds of the problem

% lbDesVar = [0, -5*ones(1,numVar-1)];  % a row vector that contains all minimum allowable values of Des. Vars
% ubDesVar = [1,  5*ones(1,numVar-1)];  % -------------- || ------------ maximum ---------------- || ---------
% lbDesVar = [-3 -30];
% ubDesVar = [30 3];
% lbDesVar = [-pi,-10];
% ubDesVar = [10,pi];
lbDesVar = -10; ubDesVar = 10;
% lbDesVar = -5*ones(1,numVar);
% ubDesVar = 5*ones(1,numVar);
% lbDesVar = 0*ones(1,numVar);
% ubDesVar = 1*ones(1,numVar);
% lbDesVar = -5; ubDesVar = 10;

%% parameters for genetic algorithm
prob_mut = 1/numVar;   % probability of mutation
% sig_mut = 0.2*mean([lbDesVar; ubDesVar]);      % standard deviation of the Gaussian mutation
sig_mut = 0.2*(ubDesVar-lbDesVar);
% sig_mut =ones(1,numVar);       % standard deviation of the Gaussian mutation

prob_cross = 0.9; % probability of cross-over
prob_sel = 1;     % probability of selection

Xover_operator = 'boundedSBX';   % crossover operator, options: 'SBX', 'boundedSBX', 'BLX'
param = 0.5;              % parameter affecting the crossover operator
                          % recommendation: from 2 to 5 for SBX; from 0 to 1 for BLX, 0.5 or 0.333 recommended values

graphics = 'on';
parallel = 'off';
                          
[parents,F_parents] =  NSGA_II(numInd,numGen,numVar,numObj, ...
    lbDesVar,ubDesVar,prob_sel,prob_mut,sig_mut,prob_cross,param,graphics,parallel,Xover_operator);

