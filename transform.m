%% linear transformation of data from one space to another
%  inputs:
%    X ..... original data (each row corresponds with one experiment and
%            each column corresponds with one variable) [numExp,Dim] = size(X)
%    Omin .. lower bound of original data (can be scalar or row vector)
%    Omax .. upper bound of original data (can be scalar or row vector)
%    Nmin .. lower bound of shifted data (can be scalar or row vector)
%    Nmax .. upper bound of shifted data (can be scalar or row vector)
%  outputs:
%    X_shifted ... shifted data (each row corresponds with one experiment and
%      each column corresponds with one variable) [numExp,Dim] = size(X)
%    
%
% example: shifting LHS design in uniform space to uniform space from -5 to
% 3 and 12 to 18 (two variables)
% X = lhsdesign(5,2);
% X_shifted = transform(X,[0 0],[1 1],[-5 12],[3 18]); 
%
% tested on Octave 6.3.0 (2021-07-11)
% author:  Adela Hlobilova, adela.hlobilova@gmail.com
% version: 3/7/2014 (originally written), 23/2/2022 (last version)

function X_shifted = transform(X,Omin,Omax,Nmin,Nmax)

if (length(Nmax) ~= length(Nmin) || length(Omax) ~= length(Omin) || length(Omax) ~= length(Nmin))
    error('Omin, Omax, Nmin, Nmax have to have the same size\n');
end

ratio = (Nmax-Nmin)./(Omax-Omin);

X_shifted = bsxfun(@minus,X,Omin);
X_shifted = bsxfun(@times,X_shifted,ratio);
X_shifted = bsxfun(@plus,X_shifted,Nmin);

end


