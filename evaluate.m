%% evaluate objectives in NSGA-II (constrained version)
% WARNING:  You have to modify this m-file for a different optimization task!
% input:    X ......... a matrix with individuals with size [numSol,numVar]
%           parallel .. 'on' or 'off'
%           Y ......... objective function values in matrix with size [numSol,numObj]
% tested on Octave 6.3.0 (2021-07-11)
% author:  Adela Hlobilova, adela.hlobilova@gmail.com
% version: 3/7/2014 (originally written), 23/2/2022 (last version)

function [Y] = evaluate(X,parallel)
% numObj = 2; numVar = 1;
%% Schaffer's test function No. 1
%  x is in range [-A,A] where A is from 10 to 10000 (the higher number, the
%  more difficulties are there)
% tested: yes, with range [-1E4 1E4], working nicely
% -----

Y(:,1) = X(:,1).^2;
Y(:,2) = (X(:,1)-2).^2;

%% Schaffer's test function No. 2
% numObj = 2; numVar = 1;
% disconected
% X is in range [-5 10]
% tested: yes, with range [-5 10], working nicely
% -----

% numInd = size(X,1);
% Y = zeros(numInd,2);
% for i=1:numInd
%    if X(i) <= 1
%        Y(i,1) = -X(i);
%    elseif (X(i) > 1) && (X(i) <= 3)
%        Y(i,1) = X(i)-2;
%    elseif (X(i) > 3) && (X(i) <= 4)
%        Y(i,1) = 4-X(i);
%    else
%        Y(i,1) = X(i) - 4;
%    end
%    
%    Y(:,2) = (X-5).^2; 
% end

%% Poloni's two objective function
% numObj = numVar = 2;
% nonconvex, disconected
% X(:,1) is in range [-pi,Inf]
% X(:,2) is in range [-Inf,pi]
% tested: yes, with range +-10, +- pi, working nicely
% -----

% A1 = 0.5*sin(1) - 2*cos(1) +   sin(2) - 1.5*cos(2);
% A2 = 1.5*sin(1) -   cos(1) + 2*sin(2) - 0.5*cos(2);
% B1 = 0.5*sin(X(:,1)) - 2*cos(X(:,1)) +   sin(X(:,2)) - 1.5*cos(X(:,2));
% B2 = 1.5*sin(X(:,1)) - cos(X(:,1)) +   2*sin(X(:,2)) - 0.5*cos(X(:,2));
% 
% Y(:,1) = (1+(A1-B1).^2 + (A2-B2).^2);
% Y(:,2) = (X(:,1)+3).^2 + (X(:,2)+1).^2;

%% Kursawe function
% numObj = 2, numVar = 3;
% nonconvex
% X is in range [-5,5]
% tested: yes, works nicely
% -----
% numInd = size(X,1);
% Y = zeros(numInd,2);
% 
% for i=1:2
%     Y(:,1) = Y(:,1) + (-10*exp(-0.2*sqrt(X(:,i).^2+X(:,i+1).^2)));
% end
% 
% for i=1:3
%     Y(:,2) = Y(:,2) + (abs(X(:,i)).^0.8 + 5*sin(X(:,i).^3));
% end


%% Zitzler-Deb-Thiele's function No. 1
% numObj = 2, numVar = 30;
% convex
% X is in range [0, 1]
% tested: yes, ok
% -----
% numInd = size(X,1);
% Y = zeros(numInd,2);
% 
% Y(:,1) = X(:,1);
% G      = 1 + (9/29 * sum(X(:,2:30),2));
% Y(:,2) = G.*(1-sqrt(Y(:,1)./G));

%% Zitzler-Deb-Thiele's function No. 2
% numObj = 2, numVar = 30;
% convex
% X is in range [0, 1]
% tested:
% -----
% numInd = size(X,1);
% Y = zeros(numInd,2);
% 
% Y(:,1) = X(:,1);
% G      = 1 + (9/29 * sum(X(:,2:30),2));
% Y(:,2) = G.*(1-(Y(:,1)./G).^2);

%% Zitzler-Deb-Thiele's function No. 3
% numObj = 2, numVar = 30;
% convex
% X is in range [0, 1]
% tested:
% -----

% numInd = size(X,1);
% Y = zeros(numInd,2);
% 
% Y(:,1) = X(:,1);
% G      = 1 + (9/29 * sum(X(:,2:30),2));
% Y(:,2) = G.*(1 - sqrt(Y(:,1)./G)-(Y(:,1)./G).*sin(10*pi*Y(:,1)));

%% Zitzler-Deb-Thiele's function No. 4
% numObj = 2, numVar = 10;
% convex
% X is in range [0, 1]
% tested:
% -----
% 
% numInd = size(X,1);
% Y = zeros(numInd,2);
% 
% Y(:,1) = X(:,1);
% G      = 91 + sum(X(:,2:end).^2-10*cos(4*pi*X(:,2:end)),2);
% Y(:,2) = G.*(1 - sqrt(Y(:,1)./G));

%% Zitzler-Deb-Thiele's function No. 6
% numObj = 2, numVar = 10;
% convex
% X is in range [0, 1]
% tested:
% -----

% numInd = size(X,1);
% Y = zeros(numInd,2);
% 
% Y(:,1) = 1-exp(-4*X(:,1)).*sin(6*pi*X(:,1)).^2;
% G      = 1 + 9*(sum(X(:,2:end),2)/9).^2;
% Y(:,2) = G.*(1 - (Y(:,1)./G).^2);

%% Viennet function
% numObj = 3; numVar = 2;
% -3<=X1; X2<=3
% -----

% numInd = size(X,1);
% Y = zeros(numInd,3);
% 
% Y(:,1) = 0.5*(X(:,1).^2 + X(:,2).^2)+sin(X(:,1).^2 + X(:,2).^2);
% Y(:,2) = (3*X(:,1)-2*X(:,2)+4).^2/8 + (X(:,1)-X(:,2)+1).^2/27 + 15;
% Y(:,3) = 1./(X(:,1).^2+X(:,2).^2+1)-1.1*exp(-(X(:,1).^2+X(:,2).^2));

%% Fonseca and Fleming function
% numObj = 2; numVar = from 2 to n;
% X is in range [-4 4]
% -----

% [numInd,numVar] = size(X);
% Y = zeros(numInd,2);
% 
% Y(:,1) = 1 - exp(-sum((X-1/sqrt(numVar)).^2,2));
% Y(:,2) = 1 - exp(-sum((X+1/sqrt(numVar)).^2,2));

%% TEST11 - Four bar truss in Ray, Tai and Seow (2001)
% numObj = 2; numCon = 4; numVar = 4;
% 
% 
% F = 10; E = 2e5; L = 200; sig = 10;
% 
% Y(:,1) = L*(2*X(:,1)+sqrt(2)*X(:,2)+sqrt(2)*X(:,3)+X(:,4));
% Y(:,2) = F*L/E*(2./X(:,1) + 2*sqrt(2)./X(:,2) + 2*sqrt(2)./X(:,3) + 2./X(:,4));
end