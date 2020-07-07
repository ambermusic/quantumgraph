function [lambda, funct] = WaveEigenvalues(nodes,lengths, number)
%% Look in WaveEquationMatrix for information on the matrix nodes.  lengths
%% is a vector corresponding to the side lengths where the sides are taken
%% in order as the columns of nodes.  Numbers is the number of eigenvalues
%% you want calculated.
[m,L,k] = WaveEquationMatrix(nodes);
funct = subs(det(m),L,lengths);
funct = simplify(funct);
%% If you want just k and not lambda, remove .^2
%% This last entry number*20/sum(lengths) is the upper bound for the
%% eigenvalues it checks. So, if the weyl law closely describes our data,
%% this should be a good enough bound with plenty of error room. You may have to increase it though if
%% something funky happens.
lambda = WaveZeros(funct,number,number*20/sum(lengths)).^2;
