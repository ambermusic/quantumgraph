function [m,L,k] = WaveEquationMatrix(nodes) %#ok<STOUT>
%% WaveEquationMatrix3 has the main description, this function is an
%% optimization where the resulting matrix's rows and columns are halved
%% (except for the arbitrary waves)
%% by taking two different optimizations:
%% 1) When an edge is at a boundary pinned to zero we take the equation for
%% the edge to be sin(k*(L_n - x)) with direction towards origin
%% 2) For loops with both ends connected to the same node we write the wave
%% along it as cos(k*(L_n/2-x)) instead of connecting to edges along the
%% loop.  Previously with a loop we would have [1,1;1,1] but now we can
%% write this as simply [2]. We also must multiple the final expression by
%% sin(L_n * x) to get the other zeros out
%% ALSO, one additional change was implemented: instead of boundary nodes defined by
%% the edge being on a row by itself, it is now just implied if the edge is
%% not connected to two nodes.  So for the balloon graph instead of
%% [1,0,0;0,1,1;1,1,1] it is now [1,2] denoting the one node where the
%% edges are connected.  One connected at both ends, the other with one.
%% EXAMPLES:
%% Balloon graph:       [1,2]
%% Dumbbell graph:      [2,1,0;0,1,2]
%% Y-graph:             [1,1,1]
%% 'Kilroy' graph:      [1,1,1,0;0,1,1,1]
%% Counterexample:      [1,1,1,1,1,1;0,1,1,1,1,1]  w/ Lengths
%%                                                   [pi,pi,pi,pi,pi,pi]
numvars = size(nodes,2);
%% Determines which edges have boundaries pinned to zero
bounds = find(sum(nodes,1)==1);
%% Determines which edges are loops
[bogus, loops] = find(nodes==2);
%% Determine which edges are arbitrary
% arb = 1:numvars;
% arb([loops' bounds]) = [];
%% Creates an array of symbolic variables for the lengths
L = sym(zeros(numvars, 1));

%% Determine orientation of sides.  We will take the first node in order
%% that each string is attached to as its origin.  We only care about this
%% for the arbitrary strings.
Origins = zeros(numvars,1);
for n = 1:numvars
    L(n,1) = sym(sprintf('L%d', n));
    Origins(n) = find(nodes(:,n),1,'first');
end

%% Square root of eigenvalue
syms k;
%% Position along edge
syms x;
waves = sym(zeros(numvars,2));
dwaves = sym(zeros(numvars,2));
%% The general form of the waves assuming no simplifications
for n = 1:numvars
    if(any(ismember(bounds,n)))
        waves(n,:) = [sin(k*(L(n)-x)),0];
        dwaves(n,:) = [-cos(k*(L(n)-x)),0];
    elseif(any(ismember(loops',n)))
        waves(n,:) = [cos(k*(L(n)/2-x)),0];
        dwaves(n,:) = [2*sin(k*(L(n)/2-x)),0];
    else
        waves(n,:) = [sin(x*k),cos(x*k)];
        dwaves(n,:) = [cos(x*k),-sin(x*k)];
    end
end
m = sym(zeros(2*numvars));
index = 1;
%% Boundary constraints // Section is no longer necessary
%% Equality constraints
for n = 1:size(nodes,1)
    edges = find(nodes(n,:));
    primary = edges(1);
    if(Origins(primary)==n)
        length = 0;
    else
        length = L(primary);
    end
    for edge = edges(2:end)
        if(Origins(edge)==n)
             m(index,2*edge-1:2*edge) = -subs(waves(edge,:),x,0);
             m(index,2*primary-1:2*primary) = subs(waves(primary,:),x,length);
        else
             m(index,2*edge-1:2*edge) = -subs(waves(edge,:),x,L(edge));
             m(index,2*primary-1:2*primary) = subs(waves(primary,:),x,length);
        end
        index = index + 1;
    end
    
    %% Kirchhoff condition
    for edge = edges
        if(Origins(edge)==n)
             m(index,2*edge-1:2*edge) = subs(dwaves(edge,:),x,0);
        else
             m(index,2*edge-1:2*edge) = -subs(dwaves(edge,:),x,L(edge));
        end
    end
    index = index + 1;
end
m(~any(m~=0,2),:) = [];
m(:,~any(m~=0,1)) = [];
if( ~isempty(loops))
    for n = loops'
        m(end,:) = m(end,:).*sin(L(n)/2*k);
    end
end