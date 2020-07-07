%% f is the function we test, x is the initial guess, steps is the number
%% of iterations.  We return 
function x = NewtonRaphson(f,df,x,steps)
for n = 1:steps
    x = x - (subs(f,x)./subs(df,x));
end
x(abs(subs(f,x))>.01) = NaN;