function k = WaveZeros(funct,number,upperbound)

dfunct = diff(funct);

%% Find where the function changes sign
vals1 = subs(funct, .001:.001:(upperbound));
vals1 = vals1>0;
vals1 = find(diff(vals1))/1000;

%% Shorten the list so we don't do more operations then needed
if(number<length(vals1))
    vals1 = vals1(1:number);
end
%% Get a better approximation of each point
vals1 = NewtonRaphson(funct,dfunct,vals1, 5);

%% Find where the derivative changes sign
vals2 = subs(dfunct, .001:.001:(upperbound));
vals2 = vals2>0;
vals2 = find(diff(vals2))/1000;
%% Get a better approximation of each point
vals2 = NewtonRaphson(dfunct,diff(dfunct),vals2,5);
%% See if the values are actually zeros of the original function
vals2 = vals2(abs(subs(funct,vals2))<.0001);

%% Check for false positives from derivative.  Comparing to see if it's too
%% close to zeros for the function
temp = vals2;
vals2 = [];
for v = temp
    if(sum(abs(vals1-v)<.01)<1)
        vals2 = [vals2,v];
    end
end

%% Combine lists
valsfinal = sort([vals1,vals2]);

valsfinal = valsfinal(1:number);

%% Check for multiple zeros by taking the derivative and checking if the
%% number is a zero there as well.
temp = valsfinal;
dtempfunct = dfunct;
while(~isempty(temp))
    tempfunct = dtempfunct;
    dtempfunct = diff(tempfunct);
    temp = temp(abs(subs(tempfunct,temp))<.01);
    temp = NewtonRaphson(tempfunct,dtempfunct,temp,5);
    valsfinal = [valsfinal,temp];
end

%% Cleanup
k = sort(valsfinal);
k = k(1:number);