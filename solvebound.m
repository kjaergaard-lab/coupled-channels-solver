function [Eout,u] = solvebound(r,Vfunc,E,options)



narg = 4;

if nargin<narg || ~isfield(options,'maxIter')
    options.maxIter = 100;
end
if nargin<narg || ~isfield(options,'tol')
    options.tol = 1e-6;
end
if nargin<narg || ~isfield(options,'debug')
    options.debug = 0;
end
if nargin<narg || ~isfield(options,'pauseTime')
    options.pauseTime = 0.1;
end


E = sort(E);
match = [calcBoundSolution(r,Vfunc,E(1),options),calcBoundSolution(r,Vfunc,E(2),options)];
side = 0;

for nn=1:options.maxIter
    %%
    if sign(match(1)) == sign(match(2))
        Enew = (E(1)+E(2))/2;
        method = 'bisect';
    else
        Enew = E(1)-match(1)*(E(1)-E(2))/(match(1)-match(2));
        method = 'secant';
    end
    
    if options.debug
        [mnew,~,debugOut] = calcBoundSolution(r,Vfunc,Enew,options);
    else
        mnew = calcBoundSolution(r,Vfunc,Enew,options);
    end
    
    if abs(mnew)<options.tol
        break;
    elseif strcmp(method,'bisect')
        amnew = abs(mnew);
        amatch = abs(match);
        if sign(mnew) ~= sign(match(1))
            if amnew<amatch(1)
                match(1) = mnew;
                E(1) = Enew;
            elseif amnew<amatch(2)
                match(2) = mnew;
                E(2) = Enew;
            else
                if amatch(1)>amatch(2)
                    match(1) = mnew;
                    E(1) = Enew;
                else
                    match(2) = mnew;
                    E(2) = Enew;
                end
            end
        elseif amnew<amatch(1) && amnew>amatch(2)
            match(1) = mnew;
            E(1) = Enew;
        elseif amnew<amatch(2) && amnew>amatch(1)
            match(2) = mnew;
            E(2) = Enew;
        else
            match(2) = mnew;
            E(2) = Enew;
        end    
    elseif strcmp(method,'secant')
        if sign(mnew) == sign(match(1))
            E(1) = Enew;
            match(1) = mnew;
            if side == -1
                match(2) = match(2)/2;
            end
            side = -1;
        elseif sign(mnew) == sign(match(2))
            E(2) = Enew;
            match(2) = mnew;
            if side == 1
                match(1) = match(1)/2;
            end
            side = 1;
        end
    end
    
    if options.debug
        figure(1);clf;
%         plot(x,u,'.-');
        plot(debugOut.rL,debugOut.uL,'.-');
        hold on;
        plot(debugOut.rR,debugOut.uR,'.-');
        hold off;
%         xlim(min(debugOut.xR)+[-0.1,0.1]);
        fprintf(1,'Iter: %02d, Error = %.5e, Energy = %.5f\n',nn,abs(mnew),Enew);
        pause(options.pauseTime);
    end
    
end

Eout = Enew;
if ~options.debug || nargout>1
    [~,u,debugOut] = calcBoundSolution(r,Vfunc,Eout,options);
end

if options.debug
    cla;
    plot(debugOut.rL,debugOut.uL,'.-');
    hold on;
    plot(debugOut.rR,debugOut.uR,'.-');
    hold off;
    fprintf(1,'Iter: %02d, Error = %.3e, Energy = %.3f\n',nn,abs(mnew),Enew);
    pause(options.pauseTime);
end


end






