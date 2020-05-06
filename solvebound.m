function [Eout,u,nodes] = solvebound(r,Vfunc,E,opt)

%% Parse arguments
if nargin<4
    opt = boundoptions;
elseif ~isa(opt,'boundoptions')
    error('Options argument ''opt'' must be of type boundoptions');
end


E = sort(E);
match = [calcBoundSolution(r,Vfunc,E(1),opt),calcBoundSolution(r,Vfunc,E(2),opt)];
side = 0;

for nn=1:opt.iter
    %%
    if sign(match(1)) == sign(match(2))
        Enew = (E(1)+E(2))/2;
        method = 'bisect';
    else
        Enew = E(1)-match(1)*(E(1)-E(2))/(match(1)-match(2));
        method = 'secant';
    end
    
    if opt.debug
        [mnew,~,debugOut] = calcBoundSolution(r,Vfunc,Enew,opt);
    else
        mnew = calcBoundSolution(r,Vfunc,Enew,opt);
    end
    
    err = abs(mnew);
    
    if err<opt.tol
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
    
    if opt.debug
        figure(1);clf;
%         plot(x,u,'.-');
        plot(debugOut.rL,debugOut.uL,'.-');
        hold on;
        plot(debugOut.rR,debugOut.uR,'.-');
        hold off;
        xlim(min(debugOut.rR)+[-1,1]);
%         xlim([8,13]);
        fprintf(1,'Iter: %02d, Error = %.5e, Energy = %.5f\n',nn,err,Enew);
        pause(opt.pauseDelay);
    end
    
end

Eout = Enew;
if ~opt.debug || nargout>1
    [~,nodes,debugOut] = calcBoundSolution(r,Vfunc,Eout,opt);
    u = debugOut.u;
    u = u./sqrt(sum(u(1:end-1).^2.*diff(r)));
end

if opt.debug
    cla;
    plot(debugOut.rL,debugOut.uL,'.-');
    hold on;
    plot(debugOut.rR,debugOut.uR,'.-');
    hold off;
    fprintf(1,'Iter: %02d, Error = %.5e, Energy = %.5f\n',nn,err,Enew);
    pause(opt.pauseDelay);
end


end






