function [Eout,u,dbg,nodes] = solvebound(r,Vfunc,Ein,ops,opt)

%% Parse arguments
if nargin<5
    opt = boundoptions;
elseif ~isa(opt,'boundoptions')
    error('Options argument ''opt'' must be of type boundoptions');
end


E = sort(Ein);
Enew = E(1);
match = [calcBoundSolution(r,Vfunc,E(1),ops,opt),calcBoundSolution(r,Vfunc,E(2),ops,opt)];
side = 0;
if opt.debug
    fprintf(1,'Iter: %02d, Error = %.5e, [E1,E2] = [%.3f,%.3f], Energy = %.5f\n',00,10,E,NaN);
end

usefsolve = false;

for nn=1:opt.iter
    %%
    Eold = Enew;
    if sign(match(1)) == sign(match(2))
        Enew = (E(1)+E(2))/2;
        method = 'bisect';
    else
        Enew = E(1)-match(1)*(E(1)-E(2))/(match(1)-match(2));
        method = 'secant';
    end
    
    if opt.debug
        [mnew,~,dbg] = calcBoundSolution(r,Vfunc,Enew,ops,opt);
    else
        mnew = calcBoundSolution(r,Vfunc,Enew,ops,opt);
    end
    
    err = abs(mnew);
    dE = abs(Enew-Eold)/abs(Enew+Eold);
    
    if err<opt.tolF && dE<opt.tolE
        break;
%     elseif dE<opt.tolE && err>opt.tolF
%         usefsolve = true;
%         break;
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
            if amatch(1)>amatch(2)
                maxIndex = 1;
            else
                maxIndex = 2;
            end
            match(maxIndex) = mnew;
            E(maxIndex) = Enew;
%             match(2) = mnew;
%             E(2) = Enew;
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
%         figure(1);clf;
%         plot(debugOut.rL,debugOut.uL,'.-');
%         hold on;
%         plot(debugOut.rR,debugOut.uR,'.-');
%         hold off;
%         xlim(min(debugOut.rR)+[-1,1]);
        fprintf(1,'Iter: %02d, Error = %.5e, [E1,E2] = [%#.3f,%#.3f], Energy = %.5f\n',nn,err,E,Enew);
        pause(opt.pauseDelay);
    end
end

if dE>opt.tolE || err>opt.tolF
    usefsolve = true;
end

if usefsolve
    if opt.debug
        fprintf(1,'Basic solver failed. Using fsolve()...\n');
        options = optimset('display','iter','tolX',opt.tolE,'maxiter',opt.iter);
    else
        options = optimset('display','off','tolX',opt.tolE,'maxiter',opt.iter);
    end
    Enew = fsolve(@(E) calcBoundSolution(r,Vfunc,E,ops,opt),max(Ein),options);
end
    

Eout = Enew;
if opt.debug || nargout>1
    [~,nodes,dbg] = calcBoundSolution(r,Vfunc,Eout,ops,opt);
    uL = zeros(size(dbg.zL,1),size(dbg.zL,3));
    uL(:,end) = dbg.wf;
    for kk=size(uL,2):-1:2
        uL(:,kk-1) = dbg.zL(:,:,kk)*uL(:,kk);
    end
    
    uR = zeros(size(dbg.zR,1),size(dbg.zR,3));
    uR(:,end) = dbg.wf;
    for kk=size(uR,2):-1:2
        uR(:,kk-1) = dbg.zR(:,:,kk)*uR(:,kk);
    end
    uR = flip(uR,2);
    
    u = [uL uR(:,2:end)];
    dr = repmat(diff(r(:)'),size(u,1),1);
    u = u./sqrt(sum(sum(abs(u(:,1:end-1)).^2.*dr,2),1));   
end

if opt.debug
%     cla;
%     plot(debugOut.rL,debugOut.uL,'.-');
%     hold on;
%     plot(debugOut.rR,debugOut.uR,'.-');
%     hold off;
    figure(1);clf;
    plot(r,u,'.-');
    fprintf(1,'---------------------------------------\n');
    fprintf(1,'Iter: %02d, Error = %.5e, Energy = %.5f\n',nn,err,Enew);
    pause(opt.pauseDelay);
end


end






