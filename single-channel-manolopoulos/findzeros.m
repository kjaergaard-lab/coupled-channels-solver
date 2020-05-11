function [zi,zdi] = findzeros(W)

zi = [];
zdi = [];
numzeros = 0;
for nn=1:numel(W)-1
    if sign(W(nn+1)) ~= sign(W(nn))
        if numzeros == 0 || nn ~= (zi(numzeros)+1)
            numzeros = numzeros+1;
            zi(numzeros) = nn; %#ok<AGROW>
            if (W(nn+1)-W(nn)) > 0
                zdi(numzeros) = 1;
            else
                zdi(numzeros) = -1;
            end
        end
    end
end



end