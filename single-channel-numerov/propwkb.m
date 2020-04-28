function idx = propwkb(x,W,idx,direction)
y = 0;
while y>-25 && (idx ~= 2 && idx ~= (numel(W)-1))
    if direction>0
        idx = idx+1;
        y = y-sqrt(W(idx))*(x(idx)-x(idx-1));
    else
        idx = idx-1;
        y = y+sqrt(W(idx))*(x(idx)-x(idx+1));
    end
end

end