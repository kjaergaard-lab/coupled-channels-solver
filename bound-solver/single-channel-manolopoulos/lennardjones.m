function V = lennardjones(x,bi,ni,bo,no)

if nargin<2
    bi = 5;
end
if nargin<3
    ni = 10;
end
if nargin<4
    bo = 10;
end
if nargin<5
    no = 6;
end

V = bi^(ni-2)./x.^ni-bo^(no-2)./x.^no;



end