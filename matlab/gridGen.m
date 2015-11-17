function G = gridGen(dim, gridLim, h, eta)

nx = dim(1); ny = dim(2);
xMax = gridLim(1); yMax = gridLim(2);

P = [];
for x = 0:xMax/nx:xMax
    diff = +h(x)+eta(x);
    n = ceil(ny*diff);
    Y = (linspace(-h(x),eta(x), n))';
    %Y = (-h(x):yMax/ny:eta(x))';
    X = repmat(x,size(Y,1),1);
    P = [P;X,Y];
end
T = delaunay(P(:,1),P(:,2));
G = triangleGrid(P,T);
%G = pebi(G);

G = computeGeometry(G);

plotGrid(G);

end