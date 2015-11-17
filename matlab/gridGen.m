function G = gridGen(dim, gridLim, h, eta)

nx = dim(1); ny = dim(2);
xMax = gridLim(1); yMax = gridLim(2);

P = [];
for x = 0:xMax/nx:xMax
    diff = h(x)+eta(x);
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

nc = G.cells.num;
remCells = false(nc,1);

for c = 1:G.cells.num
    if G.cells.centroids(c,2) < -h(G.cells.centroids(c,1)) | ...
       G.cells.centroids(c,2) > eta(G.cells.centroids(c,1))
        remCells(c) = true;
    end
end
G = removeCells(G,find(remCells));
G = pebi(G);
G = computeGeometry(G);
hold on
plotGrid(G);
plot(P(:,1),P(:,2),'.')
hold off

end