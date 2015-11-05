function gradPhi = centralDiffCartesian(G, phi)

    DIM = G.cartDims;
    Nx = DIM(1);
    Ny = DIM(2);
    N = Nx*Ny;

    X = G.cells.centroids;

    endpoints = zeros(Nx,1);
    endpoints([1,end]) = [-2,2];
    endpoints = repmat(endpoints,Ny,1);
    subDiag = -ones(Nx,1);
    subDiag(Nx-1) = -2;
    subDiag(end) = 0;
    subDiag = repmat(subDiag,Ny,1);
    supDiag = -subDiag(end:-1:1);

    A = 0.5*Nx*spdiags([subDiag,endpoints,supDiag],[-1,0,1],N,N);

    
    sort = [1:Nx:N-Nx+1]';
    sort = repmat(sort,Nx,1);
    addMat = repmat(0:Nx-1,Nx,1);
    addMat = addMat(:);
    sort = sort + addMat;

    dphidy = A*phi(sort);
    dphidy = dphidy(sort);

    dphidx= A*phi;
    
    gradPhi = [dphidx,dphidy];
%     
%     x = linspace(0,1,Ngrid);
%     [XM,YM] = meshgrid(x,x);
%     Z = griddata(X(:,1),X(:,2),sum((gradPhi).^2,2),XM,YM);
% 
%     surf(XM, YM, Z)
%     xlabel('x-axis')
%     ylabel('y-axis')
%     zlabel('z-axis')
%     
end