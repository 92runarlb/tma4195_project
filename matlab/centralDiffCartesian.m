function gradEta = centralDiffCartesian(eta, Nx, Ny, node_coords)
    N = Nx*Ny;
    
    
    mainDiag = zeros(Nx,1);
    mainDiag([1,end]) = [-2,2];
    mainDiag = repmat(mainDiag,Ny,1);
    subDiag = -ones(Nx,1);
    subDiag(Nx-1) = -2;
    subDiag(end) = 0;
    subDiag = repmat(subDiag,Ny,1);
    supDiag = -subDiag(end:-1:1);
    Ax = Nx*0.5*spdiags([subDiag,mainDiag,supDiag],[-1,0,1],N,N);

    detadx= Ax*eta;
    
    %Calculate detay
    
    mainDiag = zeros(Ny,1);
    mainDiag([1,end]) = [-2,2];
    mainDiag = repmat(mainDiag,Nx,1);
    subDiag = -ones(Ny,1);
    subDiag(Ny-1) = -2;
    subDiag(end) = 0;
    subDiag = repmat(subDiag,Nx,1);
    supDiag = -subDiag(end:-1:1);
    Ay = Ny*0.5*spdiags([subDiag,mainDiag,supDiag],[-1,0,1],N,N);
    
    %Sort y in the right order. That is, in the order of x.
    sort = [1:Nx:N-Nx+1]';
    sort = repmat(sort,Nx,1);
    addMat = repmat(0:Nx-1,Ny,1);
    addMat = addMat(:);
    sort = sort + addMat;

    detady = Ay*eta(sort);
    
    % Sort back the eta's
    detady = vec2mat(detady,Ny);
    detady = detady(:);
    gradEta = [detadx,detady];
    

    
    
end