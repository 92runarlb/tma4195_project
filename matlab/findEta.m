function eta = findEta(eta, gradPhi, h, k)
    assert(size(eta,1)==size(gradPhi,1));
    assert(size(h,1)+1==size(eta,1));
    assert(size(k,1)==1);
    
    % Set right h.
    h1 = h(1:end-1);
    h2 = h(2:end);
    h = [h1(1);sum([h1,h2],2);h2(end)];
    
    
    nx = size(eta,1);
    mainDiag = 1/k*ones(nx,1); 
    subDiag = -gradPhi(:,1)./h.*ones(nx,1);
    supDiag = -subDiag;
    supDiag = [1/h(1);supDiag(1:end-1)];
    subDiag = [subDiag(2:end);-1/h(end)];
    
    A = spdiags([subDiag,mainDiag,supDiag],[-1,0,1],nx,nx);

    rhs = 1/k*eta + gradPhi(:,2);
    
    eta = A\rhs;
end