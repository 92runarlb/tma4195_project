function [phi, gradPhiSqr] = poissonMimetic(G, h, eta, etat);
    nc = G.cells.num;
    nf = G.faces.num;
    half_faces = G.cells.faces(:, 1);
    nhf = numel(half_faces);
    
    % Find dirichlet faces. They are on the left and right boundary
    right_faces = (1 : G.faces.num)';
    right_faces = right_faces(G.faces.centroids(:, 1) == 1);
    right_phi = 0.5*(G.faces.centroids(right_faces,2).^2-1);
    
    left_faces = (1 : G.faces.num)';    
    left_faces = left_faces(G.faces.centroids(:, 1) == 0);
    
    left_phi = 0.5*(G.faces.centroids(right_faces,2).^2);
    
    dirich_faces = [right_faces;left_faces];
    dirich_phi = [right_phi;left_phi];

    % Identify bottom and top faces, to apply the Neuman boundary conditions.
    top_faces = (1 : G.faces.num)';
    top_faces = top_faces(G.faces.centroids(:, 2) == 1);
    top_neu = etat;

    bottom_faces = (1 : G.faces.num)';
    bottom_faces = bottom_faces(G.faces.centroids(:, 2) == 0);
    bottom_neu = zeros(numel(bottom_faces), 1);

    neumann_faces = [top_faces; bottom_faces];
    neumann_Gphi = zeros(nf,1);
    neumann_Gphi(neumann_faces) = [top_neu; bottom_neu];
    
    neuman_Gphi_index = 1:size(neumann_Gphi,1);
    
    
    % Add source (roughly) in the middle
    source_rhs = zeros(nc, 1);
    source_cell = floor(nc/2);
    source_rhs(source_cell) = 0/G.cells.volumes(source_cell); 

    %

    % We identify the flux unknown correspoinding to half faces where the 
    % Dirichlet conditions are applied.
    %is_dirich_faces = false(nf, 1);
    %is_dirich_faces(dirich_faces) = true;
    %is_dirich_half_faces = is_dirich_faces(half_faces);
    is_ext_faces  = (G.faces.neighbors(:, 1) == 0) | (G.faces.neighbors(:, 2) == 0) ;
    is_int_faces = ~is_ext_faces;
    nif = nnz(is_int_faces);
    
    % Find Neumann half faces.
    is_neumann_faces = false(nf,1);
    is_neumann_faces(neumann_faces) = true;
    %is_neumann_faces = is_ext_faces & ~is_dirich_faces
    is_neumann_half_faces = is_neumann_faces(half_faces);
    neuman_half_faces = (1 : nhf)';
    neuman_half_faces = neuman_half_faces(is_neumann_half_faces);
    nnhf = nnz(is_neumann_half_faces);
    
    %topHalfFaces = neuman_half_faces(G.faces.centroids(G.cells.faces(neuman_half_faces), 2) == 1);
    %bottomHalfFaces = neuman_half_faces(G.faces.centroids(G.cells.faces(neuman_half_faces), 2) == 0);
       
    neumann_rhs = neumann_Gphi(half_faces(neuman_half_faces));

    
    % Deform the grid.
    x = G.nodes.coords;
    xx = zeros(size(x));
    x1 = x(:, 1);
    xx(:, 1) = x1;
    xx(:, 2) = -h(x1) + (eta + h(x1)).*x(:, 2);
    G.nodes.coords = xx;
   
    % We compute the mimetic scalar product
    rock.perm = ones(G.cells.num, 1); % this is because MRST is a code for geosciences...
    s   = computeMimeticIP(G, rock);
    BI = s.BI; % Inverse of B (for B as defined on slide)
    
    % We assemble the matrices C and D.
    cellNo = rldecode(1:G.cells.num, diff(G.cells.facePos), 2) .';
    C = - sparse(1:numel(cellNo), cellNo, 1);
    D = sparse(1:numel(cellNo), double(half_faces), 1, numel(cellNo), G.faces.num);

    % We construct the right-hand side corresponding to the source term coming from the
    % Dirichlet boundary conditons.
    dirich_pii_rhs = zeros(nf, 1); % called pii instead of pi
    dirich_pii_rhs(dirich_faces) = dirich_phi; 
    dirich_rhs = -D*dirich_pii_rhs; 

    % Reduce the system to the unknown variables.
    D = D(:, is_int_faces);
    N = - sparse(neuman_half_faces, 1 : nnhf , 1, nhf, nnhf);
    
    % Assemble the source term and compute the right-hand side.
    source_rhs = [zeros(nhf, 1); source_rhs; zeros(nif, 1); zeros(nnhf, 1)];
    dirich_rhs = [dirich_rhs; zeros(nc + nif + nnhf, 1)];
    neumann_rhs = [zeros(nhf+nc+nif, 1); neumann_rhs];
    rhs = source_rhs + dirich_rhs + 1/40*neumann_rhs;

    % Schur reduction

    R = [[-C'; -D'; -N']*BI, eye(nc + nif + nnhf)];
    A = [[C, D, N]; zeros(nc + nif + nnhf)];
    Q = R*A; 
    rhs = R*rhs;

    % Solve the system.
    sol = Q\rhs;
    phi = sol([true(nc, 1); false(nif, 1); false(nnhf, 1)]);
    pii = sol([false(nc, 1); true(nif, 1); false(nnhf, 1)]);
    phi_neu = sol([false(nc, 1); false(nif, 1); true(nnhf, 1)]);
    gradPhiN = BI*(dirich_rhs(1:nhf) - C*phi - D*pii - N*phi_neu);
    
    gradPhiSqr = gradPhiN;
    
    plotCellData(G, phi);
    colorbar;

end