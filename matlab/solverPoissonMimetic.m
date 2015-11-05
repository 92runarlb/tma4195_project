function phi = solverPoissonMimetic(G, h, eta, etat);
    DIM = G.cartDims;
    nx = DIM(1);
    ny = DIM(2);
    nc = G.cells.num;
    nf = G.faces.num;
    half_faces = G.cells.faces(:, 1);
    nhf = numel(half_faces);

    % Identify left and rigth faces, to apply the Dirichlet boundary conditions.
    right_faces = (1 : G.faces.num)';
    right_faces = right_faces(G.faces.centroids(:, 1) == 1);
    right_p = zeros(numel(right_faces), 1);

    left_faces = (1 : G.faces.num)';
    left_faces = left_faces(G.faces.centroids(:, 1) == 0);
    left_p = zeros(numel(left_faces), 1);

    dirich_faces = [right_faces;left_faces]; 

    dirich_p_rigth = right_p;
    dirich_p_left = left_p;

    % Add source (roughly) in the middle
    source_rhs = zeros(nc, 1);
    source_cell = floor(nx/2 + nx*ny/2);
    source_rhs(source_cell) = 0/G.cells.volumes(source_cell);

    % We identify the flux unknown corresponding to half faces where the Dirichlet
    % condition is applied.
    is_dirich_faces = false(nf, 1);
    is_dirich_faces(dirich_faces) = true;
    is_dirich_half_faces = is_dirich_faces(half_faces);
    is_ext_faces  = (G.faces.neighbors(:, 1) == 0) | (G.faces.neighbors(:, 2) == 0) ;
    is_int_faces = ~is_ext_faces;
    nif = nnz(is_int_faces);
    is_neuman_faces = is_ext_faces & ~is_dirich_faces;
    is_neuman_half_faces = is_neuman_faces(half_faces);
    neuman_half_faces = (1 : nhf)';
    neuman_half_faces = neuman_half_faces(is_neuman_half_faces);
    nnhf = nnz(is_neuman_half_faces);

    neuman_rhs = zeros(nnhf,1);
    neuman_rhs(G.faces.centroids(G.cells.faces(neuman_half_faces), 2) == 1) = etat;
    neuman_rhs(G.faces.centroids(G.cells.faces(neuman_half_faces), 2) == 0) = 0;

    % Deform the grid at top and bottom using the functions h and eta.
    %h = @(x) 1;
    %epsilon = 5e-1;
    %eta = @(x) (1/(sqrt(pi)*epsilon)*exp(-(x-0.5).^2/(epsilon^2)));
    x = G.nodes.coords;
    xx = zeros(size(x));
    x1 = x(:, 1);
    xx(:, 1) = x1;
    xx(:, 2) = -h(x1) + (eta + h(x1)).*x(:, 2);
    G.nodes.coords = xx;

    % We compute the mimetic scalar product
    rock.perm = ones(G.cells.num, 1); % this is because MRST is a code for geosciences...
    s   = computeMimeticIP(G, rock);
    BI = s.BI;

    % We assemble the matrices C and D.
    cellNo = rldecode(1:G.cells.num, diff(G.cells.facePos), 2) .';
    C = - sparse(1:numel(cellNo), cellNo, 1);
    D = sparse(1:numel(cellNo), double(half_faces), 1, numel(cellNo), G.faces.num);


    % We construct the right-hand side corresponding to the source term coming from the
    % Dirichlet boundary conditons.
    dirich_pii_rhs = zeros(nf, 1); % called pii instead of pi
    dirich_pii_rhs(right_faces) = dirich_p_rigth;
    dirich_pii_rhs(left_faces) = dirich_p_left;
    dirich_rhs = -D*dirich_pii_rhs; 

    % Reduce the system to the unknown variables.
    D = D(:, is_int_faces);
    N = - sparse(neuman_half_faces, 1 : nnhf , 1, nhf, nnhf);

    % Assemble the source term and compute the right-hand side.
    source_rhs = [zeros(nhf, 1); source_rhs; zeros(nif, 1); zeros(nnhf, 1)];

    % neumann_rhs(bottom_faces) = neumann_bottom_rhs;

    neuman_rhs = [zeros(nhf+nc+nif,1); neuman_rhs];
    dirich_rhs = [dirich_rhs; zeros(nc + nif + nnhf, 1)];

    rhs = source_rhs + dirich_rhs + neuman_rhs;

    % Schur reduction

    R = [[-C'; -D'; -N']*BI, eye(nc + nif + nnhf)];
    A = [[C, D, N]; zeros(nc + nif + nnhf)];
    Q = R*A; 
    rhs = R*rhs;

    % Solve the system.
    sol = Q\rhs;

    % Recover the pressure.
    phi = sol([true(nc, 1); false(nif, 1); false(nnhf, 1)]);

    % Plot the pressure.
    figure(1); clf;
    plotCellData(G, phi);
    colorbar;
end
