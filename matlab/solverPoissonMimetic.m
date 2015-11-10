function [phiSurf, gradPhiSurf] = solverPoissonMimetic(G, h, eta, etat);
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
    right_p = 0.5*(G.faces.centroids(right_faces,2).^2-1); %ones(numel(right_faces), 1)

    left_faces = (1 : G.faces.num)';    
    left_faces = left_faces(G.faces.centroids(:, 1) == 0);
    
    left_p = 0.5*(G.faces.centroids(right_faces,2).^2); % ones(numel(right_faces), 1)

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
    
    is_topHalfFaces = G.faces.centroids(G.cells.faces(neuman_half_faces), 2) == 1;
    is_bottomHalfFaces = G.faces.centroids(G.cells.faces(neuman_half_faces), 2) == 0;
    dx = G.faces.areas(G.cells.faces(neuman_half_faces(is_topHalfFaces),1));
    %%Crap 
    topHalfFaces = neuman_half_faces(G.faces.centroids(G.cells.faces(neuman_half_faces), 2) == 1);
    bottomHalfFaces = neuman_half_faces(G.faces.centroids(G.cells.faces(neuman_half_faces), 2) == 0);
    topCells = zeros(size(topHalfFaces));
    
    deta = diff(eta(G.nodes.coords(:,2)==1));
    deta = [deta]./dx;
    nvec = sqrt(sum([deta,ones(size(deta,1),1)].^2,2));
    
    for i = 1:size(topHalfFaces,1)
       topCells(i) = sum(G.faces.neighbors(G.cells.faces(topHalfFaces(i)),:),2); 
    end

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

    G = computeGeometry(G);
    faceLengths = G.faces.areas(G.cells.faces(neuman_half_faces(is_topHalfFaces),1));
    
    neuman_rhs = zeros(nnhf,1);
    neuman_rhs(is_topHalfFaces) = etat./nvec;
    neuman_rhs(is_bottomHalfFaces) = 0;

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

    rhs = (source_rhs + dirich_rhs + neuman_rhs);

    % Schur reduction

    R = [[-C'; -D'; -N']*BI, eye(nc + nif + nnhf)];
    A = [[C, D, N]; zeros(nc + nif + nnhf)];
    Q = R*A; 
    rhs = R*rhs;

    % Solve the system.
    sol = Q\rhs;

    % Recover phi.
    phi = sol([true(nc, 1); false(nif,1); false(nnhf, 1)]);
    pii = sol([false(nc, 1); true(nif,1); false(nnhf, 1)]);
    p_neum = sol([false(nc, 1); false(nif,1); true(nnhf, 1)]);
    gradPhiHalfFace = BI*((C*phi)+(D*pii)+N*p_neum); % Should be - but 
                                                     % for some reason the
                                                     % sign is wrong
    
    plotGrid(G);
    gradPhiSurf = zeros(size(topCells));
    for i = 1:size(topCells,1)
        c = topCells(i);
        facePos = G.cells.facePos(c):G.cells.facePos(c+1)-1;
        faces = G.cells.faces(facePos,:);
        normals = G.faces.normals(faces(:,1),:);
        n1index = find(faces(:,2)==4);
        n1 = normals(n1index,:);
        n1 = n1/norm(n1,2);
        
%         for m = 1:size(normals,1)
%             n1temp = normals(m,:)/norm(normals(m,:),2);
%             if is_int_faces(faces(m));
%                 n1 = n1temp;
%                 n1index = m;
%                 break
%             end            
%         end
        minn2 = 1;
        for n = 1:size(normals,1)
            n2temp = normals(n,:)/norm(normals(n,:),2);
            if minn2>dot(n1, n2temp) && is_int_faces(faces(n));
                n2 = n2temp;
                minn2 = dot(n1,n2temp);
                n2index = n;
            end
        end
        face1 = faces(n1index);
        face2 = faces(n2index);
        
        n1
        n2
        gradPhifk1 = gradPhiHalfFace(half_faces==face1);
        gradPhifk1 = gradPhifk1(1);
        gradPhifk2 = gradPhiHalfFace(half_faces==face2);
        gradPhifk2 = gradPhifk2(1);
        
        %gradPhifk1 = gradPhifk(find(find(is_int_faces)==faces1));
        %gradPhifk2 = gradPhifk(find(find(is_int_faces)==faces2));
        gradPhiSurf(i) = gradPhifk1^2 + (1-dot(n1,n2)^2)*gradPhifk2^2
        plotGrid(G, c,'facecolor','r')
    end
        
    phif = sol([false(nc, 1); false(nif, 1); true(nnhf, 1)]);
    
    phiSurf = phif(is_topHalfFaces);
    
    % Plot the pressure.
    figure(1); clf;
    plotCellData(G, phi);
    colorbar;
end
