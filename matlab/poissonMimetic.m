function [phi, gradPhiSqr, top_cells,G] = poissonMimetic(G, h, eta, etat);
    nc = G.cells.num;
    nf = G.faces.num;
    half_faces = G.cells.faces(:, 1);
    nhf = numel(half_faces);
    
    % Find dirichlet faces. They are on the left and right boundary
    right_faces = (1 : G.faces.num)';
    right_faces = right_faces(G.faces.centroids(:, 1) == 1);
    right_phi = 0*0.5*(G.faces.centroids(right_faces,2).^2-1);
    
    left_faces = (1 : G.faces.num)';    
    left_faces = left_faces(G.faces.centroids(:, 1) == 0);
    
    left_phi = 0*0.5*(G.faces.centroids(right_faces,2).^2);
    
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
    
    
    % Set the change in dx
    dx = G.faces.areas(neumann_faces);
    top_cells = sum(G.faces.neighbors(top_faces,:),2);
    
    % Add source (roughly) in the middle
    source_rhs = zeros(nc, 1);
    source_cell = floor(nc/2);
    source_rhs(source_cell) = 0/G.cells.volumes(source_cell); 

    %    G.cells.centroids(G.cells.faces)

    % We identify the flux unknown correspoinding to half faces where the 
    % Dirichlet conditions are applied.
    %is_dirich_faces = false(nf, 1);
    %is_dirich_faces(dirich_faces) = true;
    %is_dirich_half_faces = is_dirich_faces(half_faces);
    is_ext_faces  = (G.faces.neighbors(:, 1) == 0) | (G.faces.neighbors(:, 2) == 0) ;
    is_int_faces = ~is_ext_faces;
    nif = nnz(is_int_faces);
    
    % Find Neumann half faces.
    is_top_cell_faces = false(nf,1);
    is_top_cell_faces(neumann_faces) = true;
    %is_neumann_faces = is_ext_faces & ~is_dirich_faces
    is_top_cell_half_faces = is_top_cell_faces(half_faces);
    neuman_half_faces = (1 : nhf)';
    neuman_half_faces = neuman_half_faces(is_top_cell_half_faces);
    nnhf = nnz(is_top_cell_half_faces);
    
    %topHalfFaces = neuman_half_faces(G.faces.centroids(G.cells.faces(neuman_half_faces), 2) == 1);
    %bottomHalfFaces = neuman_half_faces(G.faces.centroids(G.cells.faces(neuman_half_faces), 2) == 0);


    % Deform the grid.
    x = G.nodes.coords;
    xx = zeros(size(x));
    x1 = x(:, 1);
    xx(:, 1) = x1;
    xx(:, 2) = -h(x1) + (eta + h(x1)).*x(:, 2);
    G.nodes.coords = xx;
    G = computeGeometry(G);
    
    % Finde face areas in deformed grid
    face_areas = G.faces.areas;
    % Calculate the neuman rhs
    dl = face_areas(neumann_faces);
    deta = sqrt((dl./dx).^2 -1);
    nvecScale = sqrt(sum([deta,ones(size(deta,1),1)].^2,2));
    neumann_Gphi(neumann_faces) = neumann_Gphi(neumann_faces).*dl./nvecScale;
    nvecScale
        
    neumann_rhs = neumann_Gphi(half_faces(neuman_half_faces));

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
    rhs = source_rhs + dirich_rhs + neumann_rhs;

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
   
    %is_top_cell_half_faces = zeros(nhf,1);

    gradPhiSqr = zeros(size(etat));
    for i = 1:numel(top_cells)
        c = top_cells(i);
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
        
        gradPhifk1 = gradPhiN(half_faces==face1);
        gradPhifk1 = gradPhifk1(1)/face_areas(face1);
        gradPhifk2 = gradPhiN(half_faces==face2);
        gradPhifk2 = gradPhifk2(1)/face_areas(face2);
        
        %gradPhifk1 = gradPhifk(find(find(is_int_faces)==faces1));
        %gradPhifk2 = gradPhifk(find(find(is_int_faces)==faces2));
        gradPhiSqr(i) = gradPhifk1^2 + (1-dot(n1,n2)^2)*gradPhifk2^2;
 
        
%         c = topCells(i);
%         %is_top_faces = false(nf,1);
%         faces = G.cells.faces(G.cells.facePos(c):G.cells.facePos(c+1)-1);
%         %is_top_faces(faces) = true;
%         
%         % Find Neumann half faces.
%         is_top_cell_faces = false(nf,1);
%         is_top_cell_faces(faces) = true;
%         %is_neumann_faces = is_ext_faces & ~is_dirich_faces
%         is_top_cell_half_faces = is_top_cell_half_faces + is_top_cell_faces(half_faces);
%         gradPhiN(is_top_cell_half_faces)
    end
%    is_top_cell_half_faces = logical(is_top_cell_half_faces);
    %gradPhiSqr = gradPhiN(is_top_cell_half_faces);
    
    
    %plotGrid(G, c,'facecolor','r')
    %plotCellData(G, phi);
    %colorbar;

    
end