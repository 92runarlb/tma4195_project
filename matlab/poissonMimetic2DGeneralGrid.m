function [phi, gradPhi_top, top_cells G] = poissonMimetic2DGeneralGrid(G, h, phi_top, gridLimits, top_faces, bottom_faces, right_faces, left_faces);
%--------------------------------------------------------------------------
%
%   Solves the Laplace equation \Delta \phi = 0 on a grid G.
%
%   Input:
%
%   G:              MRST grid.
%   h, eta:         Functionspecifyting the bottom: -h(x) <= z <= eta(x).
%   phi_top:        Dirichlet condition top boundary; z = eta(x).
%   gridLimits:     Limits of original grid, to identify boundary faces.
%   boundaryCond:   Boundary conditions:
%                   1: Homogeneous Neumann on bottom, homogenepus
%                      Dirichlet on the remaining boundary.
%                   2: Homogeneous Dirichlet on top, homogeneous neumann on
%                      the remaining boundary.
%
%   Output:
%
%   phi:            function value on the cell centroids, taken as the
%                   average over the entire cell.
%   gradPhiSqr:     |\nabla \phi|^2 in all top face centroids.
%   G:              Altered grid.
%
%--------------------------------------------------------------------------
    
    %%  Get grid data

    nc = G.cells.num;
    nf = G.faces.num;
    
    
    %   Half-faces and map from half-faces to cells.
    half_faces = G.cells.faces(:, 1);
    half_face_to_cell = rldecode(1 : G.cells.num, diff(G.cells.facePos), 2)';
    nhf = numel(half_faces);
    

    % Dirichelt top, Neumann bottom, left and right.

    dirich_faces = top_faces;
    dirich_phi = phi_top;

    % Find Neumann faces, set BC's
    bottom_neu = zeros(numel(bottom_faces), 1);
    right_neu = zeros(numel(right_faces),1);
    left_neu = zeros(numel(left_faces),1);

    neumann_faces = [left_faces; right_faces; bottom_faces];
    neumann_Gphi = zeros(nf,1);
    neumann_Gphi(neumann_faces) = [left_neu; right_neu; bottom_neu];
    


    %%  Get data from grid beforE deformation

    %   Find top cells
    top_cells = sum(G.faces.neighbors(top_faces,:),2);
    
    %   Add source (roughly) in the middle
    source_rhs = zeros(nc, 1);
    source_cell = floor(nc/2);
    source_rhs(source_cell) = 0/G.cells.volumes(source_cell); 

    %   Find internal faces
    is_ext_faces  = (G.faces.neighbors(:, 1) == 0) | (G.faces.neighbors(:, 2) == 0) ;
    is_int_faces = ~is_ext_faces;
    nif = nnz(is_int_faces);
    
    %   Find Neumann half faces
    is_top_cell_faces = false(nf,1);
    is_top_cell_faces(neumann_faces) = true;
    is_top_cell_half_faces = is_top_cell_faces(half_faces);
    neuman_half_faces = (1 : nhf)';
    neuman_half_faces = neuman_half_faces(is_top_cell_half_faces);
    nnhf = nnz(is_top_cell_half_faces);

    %%  Set Neumann rhs
    
    neumann_rhs = neumann_Gphi(half_faces(neuman_half_faces));
   
    %% Compute Mimetic scalar product
    
    rock.perm = ones(G.cells.num, 1);
    s   = computeMimeticIP(G, rock);
    BI = s.BI; % Inverse of B
    
    %   Assemble the matrices C and D.
    cellNo = rldecode(1:G.cells.num, diff(G.cells.facePos), 2) .';
    C = - sparse(1:numel(cellNo), cellNo, 1);
    D = sparse(1:numel(cellNo), double(half_faces), 1, numel(cellNo), G.faces.num);

    %%  Set Dirichlet rhs

    dirich_pii_rhs = zeros(nf, 1); % called pii instead of pi
    dirich_pii_rhs(dirich_faces) = dirich_phi; 
    dirich_rhs = -D*dirich_pii_rhs; 

    %%  Reduce the system to the unknown variables.
    
    D = D(:, is_int_faces);
    N = - sparse(neuman_half_faces, 1 : nnhf , 1, nhf, nnhf);
    
    %%  Assemble the source term and compute the right-hand side.
    
    source_rhs = [zeros(nhf, 1); source_rhs; zeros(nif, 1); zeros(nnhf, 1)];
    dirich_rhs = [dirich_rhs; zeros(nc + nif + nnhf, 1)];
    neumann_rhs = [zeros(nhf+nc+nif, 1); neumann_rhs];
    rhs = source_rhs + dirich_rhs + neumann_rhs;

    %%  Schur reduction

    R = [[-C'; -D'; -N']*BI, eye(nc + nif + nnhf)];
    A = [[C, D, N]; zeros(nc + nif + nnhf)];
    Q = R*A; 
    rhs = R*rhs;

    %%  Solve the system.
    
    sol = Q\rhs;
    phi = sol([true(nc, 1); false(nif, 1); false(nnhf, 1)]);
    pii = sol([false(nc, 1); true(nif, 1); false(nnhf, 1)]);
    phi_neu = sol([false(nc, 1); false(nif, 1); true(nnhf, 1)]);
    
    %   Calculate half-face fluxes
    gradPhiN = BI*(dirich_rhs(1:nhf) - C*phi - D*pii - N*phi_neu);
   
    %%  Find |\nabla \phi|^2 in all top face centroid
    
    % Finde face areas in deformed grid
    face_areas = G.faces.areas;
    gradPhi_top = zeros(size(phi_top,1),2);
    for i = 1:numel(top_cells)
        
        c = top_cells(i);
        facePos = G.cells.facePos(c):G.cells.facePos(c+1)-1;
        faces = G.cells.faces(facePos,1);
        normals = G.faces.normals(faces(:,1),:);
        
        nf = numel(faces);
        gradPhifk = zeros(nf,1);
        for f = 1:numel(faces);
            gradPhifk(f) = gradPhiN(and(half_faces==faces(f), half_face_to_cell== c));
        end
        
        gradPhifk = gradPhifk./face_areas(faces);
        normals = normals./repmat(sqrt(sum(normals.^2,2)),1,2);
        normals = normals.*(repmat(-1+2*(G.faces.neighbors(faces,1) == c),1,2));
        n1 = [0,1]; n2 = [1,0];
        for f = 1:nf
            gradPhi_top(i,:) = gradPhi_top(i,:) + gradPhifk(f)*dot(normals(f,:),n1)*n1 + gradPhifk(f)*dot(normals(f,:),n2)*n2;
        end
    end

end