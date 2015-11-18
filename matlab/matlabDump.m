    
    %    G.cells.centroids(G.cells.faces)
% Dirichlet conditions are applied.
    %is_dirich_faces = false(nf, 1);
    %is_dirich_faces(dirich_faces) = true;
    %is_dirich_half_faces = is_dirich_faces(half_faces);
    
    
        
    %topHalfFaces = neuman_half_faces(G.faces.centroids(G.cells.faces(neuman_half_faces), 2) == 1);
    %bottomHalfFaces = neuman_half_faces(G.faces.centroids(G.cells.faces(neuman_half_faces), 2) == 0);
    
    %is_neumann_faces = is_ext_faces & ~is_dirich_faces
    %is_top_cell_half_faces = zeros(nhf,1);
    
%         for m = 1:size(normals,1)
%             n1temp = normals(m,:)/norm(normals(m,:),2);
%             if is_int_faces(faces(m));
%                 n1 = n1temp;
%                 n1index = m;
%                 break
%             end            
%         end
%         minn2 = 1;
%         for n = 1:size(normals,1)
%             n2temp = normals(n,:)/norm(normals(n,:),2);
%             if minn2>dot(n1, n2temp);
%                 n2 = n2temp;
%                 minn2 = dot(n1,n2temp);
%                 n2index = n;
%             end
%         end

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

%    is_top_cell_half_faces = logical(is_top_cell_half_faces);
    %gradPhiSqr = gradPhiN(is_top_cell_half_faces);
    
    
    %plotGrid(G, c,'facecolor','r')
    %plotCellData(G, phi);
    %colorbar;
    
    

    
    % Calculate the neuman rhs
    dl = face_areas(neumann_faces);
    deta = sqrt((dl./dx).^2 -1);
    nvecScale = sqrt(sum([deta,ones(size(deta,1),1)].^2,2));
    neumann_Gphi(neumann_faces) = neumann_Gphi(neumann_faces).*dl./nvecScale;
        dx = G.faces.areas(neumann_faces);
        
            
    %   Find face areas, to be used in the gradient approximation of eta
    