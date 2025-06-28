% efficient assembly using triplets
I = zeros(nelx*nely*64,1); J = zeros(nelx*nely*64,1); X = zeros(nelx*nely*64,1); 
F = sparse(2*(nely+1)*(nelx+1),1); U = zeros(2*(nely+1)*(nelx+1),1); 
ntriplets = 0; 
for elx = 1:nelx 
   for ely = 1:nely 
     n1 = (nely+1)*(elx-1)+ely; 
     n2 = (nely+1)* elx +ely; 
     edof = [2*n1-1 2*n1 2*n2-1 2*n2 2*n2+1 2*n2+2 2*n1+1 2*n1+2]; 
     xval = x(ely,elx)^penal; 
     for krow = 1:8 
       for kcol = 1:8 
         ntriplets = ntriplets+1; 
         I(ntriplets) = edof(krow); 
         J(ntriplets) = edof(kcol); 
         X(ntriplets) = xval*KE(krow,kcol); 
       end 
     end 
   end 
end 
K = sparse(I,J,X,2*(nelx+1)*(nely+1),2*(nelx+1)*(nely+1));


% plotting displacements
colormap(gray); axis equal; axis tight; axis off;
Deform = [U(1:2:end,1),U(2:2:end,1)];
patch('Faces',Faces,'Vertices',Grid+Deform*0.05,'FaceColor','Flat', ...
    'FaceVertexCData',1-x(:),'EdgeColor','none');
drawnow; clf;

