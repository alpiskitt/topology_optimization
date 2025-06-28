function  [] = PlotLoadsAndSupports2021(nelx, nely, fixeddofs, F, kspringdofs, PlotNodeNumbers, PlotNodes)
% Input: 
%   nelx:         Number of elements in x-direction
%   nely:         Number of elements in y-direction
%   fixeddofs:    List of fixed dofs
%   F:            Force vector
%   kspringdofs:  Vector with springs stiffnesses, same numbering as the force vector.
%   PlotNodeNumbers: 'yes' or 'no'
%   PlotNodes: 'yes' or 'no'
%
% Output:
%   A plot with the loads, the supports and the node numbers
%
% Modified 2021 by Casper Schousboe Andreasen, csan@mek.dtu.dk
% Original version by:
% chrlund@mek.dtu.dk or chrlundg@gmail.com

% Initialize the figure that will be used
figure(1001); clf;
hold on;

% Some settings
NodeNumberPlacementX = 0.25/2;
NodeNumberPlacementY = 0.25;
ScalingTriangle = 0.4;
ScalingRollers = 0.1;
LineWidthInRollers = 1;
PlacementOfCirclesInRollers = ScalingTriangle/3;
HeadSizeOfForceVectors = 10;
ForceVectorsLineWidth = 2;
LineWidthSpring = 1;
PointsInCircleInRollers = 0:pi/50:2*pi;

% Use the edofMat structure as in the top88.m code
nodenrs = reshape(1:(1+nelx)*(1+nely),1+nely,1+nelx);
edofVec = reshape(2*nodenrs(1:end-1,1:end-1)+1,nelx*nely,1);
edofMat = repmat(edofVec,1,8)+repmat([0 1 2*nely+[2 3 0 1] -2 -1],nelx*nely,1);

% Display a warning if many nodenumbers must be plotted (slow)
if nelx*nely > 200
    disp('Plotting the loads and supports for a large numbers of nodes may take a long time')
    if strcmp(PlotNodeNumbers,'yes')
        disp('"PlotNodeNumbers = no" may provide a speed up')
    end
end

%% Initialize plotting
Faces = ((edofMat(:,[1 3 5 7])-1)/2)+1;
[XGrid,YGrid]=meshgrid(0:nelx,nely:-1:0);
Grid = [XGrid(:),YGrid(:)];
nn=(nelx+1)*(nely+1); 

% Plot nodes and node numbering
if strcmp(PlotNodes,'yes')
    plot(XGrid(:),YGrid(:),'xk');
end
if strcmp(PlotNodeNumbers,'yes')
    text(XGrid(:)+NodeNumberPlacementX,YGrid(:)+NodeNumberPlacementY,num2cell(1:nn));
end
% Plot the elements
patch('Faces',Faces,'Vertices',Grid,'FaceColor','none', ...
    'FaceVertexCData',1,'EdgeColor',[1 1 1]*.5);

%% Plot load vectors
% Matrix with loads:
% node, dof, magnitude, x-coordinate, ycoordinate
Fnodes=floor(((find(F~=0))+1)/2);
F0 = [Fnodes find(F~=0) full(F((F~=0))) XGrid(Fnodes) YGrid(Fnodes)];
for i = 1:nnz(F)
    if mod(F0(i,2),2) == 1
        quiver(F0(i,4), F0(i,5), F0(i,3), 0, 'LineWidth', ForceVectorsLineWidth, 'Color','blue','MaxHeadSize', HeadSizeOfForceVectors)
    elseif mod(F0(i,2), 2) == 0
        quiver(F0(i,4), F0(i,5), 0, F0(i,3), 'LineWidth', ForceVectorsLineWidth, 'Color','blue','MaxHeadSize', HeadSizeOfForceVectors)
    end
end
 
%% Plot supports
% Matrix with supports
% node number, dof, x-coordinate, y-coordinate
Bnodes = floor((fixeddofs+1)/2)';
B0 = [Bnodes  fixeddofs'  XGrid(Bnodes) YGrid(Bnodes)];
 
for i = 1:size(B0,1)
    xcor = B0(i,3);
    ycor = B0(i,4);
     
    if mod(B0(i,2),2) == 1 % x-constraints
        % Plot triangles in rollers
        plot([xcor xcor-ScalingTriangle],[ycor ycor-ScalingTriangle],'color','r','LineWidth', LineWidthInRollers);
        plot([xcor-ScalingTriangle xcor-ScalingTriangle],[ycor+ScalingTriangle ycor-ScalingTriangle],'color','r','LineWidth', LineWidthInRollers);
        plot([xcor xcor-ScalingTriangle],[ycor ycor+ScalingTriangle],'color','r','LineWidth', LineWidthInRollers);
         
        % Plot circles in rollers
        xunit = ScalingRollers * cos(PointsInCircleInRollers) + xcor-ScalingTriangle-PlacementOfCirclesInRollers;
        yunit = ScalingRollers * sin(PointsInCircleInRollers) + ycor-ScalingTriangle+PlacementOfCirclesInRollers;
        plot(xunit, yunit, 'color', 'r','LineWidth', LineWidthInRollers);
         
        xunit = ScalingRollers * cos(PointsInCircleInRollers) + xcor-ScalingTriangle-PlacementOfCirclesInRollers;
        yunit = ScalingRollers * sin(PointsInCircleInRollers) + ycor+ScalingTriangle-PlacementOfCirclesInRollers;
        plot(xunit, yunit, 'color', 'r','LineWidth', LineWidthInRollers);
         
        xunit = ScalingRollers * cos(PointsInCircleInRollers) + xcor;
        yunit = ScalingRollers * sin(PointsInCircleInRollers) + ycor;
        plot(xunit, yunit, 'color', 'r','LineWidth', LineWidthInRollers);
    elseif mod(B0(i,2),2) == 0 % y-constraints
        % Plot triangles in rollers
        plot([xcor xcor+ScalingTriangle],[ycor ycor-ScalingTriangle],'color','r','LineWidth', LineWidthInRollers);
        plot([xcor-ScalingTriangle xcor+ScalingTriangle],[ycor-ScalingTriangle ycor-ScalingTriangle],'color','r','LineWidth', LineWidthInRollers);
        plot([xcor xcor-ScalingTriangle],[ycor ycor-ScalingTriangle],'color','r','LineWidth',LineWidthInRollers);
         
        % Plot circles in rollers
        xunit = ScalingRollers * cos(PointsInCircleInRollers) + xcor+ScalingTriangle-PlacementOfCirclesInRollers;
        yunit = ScalingRollers * sin(PointsInCircleInRollers) + ycor-ScalingTriangle-PlacementOfCirclesInRollers;
        plot(xunit, yunit, 'color', 'r','LineWidth', LineWidthInRollers);
         
        xunit = ScalingRollers * cos(PointsInCircleInRollers) + xcor-ScalingTriangle+PlacementOfCirclesInRollers;
        yunit = ScalingRollers * sin(PointsInCircleInRollers) + ycor-ScalingTriangle-PlacementOfCirclesInRollers;
        plot(xunit, yunit, 'color', 'r','LineWidth', LineWidthInRollers);
         
        xunit = ScalingRollers * cos(PointsInCircleInRollers) + xcor;
        yunit = ScalingRollers * sin(PointsInCircleInRollers) + ycor;
        plot(xunit, yunit, 'color', 'r','LineWidth', LineWidthInRollers);
    else
        error('Not showing any BC')
    end
     
end
 
 
%% Plot springs 
Knodes = floor(((find(kspringdofs~=0))+1)/2);
K0 = [Knodes find(kspringdofs~=0) full(kspringdofs((kspringdofs~=0))) XGrid(Knodes) YGrid(Knodes)];
for i = 1:nnz(kspringdofs) 
      if mod(K0(i,2),2) == 0 % y-constraints
        plot([K0(i,4) K0(i,4)], [K0(i,5) K0(i,5)-0.2], 'color', 'g','LineWidth', LineWidthSpring)
        plot([K0(i,4) K0(i,4)+0.2], [K0(i,5)-0.2 K0(i,5)-0.25], 'color', 'g','LineWidth', LineWidthSpring)
        plot([K0(i,4)+0.2 K0(i,4)-0.2], [K0(i,5)-0.25 K0(i,5)-0.35], 'color', 'g','LineWidth', LineWidthSpring)
        plot([K0(i,4)+0.2 K0(i,4)-0.2], [K0(i,5)-0.45 K0(i,5)-0.35], 'color', 'g','LineWidth', LineWidthSpring)
        plot([K0(i,4)+0.2 K0(i,4)], [K0(i,5)-0.45 K0(i,5)-0.5], 'color', 'g','LineWidth', LineWidthSpring)
        plot([K0(i,4) K0(i,4)], [K0(i,5)-0.5 K0(i,5)-0.7],'color', 'g','LineWidth', LineWidthSpring)
        plot([K0(i,4)-0.2 K0(i,4)+0.2], [K0(i,5)-0.7 K0(i,5)-0.7],'color', 'g','LineWidth', LineWidthSpring)        
        plot([K0(i,4)-0.25 K0(i,4)-0.2], [K0(i,5)-0.8 K0(i,5)-0.7],'color', 'g','LineWidth', LineWidthSpring)
        plot([K0(i,4)-0.15 K0(i,4)-0.1], [K0(i,5)-0.8 K0(i,5)-0.7],'color', 'g','LineWidth', LineWidthSpring)        
        plot([K0(i,4)-0.05 K0(i,4)], [K0(i,5)-0.8 K0(i,5)-0.7],'color', 'g','LineWidth', LineWidthSpring)                       
        plot([K0(i,4)+0.15 K0(i,4)+0.2], [K0(i,5)-0.8 K0(i,5)-0.7],'color', 'g','LineWidth', LineWidthSpring)        
        plot([K0(i,4)+0.05 K0(i,4)+0.1], [K0(i,5)-0.8 K0(i,5)-0.7],'color', 'g','LineWidth', LineWidthSpring)        
        plot([K0(i,4) K0(i,4)], [K0(i,5) K0(i,5)], '.g','LineWidth', LineWidthSpring)
         
      elseif mod(K0(i,2),2) == 1 % x-constraints
        plot([K0(i,4) K0(i,4)-0.2], [K0(i,5) K0(i,5)], 'color', 'g','LineWidth', LineWidthSpring)
        plot([K0(i,4)-0.2 K0(i,4)-0.25], [K0(i,5) K0(i,5)+0.2], 'color', 'g','LineWidth', LineWidthSpring)
        plot([K0(i,4)-0.25 K0(i,4)-0.35], [K0(i,5)+0.2 K0(i,5)-0.2], 'color', 'g','LineWidth', LineWidthSpring)
        plot([K0(i,4)-0.45 K0(i,4)-0.35], [K0(i,5)+0.2 K0(i,5)-0.2], 'color', 'g','LineWidth', LineWidthSpring)
        plot([K0(i,4)-0.45 K0(i,4)-0.5], [K0(i,5)+0.2 K0(i,5)], 'color', 'g','LineWidth', LineWidthSpring)
        plot([K0(i,4)-0.5 K0(i,4)-0.7], [K0(i,5) K0(i,5)],'color', 'g','LineWidth', LineWidthSpring)
        plot([K0(i,4)-0.7 K0(i,4)-0.7], [K0(i,5)-0.2 K0(i,5)+0.2],'color', 'g','LineWidth', LineWidthSpring)
                 
        plot([K0(i,4)-0.8 K0(i,4)-0.7], [K0(i,5)-0.25 K0(i,5)-0.2],'color', 'g','LineWidth', LineWidthSpring)
        plot([K0(i,4)-0.8 K0(i,4)-0.7], [K0(i,5)-0.15 K0(i,5)-0.1],'color', 'g','LineWidth', LineWidthSpring)        
        plot([K0(i,4)-0.8 K0(i,4)-0.7], [K0(i,5)-0.05 K0(i,5)],'color', 'g','LineWidth', LineWidthSpring)                       
        plot([K0(i,4)-0.8 K0(i,4)-0.7], [K0(i,5)+0.15 K0(i,5)+0.2],'color', 'g','LineWidth', LineWidthSpring)        
        plot([K0(i,4)-0.8 K0(i,4)-0.7], [K0(i,5)+0.05 K0(i,5)+0.1],'color', 'g','LineWidth', LineWidthSpring)        
        plot([K0(i,4) K0(i,4)], [K0(i,5) K0(i,5)]', '.g','LineWidth', LineWidthSpring)
      else
          error('No spring');
      end
 
end
MaxForceSize = max(abs(F0(:,3)));
xlim([-MaxForceSize nelx+MaxForceSize]);
ylim([-MaxForceSize nely+MaxForceSize]);


% Set figure properties
xlabel('x-coordinate')
ylabel('y-coordinate')
 


axis equal 
hold off
