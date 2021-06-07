Coor=[0,5,10]; % Initial Coordinate Point: Top of top circle of points
for i=1:19                                                           %
    theta=(pi/180)*(360-18*i);                                       %
    RotMat=[cos(theta),-sin(theta),0;sin(theta),cos(theta),0;0,0,1]; % Gather values for
    Rotated=RotMat*(Coor(1,:)');                                     % the top ring
    Rotated=Rotated';                                                % of coordinates
    Coor=[Coor;Rotated];                                             %
end
ShiftMatrix=zeros(20,3);
for i=1:20
    ShiftMatrix(i,:)=[0,0,-1];
end
for i=1:20                                         %
    CoorShift=Coor(1+(i-1)*20:i*20,:)+ShiftMatrix; % Now for the rest of the grid, in blocks of 20 points
    Coor=[Coor;CoorShift];                         %
end                                                %
CellConc=zeros(400,7);                             %
for i=1:400                                        % Set the initial concentrations of variables
    CellConc(i,:)=[1,0,0,1,0.98+0.04*rand(1),1,1]; % in each of the 400 cells
end                                                %
NoWalls=1;
WallData(1,:)=[0,0,-1,0,1];
for i=1:18                                  %
    WallData=[WallData;NoWalls,i,-1,i,i+1]; % Walls along
    NoWalls=NoWalls+1;                      % the top
end                                         %
WallData=[WallData;NoWalls,19,-1,19,0];
NoWalls=NoWalls+1;
for i=20:38                                     %
    WallData=[WallData;NoWalls,i-20,i,i,i+1]; % Second horizontal
    NoWalls=NoWalls+1;                          % set of walls
end                                             %
WallData=[WallData;NoWalls,19,39,39,20];
NoWalls=NoWalls+1;
ShiftMatrix=zeros(20,5);
for i=1:20
    ShiftMatrix(i,:)=[20,20,20,20,20];
end
for i=1:18                                                    %
    WallDataShift=WallData(NoWalls-19:NoWalls,:)+ShiftMatrix; % Shift the previous
    WallData=[WallData;WallDataShift];                        % set of wall data
    NoWalls=NoWalls+20;                                       %
end                                                           %
for i=400:418                                  %
    WallData=[WallData;NoWalls,i-20,-1,i,i+1]; % Walls down
    NoWalls=NoWalls+1;                         % the RHS
end                                            %
WallData=[WallData;NoWalls,399,-1,419,400];
NoWalls=NoWalls+1;
for i=1:20                                                 %  
    WallData=[WallData;NoWalls,20*(i-1),20*(i-1)+19,20*(i-1),20*i]; % Walls along
    NoWalls=NoWalls+1;                                     % the joining edge
end                                                        %
for i=1:20                                                             %
    WallData=[WallData;NoWalls,20*(i-1),1+20*(i-1),1+20*(i-1),1+20*i]; % Second horizontal
    NoWalls=NoWalls+1;                                                 % set of walls
end                                                                    %
ShiftMatrix=zeros(20,5);
for i=1:20
    ShiftMatrix(i,:)=[20,1,1,1,1];
end
for i=1:18                                                    %
    WallDataShift=WallData(NoWalls-19:NoWalls,:)+ShiftMatrix; % Shift the previous
    WallData=[WallData;WallDataShift];                        % set of wall data
    NoWalls=NoWalls+20;                                       %
end                                                           %
WallConc=zeros(NoWalls,7);                     %
for i=1:NoWalls                                % Set the initial concentrations
    WallConc(i,:)=[norm([Coor(WallData(i,4)+1,:)-Coor(WallData(i,5)+1,:)]),0.1,0.1,0.1,0.1,0.1,0.1]; % of variables in each of the walls
end                                            %