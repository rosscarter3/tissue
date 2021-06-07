% MatLab script for generating tissue geometry file for 2d hypocotyl cross section.
% By Behruz Bozorg, July 2017

%change cell ids to vectors
clear all
nCell=110;
nWall=nCell*6;
nVer=nCell*6;
nVar=30;

connectivity=zeros(nWall,5);
for i=1:6
    connectivity(i,1)=i-1;
end
for i=1:6
    connectivity(i,2)=0;
end
for i=1:6
    connectivity(i,3)=-1;
end
for i=1:6
    connectivity(i,4)=i-1;
end
for i=1:5
    connectivity(i,5)=i;
end

if nCell>1
  for i=2:nCell
     for j=1:6
        connectivity(j+(i-1)*6,1)=connectivity(j+(i-2)*6,1)+6;
        connectivity(j+(i-1)*6,2)=connectivity(j+(i-2)*6,2)+1;
        connectivity(j+(i-1)*6,3)=connectivity(j,3);
        connectivity(j+(i-1)*6,4)=connectivity(j+(i-2)*6,4)+6;
        connectivity(j+(i-1)*6,5)=connectivity(j+(i-2)*6,5)+6;
     end
 end
end

connectivity
coor=zeros(nVer,3);
coor(1,:)=[0 0 0];
coor(2,:)=[1 0 0];
coor(3,:)=[1 3 0];
coor(4,:)=[1 6 0];
coor(5,:)=[0 6 0];
coor(6,:)=[0 3 0];

for i=2:10
coor(1+6*(i-1),:)=[0 0 0]+[0 6 0]*(i-1);
coor(2+6*(i-1),:)=[1 0 0]+[0 6 0]*(i-1);
coor(3+6*(i-1),:)=[1 3 0]+[0 6 0]*(i-1);
coor(4+6*(i-1),:)=[1 6 0]+[0 6 0]*(i-1);
coor(5+6*(i-1),:)=[0 6 0]+[0 6 0]*(i-1);
coor(6+6*(i-1),:)=[0 3 0]+[0 6 0]*(i-1);
end

for j=2:11
    if mod(j,2)==0
          for i=1:60
              coor(((j-1)*60)+i,:)=coor(i,:)+[(j-1)*1 3 0];
          end
    else
        for i=1:60
              coor(((j-1)*60)+i,:)=coor(i,:)+[(j-1)*1 0 0];
          end
    end

end


coor;
walls=zeros(nWall,3);
for i=1:nWall
walls(i,1)=1;
end
walls;
cells=zeros(nCell,nVar);

for i=1:nCell
cells(i,1)=1;
end

for i=1:nCell
    if (mod(i,10)==0 || mod(i,10)==1)
      cells(i,29)=-1;
    end
end

for i=1:nCell
    if (i<=10 || i>100)
      cells(i,30)=1;
    end
    if ((i>10 && i<=20) || (i>90 && i<=100))
      cells(i,30)=2;
    end
    if ((i>20 && i<=30) || (i>80 && i<=90))
      cells(i,30)=3;
    end
    if ((i>30 && i<=40) || (i>70 && i<=80))
      cells(i,30)=4;
    end
    if ((i>40 && i<=50) || (i>60 && i<=70))
      cells(i,30)=5;
    end
    if (i>50 && i<=60)
      cells(i,30)=6;
    end
end

cells;

fileID = fopen('hypoCross_temp.init','w');
temp=[nCell nWall nVer]
fprintf(fileID,'%g   ',temp);
fprintf(fileID,'\n');
for i=1:nWall
    fprintf(fileID,'%g   ',connectivity(i,:));
    fprintf(fileID,'\n');
end

 fprintf(fileID,'\n');

temp=[nVer 3]
fprintf(fileID,'%g   ',temp);
fprintf(fileID,'\n');
for i=1:nVer
    fprintf(fileID,'%g   ',coor(i,:));
    fprintf(fileID,'\n');
end

fprintf(fileID,'\n');

temp=[nWall 1 2]
fprintf(fileID,'%g   ',temp);
fprintf(fileID,'\n');
for i=1:nWall
    fprintf(fileID,'%g   ',walls(i,:));
    fprintf(fileID,'\n');
end
fprintf(fileID,'\n');

temp=[nCell nVar]
fprintf(fileID,'%g   ',temp);
fprintf(fileID,'\n');
for i=1:nCell
    fprintf(fileID,'%g   ',cells(i,:));
    fprintf(fileID,'\n');
end

fclose(fileID);
