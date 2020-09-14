function [points] = read_points(filename)

fileID = fopen(filename,'r');

i = 1;
while(~feof(fileID))
   points(:,i) = str2num(fgetl(fileID))';
   i = i+1 ;
end
points = points(2:end,:);
fclose(fileID);