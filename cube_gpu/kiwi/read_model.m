fileID = fopen('model.dat');
Nbins = fread(fileID,1,'int');
length_voxel = fread(fileID,1,'float');
co = zeros(Nbins,Nbins,Nbins);
for i=1:1:Nbins
    for j=1:1:Nbins
        for k=1:1:Nbins
            co(i,j,k) = fread(fileID,1,'int16');
        end
    end
end
FLAG = fread(fileID,1,'int');
p1 = zeros(FLAG,3);
for i=1:1:FLAG
    p1(i,1) = fread(fileID,1,'int');
    p1(i,2) = fread(fileID,1,'int');
    p1(i,3) = fread(fileID,1,'int');
end
fclose(fileID);
