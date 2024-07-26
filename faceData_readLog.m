function faceData = faceData_readLog(textData)
% faceData text file has # of lines = # of frames 
% each line has four parts 
% (1) camera transform 
% (2) face transform
% (3) project transform
% (4) image resolution width x height
% (5) vertices
% (6) blendshapes 
% and these parts are separated by ~
% a vertex is x,y,z coordinate separated by : 
% One vertex and another are also separated by ~
% blendshapes are separated by ~ and blendshapeLocation and its value are
% separated by :
fid = fopen(textData); 

faceData = struct; 
a = regexp(textData,'.txt','split'); 
faceData.ID = a{1}; 

format long 
tline = fgetl(fid); 
nVertices = 1220; % Apple specific 
nBlendShapes = 52; % Apple specific
i = 1; 
while ischar(tline)
    b = regexp(tline,'~'); 
    % camera Transform
    tline_ct = tline(1:b(1)-1);
    c = strrep(strrep(tline_ct,'~',' '),':',' ');
    ct = sscanf(c,'%f',[4,4]); 
    faceData.ct{i} = ct;

    % face Transform
    tline_ft = tline(b(1)+1:b(2)-1);
    c = strrep(strrep(tline_ft,'~',' '),':',' ');
    ft = sscanf(c,'%f',[4,4]); 
    faceData.ft{i} = ft;
   
    % projection matrix
   tline_p = tline(b(2)+1:b(3)-1);
    c = strrep(strrep(tline_p,'~',' '),':',' ');
    ft = sscanf(c,'%f',[4,4]); 
    faceData.pt{i} = ft;

    %imageResolution
   tline_w = tline(b(3)+1:b(4)-1);
   faceData.imageResolution{i} = zeros(1,2);
    width = sscanf(tline_w,'%f'); 

    tline_h = tline(b(4)+1:b(5)-1);
    height= sscanf(tline_h, '%f');

    faceData.imageResolution{i}(1) = width;
    faceData.imageResolution{i}(2) = height;


    % vertices data start from 5th ~ and end at 5+nVertices ~ 
    tline_vertices = tline(b(5)+1:b(nVertices+5)-1);
    c = strrep(strrep(tline_vertices,'~',' '),':',' '); 
    vertices = sscanf(c,'%f',[3,nVertices]); 
    faceData.vertices{i} = vertices'; 

    % blendshape data start from 5+nVertices ~
    tline_bs = tline(b(nVertices+5)+1:end);
    
    d = regexp(tline_bs,'[\w.-]+','match');
    bsLoc = d(1:2:end)';
    faceData.blendshape(:,i) = cellfun(@str2num,d(2:2:end))';
    
    tline = fgetl(fid); 
    i = i+1;
end

% sort blendshape data 
% blendshape data can be different session by session 
[bsLoc_sorted,I] = sort(bsLoc);
faceData.blendshape = faceData.blendshape(I,:);
faceData.bsLoc = bsLoc_sorted; 
end