function faceData = faceData_readLog(textData)
% faceData text file has # of lines = # of frames 
% each line has four parts 
% (1) camera transform 
% (2) face transform
% (3) vertices
% (4) blendshapes 
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
    % vertices data start from 2nd ~ and end at 2+nVertices ~ 
    tline_vertices = tline(b(2)+1:b(nVertices+2)-1);
    c = strrep(strrep(tline_vertices,'~',' '),':',' '); 
    vertices = sscanf(c,'%f',[3,nVertices]); 
    faceData.vertices{i} = vertices'; 

    % blendshape data start from 2+nVertices ~
    tline_bs = tline(b(nVertices+2)+1:end);
    
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