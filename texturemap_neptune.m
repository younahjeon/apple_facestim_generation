% texturemap neptune
% texture map elias using project Transform and image resolution obtained
% from a new xcode project

fname = '/Users/younah/Desktop/Lab/FaceData/1221_neptune2/neptune.txt';
faceData = faceData_readLog_old(fname); 


nVertices = 1220;

frame_to_use = 1; 

v = faceData.vertices{frame_to_use};

nv = [v, ones(nVertices,1)];
ct =faceData.ct{frame_to_use};

ft = faceData.ft{frame_to_use};
cV = ft * nv';

% from world to camera
cV = pinv(ct) * cV;
cV = cV(1:3,:)';

% camera to image 
load('projT_iphoneX.mat')
nv = [cV, ones(nVertices,1)];
projectedCoord=  projT_iphoneX * nv'; 
projectedCoord= projectedCoord(1:3,:)./projectedCoord(4,:); % normalize
projectedCoord = projectedCoord'; 

% change to uv coordinate
x = projectedCoord(:,1) *0.5 + 0.5;
y = projectedCoord(:,2) * 0.5 + 0.5; 
y = 1 - y; % image coordinates are y-flipped 

% to pixels 
pixelCoord = zeros(nVertices,1);
pixelCoord(:,1) = x*1280;
pixelCoord(:,2) = y *720;

Options = struct;
I= imread('0000_max.png');
%I = imread('/Users/younah/Desktop/Lab/FaceData/20240807_max/0000_edited.png');
%I = imread('/Users/younah/Desktop/Lab/FaceData/ryan.png');
%I = imread('/Users/younah/Desktop/Lab/FaceData/lucas.png');
load('TriangleIndices.mat');
% 
% figure; 
% trimesh(TriangleIndices', v(:,1),v(:,2),v(:,3),'EdgeColor',[0,0,0],'FaceColor',[1,1,1],'FaceAlpha',0);
% figure; 
% 
% texturemap = patcht(TriangleIndices', v, TriangleIndices',flip(pixelCoord,2),I, Options);
% 
% % get vertex color
% 
% % average of the color of the faces that contain that vertex
% vertexColor = zeros(nVertices,3);
% for i = 1:nVertices
%     [allrows, idx] = findRows(TriangleIndices',i);
%     meancolor = mean(texturemap(:,:,:,find(idx==1)),4);
%     % average of the tesselations 
%     mean_meancolor = squeeze(mean(meancolor,[1,2]));
%     vertexColor(i,:) = mean_meancolor';
% end
% 
% vertexColor = vertexColor/255; 
% newmesh = surfaceMesh(v, TriangleIndices','VertexColors', vertexColor);
% 
% surfaceMeshShow(newmesh)
% %%


%% new UV coordinates



load('FaceParts_neptune.mat')


ind_part = unique(FaceParts.ALLs{1})';

UV = ones(size(newMesh{2},1),2);
UV(newMesh{3},:) = [x(ind_part),y(ind_part)];
%%

nVertices = size(newMesh{2},1);

nv_mat = [newMesh{2}(newMesh{3},:), ones(length(newMesh{3}),1)];
UV_mat = [UV(newMesh{3},:),ones(length(newMesh{3}),1)];
mat_proj = UV_mat'/nv_mat';

new_UV_mat = [newMesh{2}, ones(length(newMesh{2}),1)] *mat_proj';

figure; 
imshow(I)
hold on
scatter(new_UV_mat(:,1) * 1280 , new_UV_mat(:,2) *720);


%%


% flatten 3d
v_flattened = newMesh{2};
figure; 

trimesh(newMesh{1}, newMesh{2}(:,1),newMesh{2}(:,2),newMesh{2}(:,3),'EdgeColor',[0,0,0],'FaceColor',[1,1,1],'FaceAlpha',0);
xlabel('x')
ylabel('y')
zlabel('z')

% left
vind = find(newMesh{2}(:,1) <0& newMesh{2}(:,3)<-0.01);%& newMesh{2}(:,2)>-0.09);
hold on
scatter3(newMesh{2}(vind,1),newMesh{2}(vind,2),newMesh{2}(vind,3), 'blue','filled')


p_orig= newMesh{2}(vind,:);
rotAxis = [0,1,0];
rotAngle = -80;

R = vrrotvec2mat([rotAxis,rotAngle]); % deprecated in 2024 version
R = [[R;0,0,0],[0;0;0;1]];

center_delta = mean(p_orig);

pos_z_ind  = find(p_orig(:,3)>-0.02);
[val,minind] = min(p_orig(pos_z_ind,1));
maxind =  pos_z_ind(minind);
%[val,maxind] = max(p_orig(:,3));

scatter3(p_orig(maxind,1),p_orig(maxind,2),p_orig(maxind,3),200,'magenta','filled')

p = p_orig-center_delta;
p = [p,ones(size(p,1),1)];
pfinal = R*p';

p_new = pfinal(1:3,:)';
scatter3(p_new(:,1),p_new(:,2),p_new(:,3),'cyan')

scatter3(p_new(maxind,1),p_new(maxind,2),p_new(maxind,3),200,'magenta','filled')

p_delta= p_new(maxind,:) - p_orig(maxind,1:3);
p_new = p_new - p_delta;
scatter3(p_new(:,1),p_new(:,2),p_new(:,3),'cyan')

v_flattened(vind,:) = p_new;

%right
vind = find(newMesh{2}(:,1) >0 & newMesh{2}(:,3)<-0.01);

hold on
scatter3(newMesh{2}(vind,1),newMesh{2}(vind,2),newMesh{2}(vind,3), 'red','filled')

p_orig= newMesh{2}(vind,:);
rotAxis = [0,1,0];
rotAngle = 80;

R = vrrotvec2mat([rotAxis,rotAngle]); % deprecated in 2024 version
R = [[R;0,0,0],[0;0;0;1]];

center_delta = mean(p_orig);
[val,maxind] = max(p_orig(:,3));

pos_z_ind  = find(p_orig(:,3)>-0.02);
[val,minind] = max(p_orig(pos_z_ind,1));
maxind =  pos_z_ind(minind);


scatter3(p_orig(maxind,1),p_orig(maxind,2),p_orig(maxind,3),200,'magenta','filled')

p = p_orig-center_delta;
p = [p,ones(size(p,1),1)];
pfinal = R*p';

p_new = pfinal(1:3,:)';
scatter3(p_new(:,1),p_new(:,2),p_new(:,3),'r')

scatter3(p_new(maxind,1),p_new(maxind,2),p_new(maxind,3),200,'magenta','filled')

p_delta= p_new(maxind,:) - p_orig(maxind,1:3);
p_new = p_new - p_delta;
scatter3(p_new(:,1),p_new(:,2),p_new(:,3),'r')

v_flattened(vind,:) = p_new;
%%
figure; 

nv_mat = [newMesh{2}(newMesh{3},:), ones(length(newMesh{3}),1)];
UV_mat = [UV(newMesh{3},:),ones(length(newMesh{3}),1)];
mat_proj = UV_mat'/nv_mat';

new_UV_mat2 = [v_flattened, ones(length(newMesh{2}),1)] *mat_proj';
new_UV_final = new_UV_mat2;
I= imread('0000_max_edit.png');

new_UV_final(:,2) = (new_UV_mat2(:,2) * 720+160) / size(I,1);

new_UV_final(:,1) = (new_UV_mat2(:,1)*1280-5 )/size(I,2);

% I= imread('0000_max.png');
% new_UV_final(:,2) = (new_UV_mat2(:,2) * 720) / size(I,1);
% 
% new_UV_final(:,1) = (new_UV_mat2(:,1)*1280 )/size(I,2);


figure;
imshow(I)
hold on
%scatter(new_UV_mat2(:,1) * size(I,2), new_UV_mat2(:,2) * size(I,1))

scatter(new_UV_final(:,1) * size(I,2), new_UV_final(:,2) * size(I,1))
% color the boundary
%k = boundary(new_UV_final(:,1), new_UV_final(:,2));
%plot(new_UV_final(k,1) * size(I,2),new_UV_final(k,2) * size(I,1))


%%
neptune = GLTF();
% Add material with image texture
material_idx=neptune.addMaterial('baseColorTexture',"0000_max_edit.png");
% Add mesh with UV coordinates

mesh_idx=neptune.addMesh(newMesh{2},'indices',newMesh{1},'material',material_idx,'TEXCOORD',new_UV_final(:,1:2));


% Add node
neptune.addNode('mesh',mesh_idx);

% Write GLTF file
neptune.writeGLTF("neptune.gltf");


%%

%% all vertices on the right side

vind = find(v(:,1) <-0.01);
length(vind)
% get vertices that are on the opposite side
vind_opp = [];
for i = 1:length(vind)
    new_v = v(vind(i),:);
    new_v(1) = -new_v(1);
    [val,ind] = min(sum(abs(v - new_v),2));
    vind_opp = [vind_opp,ind];
end

f = TriangleIndices';
h = figure;
%trimesh(TriangleIndices',v(:,1),v(:,2),v(:,3),'FaceVertexCData',[1,1,1],'EdgeColor',[0,0,0],'FaceColor','flat');
trimesh(f,v(:,1),v(:,2),v(:,3),'FaceVertexCData',[1,1,1],'EdgeColor',[0,0,0],'FaceColor','flat');
xlabel('x')
ylabel('y')
zlabel('z')
hold on
scatter3(v(vind,1),v(vind,2),v(vind,3), 'red','filled')
scatter3(v(vind_opp,1),v(vind_opp,2),v(vind_opp,3), 'blue','filled')

% get vertexcolor

new_vertexColor = vertexColor;

new_vertexColor(vind_opp,:) = vertexColor(vind,:);

%newmesh = surfaceMesh(v, TriangleIndices','VertexColors', new_vertexColor);
newmesh = surfaceMesh(v, f,'VertexColors', new_vertexColor);
surfaceMeshShow(newmesh)

%% 
new_vertexColor2 = new_vertexColor; 


load('FaceParts_neptune.mat')


ind_part = unique(FaceParts.ALLs{1})';
part_color = new_vertexColor2(ind_part,:);

%nVertices = size(eliasv,1);

nVertices = size(newMesh{2},1);
all_vc = zeros(nVertices,3); 


%all_vc(eliasvind,:) = part_color;
all_vc(newMesh{3},:) = part_color;

newmesh = surfaceMesh(newMesh{2}, newMesh{1}, 'VertexColors', all_vc); 

surfaceMeshShow(newmesh)
%%
%% pick target color
h = figure;
%trimesh(TriangleIndices',v(:,1),v(:,2),v(:,3),'FaceVertexCData',[1,1,1],'EdgeColor',[0,0,0],'FaceColor','flat');
trimesh(f,v(:,1),v(:,2),v(:,3),'FaceVertexCData',[1,1,1],'EdgeColor',[0,0,0],'FaceColor','flat');

datacursormode(h)
dcm_obj = datacursormode(h);
set(dcm_obj,'UpdateFcn',@myupdatefcn)
%%
info_struct = getCursorInfo(dcm_obj);

[val,ind] = min(sum(abs(v - info_struct.Position),2));
targetInd = ind;
%%

nVertices = size(newMesh{2},1);
all_vc = nan(nVertices,3); 
all_vc(newMesh{3},:) = part_color;
targetColor = mean(part_color);
targetColor = new_vertexColor2(322,:);
targetColor = new_vertexColor2(targetInd,:);

    
% targetColor is the skin color
% take the average of colors in nose area

x = [0 1 1 0] ; y = [0 0 1 1] ;
figure
fill(x,y,targetColor)


%%

tbInd = newMesh{7};
new_coloredV = newMesh{3};


n_step = 1;
tb_per_step = {}; 
len_vind = [];
while length(new_coloredV) <nVertices
%while n_step<2
    
    if n_step > 1
        tbInd_new = [];
        for i = 1:length(tbInd)
            [rows, idx] = findRows(newMesh{1},tbInd(i));
            vind = unique(rows);
            tbInd_new = [tbInd_new, vind(~ismember(vind,new_coloredV))']; 
        end
    else
        tbInd_new = tbInd;
    end
    
    tb_to_include= [];
    for i = 1:length(tbInd_new)
        [rows,idx] = findRows(newMesh{1}, tbInd_new(i));
    
        vind = unique(rows);
        vind = vind(ismember(vind,new_coloredV));
        len_vind = [len_vind,length(vind)];
        if ~isempty(vind)
            mean_color = nanmean(cat(1,all_vc(vind,:),targetColor));
            %mean_color = nanmean(all_vc(vind,:),1);
            all_vc(tbInd_new(i),:) = mean_color;
            tb_to_include = [tb_to_include, tbInd_new(i)];
        end

    end
  
    tbInd_new = unique(tb_to_include);
    
    new_coloredV  = unique([new_coloredV, tbInd_new]);
    length(new_coloredV)
    
    tbInd = tbInd_new; 
    n_step = n_step +1;

    tb_per_step = [tb_per_step, tbInd];    
end


%%

all_vc_nonans = all_vc;
nanind = find(isnan(all_vc_nonans(:,1)));
all_vc_nonans(nanind,:) = zeros(length(nanind),3);

newmesh = surfaceMesh(newMesh{2}, newMesh{1}, 'VertexColors', all_vc_nonans); 

surfaceMeshShow(newmesh)


vertexColor = all_vc_nonans;
save('neptune_pt_vc_lucas','vertexColor')
%%
writeSurfaceMesh(newmesh, 'neptune_lucas.ply')






%%

load('meshinfo_neptune.mat')
load('neptune_pt_vc_max.mat')
all_vc_nonans = vertexColor;
newmesh = surfaceMesh(newMesh{2}, newMesh{1}, 'VertexColors', all_vc_nonans); 

surfaceMeshShow(newmesh)


