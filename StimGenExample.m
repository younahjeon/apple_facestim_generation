% Example face generation
% All figures are saved in the figures folder.

% LOAD PART MESH AND BLANK HEAD
% Load FaceParts_elias, which is a structure that contains information
% about different face parts of sophie 
% (1) triangular face indices, (2) coordinates,
% (3) indices of the outline of the part, 
% (4) and reference index/indices

global targetCoord baseF baseV partV partF partB refInd baseTR baseVN baseVA baseFN partIndices 
baseMesh = 'PotatoHead';
id_string = 'example'; 
part_string = 'ALLs'; % all parts but smaller 
emotion ='';
scale = 0.9;
%emotion = 'smile';

[baseF, baseV, partF, partV_orig, partB,refInd,baseTR,baseVN,baseVA,baseFN,...
    partIndices, partTR, partVN, partVA, partFN, partRotAxes] = loadMaterials(id_string,part_string,baseMesh,scale,emotion);

% loadMaterials produces two figures, one of the blank head and the second
% of the part. The reference indices of the part are colored red in the
% second figure.
%%
% 1. rotate the part so that the reference point has the same vertex normal
% as the target vertex normal 
refVN = partVN(refInd(2),:);
targetVN = [0,0,1]; 

partV = rotateMesh(partV_orig,'vertexnormal',refVN,targetVN); 

%%
% 2. Define the location on the blank head 
% where the reference point will be attached to. You can additionally
% define the reference point of the part. For example, ALLs contains two
% reference points

global refInd_optim
refInd_optim = refInd(1);

[targetCoord targetVN]= targetSelection(baseF,baseV);

% or you can provide targetCoord as below
targetCoord = [-0.002693421673030   0.049340374767780   0.054754320532084];

%%
% 2. do mesh surgery 
smoothingArea = 100; 
% smoothingArea = 100; % if you want to stop at this step
t = 1; % step of FaceMorph
newMesh = {}
[newMesh,optim_param]= MeshSurgery(t,part_string,smoothingArea,newMesh);

TR = triangulation(newMesh{1,t},newMesh{2,t}/max(newMesh{2,t},[],'all'));
%filename = strcat(id_string,'.stl');
%stlwrite(TR,filename)


%%
% animate the resulting face with arbitrary blendshapes

t_to_animate = 1; % 3 seconds
fps = 10; 
n_frames = t_to_animate * fps; 
writeMovie = 1;

v = VideoWriter('movie.mp4','MPEG-4');
open(v)
load('BETA.mat')

nVertices = 1220;
%X = [0.104969700000000;0.104980100000000;0.114777400000000;0.00303261400000000;0.00296590700000000;0.00666228400000000;0.237118100000000;0.238192200000000;-2.18774000000000e-12;-2.18845600000000e-12;0.0560495000000000;0.0560571100000000;0.102377000000000;0.00716774000000000;0;0;0;0;0.685995100000000;0.685998300000000;0.0816218400000000;0.0811506100000000;0.143178100000000;1.53680700000000e-17;0.151617800000000;0.0437056300000000;0.0116943600000000;0.110847700000000;0.110109400000000;0;0;0.0930285600000000;0.0510150900000000;0.792670200000000;0.786838200000000;0.0177960500000000;0.0201759400000000;0.0215139600000000;1.64507500000000e-05;0.0232451400000000;0.0340266800000000;0.00674136400000000;0.287926900000000;0.917157500000000;0.917360100000000;0.212382300000000;0.198032400000000;0.244017900000000;0.245387400000000;0.279677200000000;0.270702900000000;0.000261358700000000]
X = zeros(52,1);
X(9) = 1;
X(10) = 1;
X(length(X)+1) = 0;
delta = X'* BETA; 
delta = reshape(delta,nVertices,3); 
load('FacePartsInd.mat'); 
ALLs_ind = FacePartsInd.ALLs;
delta_sub = delta(ALLs_ind,:);

neutral = newMesh{2};

vertices_orig = neutral;

vertices_target = vertices_orig;
vertices_target(newMesh{3},:) = vertices_orig(newMesh{3},:) + delta_sub;

vertices_inc = (vertices_target - vertices_orig)/ n_frames; 
figure;
%set(gcf, 'renderer', 'zbuffer');
for i = 1:n_frames
    vertices = vertices_orig + (i-1) * vertices_inc;
   
    trimesh(newMesh{1},vertices(:,1),vertices(:,2),vertices(:,3),...
         'FaceColor',[0.65,0.65,0.65],'EdgeColor','none');
    
    cameraViewAngle = 10; % field of view (0 to 180) The greater the number, smaller objects will show more
    cameraTarget = [32 10.785 14.3812]; % how far your target is from the camera. determines the degree of visibility
    cameraPos = [0 0 0];
    set(gca,'CameraViewAngle',cameraViewAngle,'CameraTarget',cameraTarget,'CameraPosition',cameraPos);
    axis equal
    view (0,90)
    
    h1 = camlight(45,25); lighting gouraud;
    set(h1,'color',[.7 .7 .7]);
    h3 = camlight('headlight'); lighting gouraud
    set(h3,'color',[0.2 0.2 0.2]);
    h2=camlight(-45,25); lighting gouraud
    set(h2,'color',[0.7 0.7 0.7]);
    material([.3 .8 .1 10 1]);
    set(gcf,'color','w');
    set(gcf,'menubar','none');
    axis off
    if i ==1
    set(gcf,'position',[1 1 1.15 1.15].*get(gcf,'position'));
    view (0,90)
    end
    title(sprintf('animation time = %2.2f seconds \nframe = %d',t_to_animate,i));
    drawnow 
    pause(0.01)
    if writeMovie
        frame = getframe(gcf);
        writeVideo(v,frame);
    end
end
close(v); % crucial! 