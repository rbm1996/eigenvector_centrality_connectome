%% example script : Eigenvector centrality
subdir = genpath("utils");
addpath(subdir);


%% load connectome
tableName = ...
'utils/0000000_connectome_edgelist_schaefer200-yeo17_rmap_s1dild.csv';

matrixSize = 216;
edgeList = table2array(readtable(tableName));
i = edgeList(:,2) +1;
j = edgeList(:,3) +1;
w = edgeList(:,4);
% build matrix
matrix = full(sparse(i,j,w , matrixSize  , matrixSize));

% add upper part
matrix = matrix+ matrix';

%%% get only real parcels
% additional parcel created by multiAtlasTransferTool
% ignore FreeSurfer medial Wall 
% lh.Background+FreeSurfer_Defined_Medial_Wall.label

realIdx = [(2 : 100) (102 :216)]; % this valid for Scaefer200 % change as needed

matrix = matrix(realIdx , realIdx);
n = size(matrix,2);
%%%%%%% the resulting network is undirected with zero diagonal

%%%% example replace connectome with random graph 
% matrix = rand(n,n); % very dense network!!


%%% !!!!!!! One can load a network with all different ways possible as long
%%% as the variable "matrix" is a matlab array of size nxn
leftIdx = [1:100 201:207];
rightIdx = [101:200 208:214];

nl = length(leftIdx);
nr = length(rightIdx);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
A = matrix;
A=A-diag(diag(A));
% get the smallest link
lL = min(A(A ~= 0));

Aviz = log10(1 + A / lL);
figure
imagesc(Aviz); %only for viz but not for computation
axis square
colorbar

%% Check the EC
%%% get other metrics
strength = sum(A);
strLeft = strength(leftIdx)';
strRight = strength(rightIdx)';

[ec , maxEig] =eigenvector_centrality_und_fab(A);

normEc = ec/sum(ec);

normStr = strength/sum(strength);


ecLeft = normEc(leftIdx);
ecRight = normEc(rightIdx);

red = [0.8500 0.3250 0.0980];
blue = [0 0.4470 0.7410];

figure;
bb = bar([sum(ecLeft) sum(ecRight)]);

bb.FaceColor = 'flat';
bb.CData(1,:) = red;
bb.CData(2,:) = blue;

yline(0.5);
ax = gca;
ax.XTickLabel = {'Left' , 'Right'};
ax.FontSize = 14;
ylabel('EC')
%% Vizualization

%%%% load nodes coordinates 
ROI_info_Table = readtable( 'utils/Schaefer200_allinfo.csv');

coord = table2array(ROI_info_Table(:,7:9));
coord = coord - ones(n,1)* mean(coord);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[CIJtree,CIJclus] = backbone_wu(matrix,6);


figure;
G = graph(CIJtree);

markerSize = 0.1 + 3* log10(1 + ec / min(ec(ec>0)));
markerSize = 1.5 + 10*(ec - min(ec(ec>0)))/ (max(ec - min(ec(ec>0))));

red = [0.8500 0.3250 0.0980];
blue = [0 0.4470 0.7410];

color = zeros(n,3);
color(leftIdx,:) = ones(nl,1) *red;
color(rightIdx,:) = ones(nr,1) *blue;


lab = [];
for k = 1:n
    lab = [lab ; ''];
end

plot(G , 'MarkerSize' , markerSize, 'NodeColor' , color ,  'XData' , coord(:,1) , 'NodeLabel' , lab ,...
    'YData' , coord(:,2) , 'ZData' , coord(:,3) , 'EdgeAlpha' , 0.1 , 'EdgeColor' , 'k' , 'LineWidth',2)

axis equal tight off
title('Eigenvector centrality')
