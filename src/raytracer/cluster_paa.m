function [PAA_info]  = cluster_paa(nodeLoc, nodePAA_position, nodePAA_Orientation, varargin)
%CLUSTER_PAA returns the struct PAA_info given in input the position of the
%nodes nodeLoc and the position of the PAAs each node nodePAA_position.
%
%   [PAA_info]  =  CLUSTER_PAA(nodeLoc, nodePAA_position)
%   **nodeLoc is the TxNx3 matrix containing the 3D coordinates of the N
%   nodes at different T time instances
%   **nodePAA_position is the Nx1 cell array containing the 3D coordinates of
%   each PAA in the N nodes at different T time instances. If PAAs are not 
%   defined input is the Nx1 cell array of empty vectors {[], [], ..}
%
%   PAA_info{i}:
%     **nPAA_node: number of PAA in node i
%     **node_clusters: cell array of node clustered
%     **centroids: centroid of each cluster
%     **generationMethod: channel generation method. 0: common channel 1:
%     common deterministic part 2: indipendent channels
%     **PAA_loc: matrix of centroid positions
%     **nodePAAInfo: flags for interfacing with Raytracer
%     **PAA_position: matrix of unique centroid positions. If no PAAs are
%     defined in input, the output is the location of the node
%     **nPAA_centroids: number of unique centroids
%
%   [PAA_info, nodeLoc2,nodePAAInfo]  =  CLUSTER_PAA(nodeLoc, nodePAA_position, option, value)
%
%   option 'freq': specify central frequncy in Hz
%          'l_cor': specify correlation distance in wavelengths
%
%   Copyright 2019-2020 NIST/CLT (steve.blandino@nist.gov)

%#codegen

%% Input processing
var_struct = {'fc','l_cor'};
for k = 1:2:length(varargin)
    if (~any(strcmp(varargin{k}, var_struct)))
        warning(['Cannot specify "', varargin{k}, '" as input value - it will be discarted']);
    end
    eval([varargin{k},' = varargin{k+1};'])
end

if ~exist('fc','var'),     fc=60e9;        end
if ~exist('l_cor','var'),     l_cor = 50;   end

assert(size(nodeLoc,ndims(nodeLoc)-1) == 3, 'Provide (x,y,z) coordinates as input of cluster_paa')
assert(size(nodeLoc,ndims(nodeLoc)) == length(nodePAA_position), 'Provide correct input: each node should have a correspondent PAAs description')

C = 3e8;
wavelength =C/fc;
squeezeAndReshape = @(x) reshape(squeeze(x), [], 3); % Helper anonymus for using unique when multiple time divisions
findRow = @(x,y) sum(ismember(x,y),2)>0;

if ismatrix(nodeLoc)
    nodeLocTmp = zeros([1,flip(size(nodeLoc))]);
    nodeLocTmp(1,:,:) = nodeLoc.';
    nodeLoc = nodeLocTmp;
else
    nodeLoc = permute(nodeLoc, [1 3 2]);
end

%% Process PAA positions
numberOfNodes = length(nodePAA_position);
PAA_distance_paaCombinations = cell([1,numberOfNodes]);
PAA_info = cell([1,numberOfNodes]);
numberTimeDivision = size(nodeLoc,1);

%% Clustering: loop on nodes
for node_id = 1:numberOfNodes
    nPAA = size(nodePAA_position{node_id},1); % Number of PAAs at node node_id
    if nPAA == 1, nPAA =0;end
    PAA_info{node_id}.nPAA_node = nPAA;
%     PAA_info{node_id}.orientation = nodePAA_Orientation{node_id};
    PAA_info{node_id}.centroids = 1;
    PAA_node_to_cluster = 1:nPAA; % Indexes of nodes to cluster
    node_clustered = [];          % Nodes in the cluster
    PAALoc = [];
    idx_paa=0;
    paaIdComb  = nchoosek(1:nPAA ,2); % Indexes of all PAA combinations
    
    PAA_distance_paaCombinations{node_id} = sqrt(...
        sum(abs(...
        nodePAA_position{node_id}(paaIdComb(:,1),:)-...
        nodePAA_position{node_id}(paaIdComb(:,2),:)...
        ).^2, ...
        2)...
        ); % Compute pairwise distances
    
    % Find possible clusters of PAA closer than l/2
    if any(PAA_distance_paaCombinations{node_id}<= wavelength/2+eps)
        idx_paa=idx_paa+1;
        idx = paaIdComb(PAA_distance_paaCombinations{node_id}<= wavelength/2+eps,:);% This couples can be combined
        node_to_cluster = unique(idx(:)); % Remove redundancy
        while numel(node_to_cluster)
            [~, popular_node_id] = max(histc(idx(:), node_to_cluster)); % Find node with most connections
            PAA_info{node_id}.node_clusters{idx_paa} = reshape(...
                unique(idx(findRow(idx,node_to_cluster(popular_node_id)), :)),...
                [], 1); % Cluster of nodes connected with popular node
            PAA_info{node_id}.centroids(idx_paa) = node_to_cluster(popular_node_id);
            PAALoc(:,idx_paa, :) =squeezeAndReshape(nodeLoc(:,node_id, :)) + ...
                repmat(nodePAA_position{node_id}(PAA_info{node_id}.centroids(idx_paa),:), ...
                [numberTimeDivision, 1,1]);%#ok<*EMVDF> % Add PAA shifts to PAA locations
            PAA_info{node_id}.generationMethod(idx_paa) = 0;
            node_clustered = [node_clustered; PAA_info{node_id}.node_clusters{idx_paa}]; %#ok<AGROW>
            PAA_info{node_id}.centroids_shift{idx_paa} = repmat(nodePAA_position{node_id}(PAA_info{node_id}.centroids(idx_paa),:), ...
                length(PAA_info{node_id}.node_clusters{idx_paa}),1)-nodePAA_position{node_id}(PAA_info{node_id}.node_clusters{idx_paa} ,:); % Save shift from centroid
            PAA_info{node_id}.orientation{idx_paa} = nodePAA_Orientation{node_id}(PAA_info{node_id}.node_clusters{idx_paa}, :);
            node_to_cluster = node_to_cluster(~ismember(node_to_cluster,  PAA_info{node_id}.node_clusters{idx_paa})); % Node not considered in cluster
            idx_paa = idx_paa+1;
        end
        idx_paa = idx_paa-1;
        % Check if cluster can be further merged
        if any(histc(node_clustered, unique(node_clustered))>1)
            unique_node_clustered = unique(node_clustered);
            [rep, ~] = histc(node_clustered, unique_node_clustered);
            merge_paa = cell2mat(cellfun(@(x) any(ismember(x,unique_node_clustered(rep>1))), PAA_info{node_id}.node_clusters, 'UniformOutput', 0));
            node2merge = cellfun(@(x) x(~ismember(x,unique_node_clustered(rep>1))) , PAA_info{node_id}.node_clusters(merge_paa), 'UniformOutput', 0);
            node_merged = sort([ vertcat(node2merge{:}); unique_node_clustered(rep>1)]);
            cluster_idx = find(merge_paa, 1 );
            merge_paa = find(merge_paa);
            cluster2del =  merge_paa((merge_paa~=cluster_idx));
            PAA_info{node_id}.node_clusters{cluster_idx} = node_merged;
            getCentralValue = @(x) x(ceil(end/2));
            PAA_info{node_id}.centroids(cluster_idx) = getCentralValue(unique_node_clustered(rep>1));
            PAA_info{node_id}.generationMethod(cluster_idx) = 0;
            %                 PAA_info{node_id}.node_clusters{cluster2del} = cellfun(@(x) [], PAA_info{node_id}.node_clusters(cluster2del), 'UniformOutput', 0);
            PAA_info{node_id}.node_clusters = PAA_info{node_id}.node_clusters(~ismember(1:numel(PAA_info{node_id}.node_clusters), cluster2del)); % Double-check this line
            PAA_info{node_id}.centroids(cluster2del) = [];
            PAA_info{node_id}.generationMethod(cluster2del) = [];
            PAALoc(:,cluster2del, :) = [];
            PAA_info{node_id}.centroids_shift{cluster_idx} = repmat(nodePAA_position{node_id}(PAA_info{node_id}.centroids(cluster_idx),:), length(node_merged),1)-nodePAA_position{node_id}(PAA_info{node_id}.node_clusters{cluster_idx} ,:);
            PAA_info{node_id}.centroids_shift(cluster2del) = [];
            idx_paa = idx_paa-1;
        end
    end
    
    % Find possible clusters of PAA between l/2 and l_cor*lambda
    indepPAAidx = PAA_node_to_cluster(~ismember(PAA_node_to_cluster,node_clustered) ); % Node not clustered yet
    row_idx = find(findRow(paaIdComb, indepPAAidx));
    
    if  any(PAA_distance_paaCombinations{node_id}(row_idx,:)<= l_cor*wavelength/2 ...
            & PAA_distance_paaCombinations{node_id}(row_idx,:)> wavelength/2)
        idx = paaIdComb(row_idx(PAA_distance_paaCombinations{node_id}(row_idx,:)<= l_cor*wavelength/2 & PAA_distance_paaCombinations{node_id}(row_idx,:)> wavelength/2),:);% This couples can be combined
        if any(~ismember(idx(:),node_clustered))
            row_with_centroid = idx(findRow(idx, PAA_info{node_id}.centroids),:);  % This couples can be combined
            for ct = 1:numel(PAA_info{node_id}.centroids)
                idx_paa = idx_paa+1;
                cluster_tmp = unique(row_with_centroid(findRow(row_with_centroid, PAA_info{node_id}.centroids(ct)),:));
%                 PAA_info{node_id}.node_clusters{idx_paa}  = cluster_tmp(cluster_tmp ~= PAA_info{node_id}.centroids(ct));
                PAA_info{node_id}.node_clusters{idx_paa}  = cluster_tmp(~ismember(cluster_tmp,node_clustered));
                PAA_info{node_id}.centroids(idx_paa) =PAA_info{node_id}.centroids(ct);
                PAA_info{node_id}.centroids_shift{idx_paa}  = ...
                repmat(nodePAA_position{node_id}(PAA_info{node_id}.centroids(idx_paa),:),...
                length(PAA_info{node_id}.node_clusters{idx_paa}),1)-nodePAA_position{node_id}(PAA_info{node_id}.node_clusters{idx_paa} ,:);
                PAA_info{node_id}.orientation{idx_paa} =  nodePAA_Orientation{node_id}(PAA_info{node_id}.node_clusters{idx_paa}, :);
                PAALoc(:,idx_paa, :) =squeezeAndReshape(nodeLoc(:,node_id, :)) + nodePAA_position{node_id}(PAA_info{node_id}.centroids(idx_paa),:);
                PAA_info{node_id}.generationMethod(idx_paa) = 1;
                node_clustered = [node_clustered;PAA_info{node_id}.node_clusters{idx_paa}]; %#ok<AGROW>
            end
        end
    end
    
    
    % No cluster found
    if any(~ismember(1:nPAA,node_clustered))
        vec_paa = 1:nPAA;
        missing_paa = vec_paa(~ismember(1:nPAA,node_clustered));
        
        for m_id = 1:numel(missing_paa)
            PAA_info{node_id}.node_clusters{idx_paa+m_id} = missing_paa(m_id);
            PAA_info{node_id}.centroids(idx_paa+m_id) = missing_paa(m_id);
            PAA_info{node_id}.centroids_shift{idx_paa+m_id}  = zeros(1,3);
            PAA_info{node_id}.orientation{idx_paa+m_id} = nodePAA_Orientation{node_id}(missing_paa(m_id), :);
            PAALoc(:,idx_paa+m_id, :) =squeezeAndReshape(nodeLoc(:,node_id, :)) + nodePAA_position{node_id}(PAA_info{node_id}.centroids(idx_paa+m_id),:);
            PAA_info{node_id}.generationMethod(idx_paa+m_id)=2;
        end
    end
    
    %% Finalize struct
    if PAA_info{node_id}.nPAA_node ~=0
        PAA_info{node_id}.PAA_loc = PAALoc;
    else
        PAA_info{node_id}.nPAA_node = 1;
        PAA_info{node_id}.node_clusters = {1};
        PAA_info{node_id}.centroids = 1;
        PAA_info{node_id}.centroids_shift =nodePAA_position(node_id);
        PAA_info{node_id}.PAA_loc =nodeLoc(:,node_id, :);
        PAA_info{node_id}.orientation{1} = nodePAA_Orientation{node_id};
    end    
end

%% Interface with channel model
paa_id_init = 0;
% tot_number_paa = sum(arrayfun(@(x) length(unique(x{:}.centroids)), PAA_info));
% nodeLoc2 = zeros(tot_number_paa,3);
% st_id = 1;
for node_id  = 1:numberOfNodes
    [~, id ] = unique(squeezeAndReshape(PAA_info{node_id}.PAA_loc(1,:,:)), 'rows', 'stable');
    unique_PAA_location = PAA_info{node_id}.PAA_loc(:,id,:);
%     en_id = st_id + size(unique_PAA_location,1)-1;
%     nodeLoc2(st_id:en_id,:) = unique_PAA_location;
%     st_id = en_id + 1;
    [unique_centroids,position_centroids] = unique(PAA_info{node_id}.centroids);
    [~,position_centroids_sorted] = sort(position_centroids);
    centroid_rep = histc(PAA_info{node_id}.centroids, unique_centroids);
    centroid_rep = centroid_rep(position_centroids_sorted);
    size_paa= cellfun(@(x) length(x), PAA_info{node_id}.node_clusters);
%     nodePAAInfo = cell(numel(unique_centroids),1);
    unique_centroids=unique_centroids(position_centroids_sorted);
    for paa_id = 1:length(centroid_rep)
        PAA_info{node_id}.nodePAAInfo{paa_id_init+paa_id,1}.centroid_id = unique_centroids(paa_id);
        PAA_info{node_id}.nodePAAInfo{paa_id_init+paa_id,1}.tot_channel = sum(size_paa(PAA_info{node_id}.centroids==unique_centroids(paa_id)));
        PAA_info{node_id}.nodePAAInfo{paa_id_init+paa_id,1}.indep_stoch_channel = centroid_rep(paa_id);
        PAA_info{node_id}.nodePAAInfo{paa_id_init+paa_id,1}.rotated_channel = size_paa(PAA_info{node_id}.centroids==unique_centroids(paa_id));
        PAA_info{node_id}.nodePAAInfo{paa_id_init+paa_id,1}.paa_id = cell2mat(PAA_info{node_id}.node_clusters(PAA_info{node_id}.centroids==unique_centroids(paa_id)).');
    end
    paa_id_init = 0;%paa_id_init+ length(centroid_rep);
    PAA_info{node_id}.centroid_position = unique_PAA_location;
    PAA_info{node_id}.nPAA_centroids = size(PAA_info{node_id}.centroid_position,2);
    PAA_info{node_id}.node_centroid(1:numberTimeDivision, 1:3)  =   squeezeAndReshape(nodeLoc(:, node_id,:));
end

end