function ch_out =  generateChannelPaa(ch_in, infoPAA)
%GENERATECHANNEL_PAA returns the QD channel for each PAA_TX - PAA_RX
%combination.
%
%   [CH_OUT]  =  GENERATECHANNELPAA(ch_in, info)
%   **ch_in is the diagonal NxN cell array where N is the number of nodes.
%   Each cell includes structs with fields:
%      - channel between PAAs.
%      - rotation informations
%   **info is the supporting structure with PAAs information and indices
%   generated in cluster_paa.
%
%   CH_OUT is the diagonal NxN cell array where N is the number of nodes.
%   The (i,j) cell includes the Nray x Nprop x (PAA_TX x PAA_RX) matrix
%   relative to the i_th TX node and j_th RX node. Nray is the number of
%   rays generated in the QD, Nprop is the number of properties and
%   PAA_TX x PAA_RX is the number of channel connecting the i_th TX node
%   and j_th RX node
%
%   Copyright 2019-2020 NIST/CLT (steve.blandino@nist.gov)

%% Input processing 
numNodes = length(infoPAA);
nodesVector = 1:numNodes;
ch_out = cell(size(ch_in));
nvar =21;
paa_comb_struct = {};

%% Generate channel for each PAA
for nt = nodesVector % Loop on tx nodes
    for nr = nodesVector(nodesVector~=nt)% Loop on rx nodes
        chMIMOtx_rx = []; % Channel between one tx and one rx
        paa_comb = [];
        i =0;
        for c_t = 1:infoPAA{nt}.nPAA_centroids % Loop on transmitter centroid
            for c_r = 1:infoPAA{nr}.nPAA_centroids % Loop on receiver centroid
                chMIMOcentroid = []; % Channel between tx and rx centroid
                paaCombtmp = [];
                frmRotMpInfo = eval(['ch_in{nt,nr}.frmRotMpInfopaaTx', num2str(c_t-1), 'paaRx', num2str(c_r-1),';']);
                nIidTx = infoPAA{nt}.nodePAAInfo{c_t}.indep_stoch_channel;%numel(infoPAA{nt}.nodePAAInfo{c_t}.rotated_channel);
                nIidRx = infoPAA{nr}.nodePAAInfo{c_r}.indep_stoch_channel;%numel(infoPAA{nr}.nodePAAInfo{c_r}.rotated_channel)
                for iid_tx = 1:nIidTx
                    % Loop on TX PAA genarated with the same centroid 
                    for iid_rx = 1:nIidRx
                        % Loop on RX PAA generated with the same centroid
                        ch_siso_tmp = eval(['ch_in{nt,nr}.paaTx', num2str(c_t-1), 'paaRx', num2str(c_r-1),';']);
                        if ~isempty(ch_siso_tmp)
                            ch_siso =ch_siso_tmp(:,:,(iid_tx-1)*nIidRx+iid_rx);
                        else
                            ch_siso = [];
                        end
                        
                        % Pointer struct. Indeces of PAAs 
                        ptr.nt = nt;    %TX NODE ID
                        ptr.nr = nr;    %RX NODE ID
                        ptr.paatx = c_t; %TX PAA centroid pointer
                        ptr.paarx = c_r; %RX PAA centroid pointer
                        ptr.iid_tx = iid_tx; %TX PAA rotated channel pointer
                        ptr.iid_rx = iid_rx; %RX PAA rotated channel pointer
                        
                        % Get channel between tx and rx cluster (a
                        % centroid can be the center of different clusters.
                        % Eg cluster 1: PAA generated with the same channel
                        % but different rotation. cluster 2: PAA channels
                        % generated with different stochastic component)
                        if isempty(ch_siso)
                            chMIMOcluster = [];
                        else
                            [chMIMOcluster, paaCombtmp] = ddir2MIMO(ch_siso,infoPAA, frmRotMpInfo, ptr);
                        end
                        chMIMOcentroid = cat(3, chMIMOcentroid, chMIMOcluster);
                        paa_comb = [paa_comb; paaCombtmp];
                    end
                end
                i = i+1;
                chMIMOtx_rx{i} =chMIMOcentroid; %cat(3, ch_t_r, ch_t);
                paa_comb_struct{i} = paa_comb;
            end
        end
        if isempty(chMIMOtx_rx)
            chMIMOtx_rx = [];
        else
            M = max(cellfun(@(x) size(x,1), chMIMOtx_rx));
            chNanPad= cellfun(@(x) appendNan(x,M,nvar),chMIMOtx_rx,'UniformOutput',false);
            chMIMOtx_rx = cat(3,chNanPad{:});
        end
        if ~isempty(paa_comb)
            [~, index_sorted] = sortrows(paa_comb,1);
            ch_out{nt, nr} = chMIMOtx_rx(:,:, index_sorted);
        else
            ch_out{nt, nr} = chMIMOtx_rx;
        end

    end
end

end
%% Append NAN
function x = appendNan(x,M,nvar)
if size(x,1)<M
    x(end+1:M,1:nvar,:) = nan;
end
end