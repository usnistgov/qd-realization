function ch_out =  generateChannelPaa(ch_in, infoPAA)
%GENERATECHANNEL_PAA returns the QD channel for each PAA_TX - PAA_RX
%combination.
%
%   [CH_OUT]  =  GENERATECHANNEL_PAA(ch_in, info)
%   **ch_in is the diagonal NxN cell array where N is the number of nodes.
%   Each cell includes struct with fields the channel between PAAs.
%   **info is the supporting structure with PAAs information and indeces
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


nodes = length(infoPAA);
all_nodes = 1:nodes;
ch_out = cell(size(ch_in));
nvar =21;
paa_comb_struct = {};
for nt = all_nodes % Loop on tx nodes
    for nr = all_nodes(all_nodes~=nt)% Loop on rx nodes
        chMIMOtx_rx = []; % Channel between one tx and one rx
        paa_comb = [];
        i =0;
        for c_t = 1:infoPAA{nt}.nPAA_centroids % Loop on transmitter centroid
            for c_r = 1:infoPAA{nr}.nPAA_centroids % Loop on receiver centroid
                chMIMOtx = [];
                paaCombtmp = [];
                frmRotMpInfo = eval(['ch_in{nt,nr}.frmRotMpInfopaaTx', num2str(c_t-1), 'paaRx', num2str(c_r-1),';']);
                for rot_tx = 1:numel(infoPAA{nt}.nodePAAInfo{c_t}.rotated_channel)
                    % Loop on TX PAA genarated with the same centroid ( to
                    % Generate channel obtained by phase rotation)
                    for rot_rx = 1:numel(infoPAA{nr}.nodePAAInfo{c_r}.rotated_channel)
                        % Loop on RX PAA generated with the same centroid
                        ch_siso_tmp = eval(['ch_in{nt,nr}.paaTx', num2str(c_t-1), 'paaRx', num2str(c_r-1),';']);
                        if ~isempty(ch_siso_tmp)
                            ch_siso =ch_siso_tmp(:,:,rot_tx);
                        else
                            ch_siso = [];
                        end
                        
                        % Pointer struct. Indeces of PAAs 
                        ptr.nt = nt;    %TX NODE ID
                        ptr.nr = nr;    %RX NODE ID
                        ptr.paatx = c_t; %TX PAA centroid pointer
                        ptr.paarx = c_r; %RX PAA centroid pointer
                        ptr.rot_tx = rot_tx; %TX PAA rotated channel pointer
                        ptr.rot_rx = rot_rx; %RX PAA rotated channel pointer
                        if isempty(ch_siso)
                            chMIMOtmp_tx = [];
                        else
                            [chMIMOtmp_tx, paaCombtmp] = ddir2MIMO(ch_siso,infoPAA, frmRotMpInfo, ptr);
                        end
                        chMIMOtx = cat(3, chMIMOtx, chMIMOtmp_tx);
                        paa_comb = [paa_comb; paaCombtmp];
%                         size(ch_t)
                    end
                end
                i = i+1;
                chMIMOtx_rx{i} =chMIMOtx; %cat(3, ch_t_r, ch_t);
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

function x = appendNan(x,M,nvar)
if size(x,1)<M
    x(end+1:M,1:nvar,:) = nan;
end
end