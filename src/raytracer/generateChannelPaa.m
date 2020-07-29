function ch_out =  generateChannelPaa(ch_in, info)
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


nodes = length(info);
all_nodes = 1:nodes;
ch_out = cell(size(ch_in));
nvar =21;
for nt = all_nodes
    for nr = all_nodes(all_nodes~=nt)
        ch_t_r = [];
        i =0;
        for c_t = 1:info{nt}.nPAA_centroids
            for c_r = 1:info{nr}.nPAA_centroids
                ch_t = [];
                for rot_tx = 1:numel(info{nt}.nodePAAInfo{c_t}.rotated_channel)
                    for rot_rx = 1:numel(info{nr}.nodePAAInfo{c_r}.rotated_channel)
                        t = eval(['ch_in{nt,nr}.paaTx', num2str(c_t-1), 'paaRx', num2str(c_r-1),'(:,:,rot_tx);']);
                        ptr.nt = nt;
                        ptr.nr = nr;
                        ptr.paatx = c_t;
                        ptr.paarx = c_r;
                        ptr.rot_tx = rot_tx;
                        ptr.rot_rx = rot_rx;
                        if isempty(t)
                            ch_t_tmp = [];
                        else
                            ch_t_tmp = ddir2MIMO(t,info, ptr);
                        end
                        ch_t = cat(3, ch_t, ch_t_tmp);
%                         size(ch_t)
                    end
                end
                i = i+1;
                ch_t_r{i} =ch_t; %cat(3, ch_t_r, ch_t);
            end
        end
        if isempty(ch_t_r)
            ch_t_r = [];
        else
            M = max(cellfun(@(x) size(x,1), ch_t_r));
            gg= cellfun(@(x) appendNan(x,M,nvar),ch_t_r,'UniformOutput',false);
            ch_t_r = cat(3,gg{:});
        end
        ch_out{nt, nr} = ch_t_r;
    end
end

end

function x = appendNan(x,M,nvar)
if size(x,1)<M
    x(end+1:M,1:nvar,:) = nan;
end
end