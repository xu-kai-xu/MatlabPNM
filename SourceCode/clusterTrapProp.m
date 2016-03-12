function [pore_data, throat_data] = clusterTrapProp(currentPc ,pore_data, throat_data, NodeL, ...
    LinkL, cluster_A_nums , cluster_B_nums)

global n_p n_t

for ii = 1:n_p
    if NodeL(ii) ~= 0 && ~any(NodeL(ii) == cluster_A_nums(:))...
            && ~any(NodeL(ii) == cluster_B_nums(:)) && isnan(pore_data(ii,36)) && pore_data(ii,5) == 1
        pore_data(ii,36) = currentPc;
    end
end

for ii = 1:n_t
    if LinkL(ii) ~= 0 && ~any(LinkL(ii) == cluster_A_nums(:))...
            && ~any(LinkL(ii) == cluster_B_nums(:)) && isnan(throat_data(ii,42)) && throat_data(ii,5) == 1
        throat_data(ii,42) = currentPc;
    end
end

sum(pore_data(:,5));
sum(~isnan(pore_data(:,36)));
a = ~pore_data(:,5) & ~isnan(pore_data(:,36));