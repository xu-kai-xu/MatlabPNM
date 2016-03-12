function [NumberOfClusters, NodeL, LinkL,cluster_A_nums,cluster_B_nums] =...
    clustering(pore_data,throat_data)

global n_p n_t nx nz

NodeS = zeros(length(pore_data),1);
LinksOfNode = zeros(length(pore_data),26);
NodeNext = zeros(length(pore_data),26);
for i = 1:length(pore_data)
    NodeS(i,1) = pore_data(i,5); % nodes with oil(1),nodes with water(0)
    LinksOfNode(i,1:10) = pore_data(i,41:50);
    for j = 1:length(nonzeros(LinksOfNode(i,1:10)))
        if throat_data(LinksOfNode(i,j),21) == pore_data(i,1)
            NodeNext(i,j) = throat_data(LinksOfNode(i,j),22);
        elseif throat_data(LinksOfNode(i,j),21) ~= pore_data(i,1)
            NodeNext(i,j) = throat_data(LinksOfNode(i,j),21);
        end
    end
end

LinkS = zeros(n_t,1);
for i =1:n_t
    LinkS(i,1) = throat_data(i,5); % throats with oil(1), throats with water(0)
end


OFlag = 1; %oil clusters has numbers 1:NumberOfClusters %water clusters are 0
[NumberOfClusters, NodeL, LinkL] = HKNonLattice(NodeS, LinkS,NodeNext, LinksOfNode, OFlag);

inlet_cluster_indx = zeros(nx*nz,2);
a = 1;
for i = 1:n_p
    if pore_data(i,3) == -1
        inlet_cluster_indx(a,1) = pore_data(i,1);
        inlet_cluster_indx(a,2) = NodeL(i,1);
        a = a + 1;
    end
end

outlet_cluster_indx = zeros(nx*nz,2);
a = 1;
for i = 1:n_p
    if pore_data(i,3) == -2
        outlet_cluster_indx(a,1) = pore_data(i,1);
        outlet_cluster_indx(a,2) = NodeL(i,1);
        a = a + 1;
    end
end

cluster_A_nums = 0;
for i = 1:length(outlet_cluster_indx)
    if outlet_cluster_indx(i,2) ~= 0
        for j = 1:length(inlet_cluster_indx(:,2))
            if outlet_cluster_indx(i,2) == inlet_cluster_indx(j,2);
                if ~any(outlet_cluster_indx(i,2) == cluster_A_nums(:,1))
                    cluster_A_nums(end+1,1) = outlet_cluster_indx(i,2);
                    break
                end
            end
        end
    end
end
cluster_A_nums(1) = [];

cluster_B_nums = 0;
for i = 1:length(outlet_cluster_indx)
    if outlet_cluster_indx(i,2) ~= 0
        if ~any(outlet_cluster_indx(i,2) == cluster_A_nums(:))
%         for j = 1:length(cluster_A_nums)
            if ~any(outlet_cluster_indx(i,2) == cluster_B_nums(:))
                cluster_B_nums(end+1,1) = outlet_cluster_indx(i,2);
            end
%         end
        end
    end
end
cluster_B_nums(1) = [];