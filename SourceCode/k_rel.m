function [oil_rel_perm,Water_rel_perm] = k_rel(p_in,p_out)

global pore_data throat_data nz ny nx n_t
global q_tot inletLayer outletLayer


spec_pore_data = zeros(nx*ny*nz,size(pore_data,2));
a = 1;
for i = 1:length(pore_data(:,1))
    if pore_data(i,3) > 0
        spec_pore_data(a,1:end) = pore_data(i,1:end);
        a = a + 1;
    end
end

A_water = zeros(nx*ny*nz,nx*ny*nz);
A_oil = zeros(nx*ny*nz,nx*ny*nz);
for i = 1:nx*ny*nz
    for j = 1:nx*ny*nz
        if i == j
            throat_index_of_pore = nonzeros(spec_pore_data(i,41:end));
            for ii = 1:length(throat_index_of_pore)
                if throat_data(throat_index_of_pore(ii),21) == spec_pore_data(i,1)
                    next_pore = throat_data(throat_index_of_pore(ii),22);
                else
                    next_pore = throat_data(throat_index_of_pore(ii),21);
                end
                next_pore = find(pore_data(:,1) == next_pore);
                if pore_data(next_pore,3) == -1 || pore_data(next_pore,3) == -2
                    A_water(i,j) = A_water(i,j) + (throat_data(throat_index_of_pore(ii),14) /...
                        throat_data(throat_index_of_pore(ii),31) + ...
                        0.5*(spec_pore_data(i,13)/spec_pore_data(i,25) +...
                        pore_data(next_pore,13)/pore_data(next_pore,25))).^-1;
                    A_oil(i,j) = A_oil(i,j) + (throat_data(throat_index_of_pore(ii),14) /...
                        throat_data(throat_index_of_pore(ii),31) + ...
                        0.5*(spec_pore_data(i,13)/spec_pore_data(i,25) +...
                        pore_data(next_pore,13)/pore_data(next_pore,25))).^-1;
                else
                A_water(i,j) = A_water(i,j) + (throat_data(throat_index_of_pore(ii),14) /...
                    throat_data(throat_index_of_pore(ii),32) + ...
                    0.5*(spec_pore_data(i,13)/spec_pore_data(i,26) +...
                    pore_data(next_pore,13)/pore_data(next_pore,26))).^-1;
                A_oil(i,j) = A_oil(i,j) + (throat_data(throat_index_of_pore(ii),14) /...
                    throat_data(throat_index_of_pore(ii),33) + ...
                    0.5*(spec_pore_data(i,13)/spec_pore_data(i,27) +...
                    pore_data(next_pore,13)/pore_data(next_pore,27))).^-1;
                end
            end
        elseif i ~= j
            middle_throat_index = find(throat_data(:,21) == spec_pore_data(i,1)...
                & throat_data(:,22) == spec_pore_data(j,1) |...
                throat_data(:,21) == spec_pore_data(j,1)...
                & throat_data(:,22) == spec_pore_data(i,1) );
            if ~isempty(middle_throat_index)
                A_water(i,j) = -(throat_data(middle_throat_index,14) /...
                    throat_data(middle_throat_index,32) +...
                    0.5*(spec_pore_data(i,13)/spec_pore_data(i,26) +...
                    spec_pore_data(j,13)/spec_pore_data(j,26))).^-1;
                A_oil(i,j) = -(throat_data(middle_throat_index,14) /...
                    throat_data(middle_throat_index,33) +...
                    0.5*(spec_pore_data(i,13)/spec_pore_data(i,27) +...
                    spec_pore_data(j,13)/spec_pore_data(j,27))).^-1;
            end 
        end
    end
end

D_water = zeros(nx*ny*nz,1);
D_oil = zeros(nx*ny*nz,1);
for i = 1:nx*ny*nz
    throat_neigh_spec_pore_index = nonzeros(spec_pore_data(i,41:end));
    for ii = 1:length(throat_neigh_spec_pore_index)
        if throat_data(throat_neigh_spec_pore_index(ii),21) == spec_pore_data(i,1)
            next_pore = throat_data(throat_neigh_spec_pore_index(ii),22);
        else
            next_pore = throat_data(throat_neigh_spec_pore_index(ii),21);
        end
        if pore_data(next_pore,3) == -1 
            inlet_pore = next_pore;
            D_water(i,1) = p_in*((throat_data(throat_neigh_spec_pore_index(ii),14)/...
                throat_data(throat_neigh_spec_pore_index(ii),31))+ ...
                0.5*((pore_data(inlet_pore,13)/pore_data(inlet_pore,25))+...
                (spec_pore_data(i,13)/spec_pore_data(i,25))))^-1;
            D_oil(i,1) = p_in*((throat_data(throat_neigh_spec_pore_index(ii),14)/...
                throat_data(throat_neigh_spec_pore_index(ii),31))+ ...
                0.5*((pore_data(inlet_pore,13)/pore_data(inlet_pore,25))+...
                (spec_pore_data(i,13)/spec_pore_data(i,25))))^-1;

        elseif pore_data(next_pore,3) == -2
            outlet_pore = next_pore;
            D_water(i,1) = p_out*((throat_data(throat_neigh_spec_pore_index(ii),14)/...
                throat_data(throat_neigh_spec_pore_index(ii),31))+ ...
                0.5*((pore_data(outlet_pore,13)/pore_data(outlet_pore,25))+...
                (spec_pore_data(i,13)/spec_pore_data(i,25))))^-1;
            D_oil(i,1) = p_out*((throat_data(throat_neigh_spec_pore_index(ii),14)/...
                throat_data(throat_neigh_spec_pore_index(ii),31))+ ...
                0.5*((pore_data(outlet_pore,13)/pore_data(outlet_pore,25))+...
                (spec_pore_data(i,13)/spec_pore_data(i,25))))^-1;           
        end
    end
end

[pore_press_water,flag] = pcg(A_water,D_water,1e-10,30000);
clear flag
[pore_press_oil,flag] = pcg(A_oil,D_oil,1e-10,30000);
clear flag

a = 1;
for i =1:nx*(ny+2)*nz
    if pore_data(i,3) == -1
        pore_data(i,19) = p_in;
        pore_data(i,20) = p_in;
    elseif pore_data(i,3) == -2
        pore_data(i,19) = p_out;
        pore_data(i,20) = p_out;
    else
        pore_data(i,19) = pore_press_water(a);
        pore_data(i,20) = pore_press_oil(a);
        a = a + 1;
    end
end

pres_av_w = zeros(1,ny);
for ii = 1:ny
    cntr = pore_data(:,3) == ii;
    pres_av_w(ii) = sum(pore_data(cntr,19))/length(nonzeros(cntr));      
end
pres_av_o = zeros(1,ny);
for ii = 1:ny
    cntr = pore_data(:,3) == ii;
    pres_av_o(ii) = sum(pore_data(cntr,20))/length(nonzeros(cntr));      
end


q_tot_w = 0;
q_tot_o = 0;
for i = 1:n_t
    if pore_data((throat_data(i,21)),3) == inletLayer && pore_data((throat_data(i,22)),3) == inletLayer + 1 || ...
            pore_data((throat_data(i,21)),3) == inletLayer + 1 && pore_data((throat_data(i,22)),3) == inletLayer
        
        q_tot_w = q_tot_w + (pore_data((throat_data(i,21)),19) - pore_data((throat_data(i,22)),19))*...
            ((throat_data(i,14)/throat_data(i,32))+0.5*((pore_data((throat_data(i,21)),13)/...
            (pore_data((throat_data(i,21)),26)))+(pore_data((throat_data(i,22)),13)/...
            pore_data((throat_data(i,22)),26))))^-1;
        q_tot_o = q_tot_o + (pore_data((throat_data(i,21)),20) - pore_data((throat_data(i,22)),20))*...
            ((throat_data(i,14)/throat_data(i,33))+0.5*((pore_data((throat_data(i,21)),13)/...
            (pore_data((throat_data(i,21)),27)))+(pore_data((throat_data(i,22)),13)/...
            pore_data((throat_data(i,22)),27))))^-1;
        
%     elseif pore_data((throat_data(i,21)),3) == outletLayer && pore_data((throat_data(i,22)),3) == -2 || ...
%             pore_data((throat_data(i,21)),3) == -2 && pore_data((throat_data(i,22)),3) == outletLayer
%         
%         q_tot_w = q_tot_w + (pore_data((throat_data(i,21)),19) - pore_data((throat_data(i,22)),19))*...
%             ((throat_data(i,14)/throat_data(i,32))+0.5*((pore_data((throat_data(i,21)),13)/...
%             (pore_data((throat_data(i,21)),26)))+(pore_data((throat_data(i,22)),13)/...
%             pore_data((throat_data(i,22)),26))))^-1;
%         q_tot_o = q_tot_o + (pore_data((throat_data(i,21)),20) - pore_data((throat_data(i,22)),20))*...
%             ((throat_data(i,14)/throat_data(i,33))+0.5*((pore_data((throat_data(i,21)),13)/...
%             (pore_data((throat_data(i,21)),27)))+(pore_data((throat_data(i,22)),13)/...
%             pore_data((throat_data(i,22)),27))))^-1;
       
    end
end
% q_tot_w = q_tot_w / 2;
% q_tot_o = q_tot_o / 2;

oil_rel_perm = q_tot_o / q_tot;
if oil_rel_perm > 1
    oil_rel_perm = 1;
elseif oil_rel_perm < 0 
    oil_rel_perm  = 0;
end

Water_rel_perm = q_tot_w / q_tot;
if Water_rel_perm > 1
    Water_rel_perm = 1;
elseif Water_rel_perm < 0 
    Water_rel_perm  = 0;
end