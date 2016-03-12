function absolute_permeability = k_abs(p_in ,p_out)

global pore_data throat_data nz ny nx water_viscosity n_t
global q_tot inletLayer outletLayer
global modelLength modelCrossArea

inletLayer = floor(.5*ny)+1;
outletLayer = ny;
modelLength = max(pore_data(:,37));
modelCrossArea = max(pore_data(:,38)) * max(pore_data(end,39));

spec_pore_data = zeros(nx*ny*nz,size(pore_data,2));
a = 1;
for i = 1:length(pore_data(:,1))
    if pore_data(i,3) > 0
        spec_pore_data(a,1:end) = pore_data(i,1:end);
        a = a + 1;
    end
end


A = zeros(nx*ny*nz,nx*ny*nz);
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
                    A(i,j) = A(i,j) + (throat_data(throat_index_of_pore(ii),14) /...
                        throat_data(throat_index_of_pore(ii),31) + ...
                        0.5*(spec_pore_data(i,13)/spec_pore_data(i,25) +...
                        pore_data(next_pore,13)/pore_data(next_pore,25))).^-1;
                else
                A(i,j) = A(i,j) + (throat_data(throat_index_of_pore(ii),14) /...
                    throat_data(throat_index_of_pore(ii),31) + ...
                    0.5*(spec_pore_data(i,13)/spec_pore_data(i,25) +...
                    pore_data(next_pore,13)/pore_data(next_pore,25))).^-1;
                end
            end
        elseif i ~= j
            middle_throat_index = find(throat_data(:,21) == spec_pore_data(i,1)...
                & throat_data(:,22) == spec_pore_data(j,1) |...
                throat_data(:,21) == spec_pore_data(j,1)...
                & throat_data(:,22) == spec_pore_data(i,1) );
            if ~isempty(middle_throat_index)
                A(i,j) = -(throat_data(middle_throat_index,14) /...
                    throat_data(middle_throat_index,31) +...
                    0.5*(spec_pore_data(i,13)/spec_pore_data(i,25) +...
                    spec_pore_data(j,13)/spec_pore_data(j,25))).^-1;
            end 
        end
    end
end


D = zeros(nx*ny*nz,1);
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
            D(i,1) = p_in*((throat_data(throat_neigh_spec_pore_index(ii),14)/...
                throat_data(throat_neigh_spec_pore_index(ii),31))+ ...
                0.5*((pore_data(inlet_pore,13)/pore_data(inlet_pore,25))+...
                (spec_pore_data(i,13)/spec_pore_data(i,25))))^-1;
        elseif pore_data(next_pore,3) == -2
            outlet_pore = next_pore;
            D(i,1) = p_out*((throat_data(throat_neigh_spec_pore_index(ii),14)/...
                throat_data(throat_neigh_spec_pore_index(ii),31))+ ...
                0.5*((pore_data(outlet_pore,13)/pore_data(outlet_pore,25))+...
                (spec_pore_data(i,13)/spec_pore_data(i,25))))^-1;
            
        end
    end
end

% [pore_press,flag] = pcg(A,D,1e-10,1000);
% clear flag;
% [pore_press,flag] = minres(A,D,1e-10,1000);
[pore_press,flag] = pcg(A,D,1e-10,30000);
clear flag;
pore_press(pore_press > p_in) = p_in;

    
a = 1;
for i =1:nx*(ny+2)*nz
    if pore_data(i,3) == -1
        pore_data(i,19) = p_in;
    elseif pore_data(i,3) == -2
        pore_data(i,19) = p_out;
    else
        pore_data(i,19) = pore_press(a);
        a = a + 1;
    end
end
%
pres_av = zeros(1,ny);
for ii = 1:ny
    if ii == 0
        continue
    end
    cntr = pore_data(:,3) == ii;
    pres_av(ii) = sum(pore_data(cntr,19))/length(nonzeros(cntr));      
end
%% total flow rate calculation
q_tot = 0;
L = 0;
A = 0;
c = 0;
% P_A_in = 0;
% P_A_out = 0;
% A_in = 0;
% A_out = 0;
% layerNumber = outletLayer - inletLayer + 1;
for i = 1:n_t
%     for jj = inletLayer:outletLayer
%         if jj == ny
%             outlay = -2;
%         else
%             outlay = jj + 1;
%         end
        if pore_data((throat_data(i,21)),3) == inletLayer && pore_data((throat_data(i,22)),3) == inletLayer+1 || ...
                pore_data((throat_data(i,21)),3) == inletLayer+1 && pore_data((throat_data(i,22)),3) == inletLayer
        
            L = L + throat_data(i,20);
            A = A + throat_data(i,15); 
%         P_A_in = P_A_in + pore_data((throat_data(i,21)),19)*pore_data((throat_data(i,21)),14);
%         A_in = A_in + pore_data((throat_data(i,21)),14);
%         P_A_out = P_A_out + pore_data((throat_data(i,22)),19)*pore_data((throat_data(i,22)),14);
%         A_out = A_out + pore_data((throat_data(i,22)),14);
        
            q_tot = q_tot + (pore_data((throat_data(i,21)),19) - pore_data((throat_data(i,22)),19))*...
                ((throat_data(i,14)/throat_data(i,31))+0.5*((pore_data((throat_data(i,21)),13)/...
                (pore_data((throat_data(i,21)),25)))+(pore_data((throat_data(i,22)),13)/...
                pore_data((throat_data(i,22)),25))))^-1;
            
%         elseif pore_data((throat_data(i,21)),3) == outletLayer && pore_data((throat_data(i,22)),3) == -2 || ...
%                 pore_data((throat_data(i,21)),3) == -2 && pore_data((throat_data(i,22)),3) == outletLayer
%             q_tot = q_tot + (pore_data((throat_data(i,21)),19) - pore_data((throat_data(i,22)),19))*...
%                 ((throat_data(i,14)/throat_data(i,31))+0.5*((pore_data((throat_data(i,21)),13)/...
%                 (pore_data((throat_data(i,21)),25)))+(pore_data((throat_data(i,22)),13)/...
%                 pore_data((throat_data(i,22)),25))))^-1;
            c = c + 1;
        end 
%     end
end
% q_tot = q_tot / layerNumber;
% q_tot = q_tot /2;

L = L/c;
A = A/c;
% P_in = P_A_in / A_in;
% P_out = P_A_out / A_out;
%%
absolute_permeability = water_viscosity*q_tot*modelLength/modelCrossArea/(p_in - p_out); %(m2)