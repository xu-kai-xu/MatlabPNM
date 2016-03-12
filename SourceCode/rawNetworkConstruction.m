function [pore_data , throat_data] = rawNetworkConstruction(numX , numY , numZ ,...
    minThroatR , maxThroatR , etaThroatR , deltaThroatR , ...
    minThroatL , maxThroatL , etaThroatL , deltaThroatL , ...
    minAspectRatio , maxAspectRatio , etaAspectRatio , deltaAspectRatio)
%%
global nx ny nz n_p n_t
nx =  numX; ny = numY; nz =  numZ;
n_p = nx * (ny + 2) * nz;
n_t = round((2*nx - 1)*(2*ny + 3)*(2*nz - 1)/2);
%%
r_throat_min = minThroatR;
r_throat_max = maxThroatR;
delta_throat = deltaThroatR;
eta_throat   = etaThroatR ;
x = rand(n_t,1);
r_throat = weibull(r_throat_max,r_throat_min,delta_throat,eta_throat,x);

L_throat_min   = minThroatL;
L_throat_max   = maxThroatL;
delta_L_throat = deltaThroatL;
eta_L_throat   = etaThroatL;
x = rand(n_t,1);
L_throat = weibull(L_throat_max,L_throat_min,delta_L_throat,eta_L_throat,x);
%% pore & throat connectness
% I will imprive it later, to change the coordination number
x = rand(n_p,1);
C_pore = x >= 0; % all of the pores are connected
x = rand(n_t,1);
C_throat = x >= 0; % all of the throats are connected
%% constructing the Raw network
network_model = cell([2*nx - 1,2*ny + 3,2*nz - 1]);
for k = 1:2*nz - 1
    for j = 1:2*ny + 3
        for i = 1:2*nx - 1
            network_model(i,j,k)={[0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ...
               0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]};
        end
    end
end

b=1;
c=1;    %Counter for throat index
for k=1:1:2*nz - 1
    if rem(k,2) ~= 0; % If k is odd
        for j=2:1:2*ny + 2 % choon pore haye inlet nabayad be soorate amoodi be ham vasl bashand
            if j == 2 || j == 2*ny + 2 %j == 2 || j == 4 || j == 2*ny || j == 2*ny + 2
                for i = 1:2:2*nx - 1
                    network_model(i,j,k)={[c 2 0 C_throat(b,1) 0 0 ...
                        0 0 0 0 0 0 r_throat(b,1) L_throat(b,1) 0 0 0 0 0 ...
                        0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]};
                    c = c + 1;
                    b = b + 1; 
                end
            elseif j > 2 && j < 2*ny+2 %j > 4 && j < 2*ny
                if rem(j,2) ~= 0; % If j is odd
                    for i = 2:2:2*nx - 2
                    network_model(i,j,k)={[c 2 0 C_throat(b,1) 0 0 ...
                        0 0 0 0 0 0 r_throat(b,1) L_throat(b,1)...
                        0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]};                        c = c + 1;
                        b = b + 1; 
                    end
                else
                    for i=1:2:2*nx - 1
                    network_model(i,j,k)={[c 2 0 C_throat(b,1) 0 0 ...
                        0 0 0 0 0 0 r_throat(b,1) L_throat(b,1)...
                        0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]};                        c = c + 1;
                        b = b + 1;
                    end
                end
            end
        end
    else
        for j=3:2:2*ny - 1 % choon dar k haye zoj ta khode network asli nabayad throti dashte bashim
            for i=1:2:2*nx-1;
                    network_model(i,j,k)={[c 2 0 C_throat(b,1) 0 0 ...
                        0 0 0 0 0 0 r_throat(b,1) L_throat(b,1)...
                        0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]};                c = c + 1;
                b = b + 1;
            end
        end
    end
end
n_t = b - 1;

%% Pore construction
aspect_ratio_min = minAspectRatio;
aspect_ratio_max = maxAspectRatio;
delta_aspect_ratio = deltaAspectRatio;
eta_aspect_ratio = etaAspectRatio;
x = rand(n_p,1);
aspect_ratio = weibull(aspect_ratio_max,aspect_ratio_min,delta_aspect_ratio,eta_aspect_ratio,x);

r_pore = zeros(n_p,1);
c = 1; % counter for pore index
for k = 1:2:2*nz - 1
    for j = 1:2:2*ny + 3
        for i = 1:2:2*nx - 1
            [Neigh_El] = Neigh_Whole(i,j,k,nx,ny,nz);
            sum_r_Neigh_throats = 0;
            max_r_Neigh_throats = 0;
            for ii = 1:length(Neigh_El.x)
                sum_r_Neigh_throats = sum_r_Neigh_throats + ...
                    network_model{Neigh_El.x(ii),Neigh_El.y(ii),Neigh_El.z(ii)}(13);
                if network_model{Neigh_El.x(ii),Neigh_El.y(ii),Neigh_El.z(ii)}(13) > max_r_Neigh_throats
                    max_r_Neigh_throats = network_model{Neigh_El.x(ii),Neigh_El.y(ii),Neigh_El.z(ii)}(13);
                end
            end
                r_pore(c,1) = max(aspect_ratio(c,1).*sum_r_Neigh_throats/length(Neigh_El.x),...
                    max_r_Neigh_throats);
                c = c + 1;
        end
    end
end
%% Detrmination of the pore Length
L_pore = 2 * r_pore;
%%
c = 1;
for k=1:2:2*nz-1
    for j=1:2:2*ny+3
        for i=1:2:2*nx-1
            if j == 1
                a = -1;
            elseif j == 2*ny+3
                a = -2;
            else 
                a = floor(j/2) ; % J layer is in the middle region
            end
            network_model(i,j,k)={[c 1 a C_pore(c,1) 0 0 ...
                0 0 0 0 0 r_pore(c,1) L_pore(c,1) 0 0 ...
                0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]};
            c = c + 1; 
            b = b + 1;
        end
    end
end
%% Finding the pore indexes connected to each throat
for k = 1:2*nz-1
    for j = 1:2*ny+3
        for i = 1:2*nx-1
            if network_model{i,j,k}(2) == 2
                throat_neigh = Neigh_Whole(i,j,k,nx,ny,nz);
                for ii = 1:length(throat_neigh.x)
                    if network_model{throat_neigh.x(ii),...
                            throat_neigh.y(ii),throat_neigh.z(ii)}(4) == 1
                        if isempty(~nonzeros(network_model{i,j,k}(21)))
                        network_model{i,j,k}(21) = network_model...
                            {throat_neigh.x(ii),throat_neigh.y(ii),throat_neigh.z(ii)}(1);
                        network_model{i,j,k}(23) = network_model...
                            {throat_neigh.x(ii),throat_neigh.y(ii),throat_neigh.z(ii)}(13);
                        elseif isempty(~nonzeros(network_model{i,j,k}(22)))
                        network_model{i,j,k}(22) = network_model...
                            {throat_neigh.x(ii),throat_neigh.y(ii),throat_neigh.z(ii)}(1);
                        network_model{i,j,k}(24) = network_model...
                            {throat_neigh.x(ii),throat_neigh.y(ii),throat_neigh.z(ii)}(13);
                        end
                    end
                end
                network_model{i,j,k}(20) = network_model{i,j,k}(14) + ...
                    network_model{i,j,k}(23) + network_model{i,j,k}(24);
            end
        end
    end
end
%% Characterizing the connecting link indexes
for k = 1:2*nz - 1
    for j = 1:2*ny + 3
        for i = 1:2*nx - 1
            if network_model{i,j,k}(4) == 1
                if rem(i,2) == 1 && rem(j,2) == 1 && rem(k,2) == 1 
                    throats_conc_pore = Neigh_Whole(i,j,k,nx,ny,nz);
                    jj = 41;
                    for ii = 1:length(throats_conc_pore.x)
                        if network_model{throats_conc_pore.x(ii),throats_conc_pore.y(ii),...
                                throats_conc_pore.z(ii)}(4) == 1
                        network_model{i,j,k}(jj) = network_model{throats_conc_pore.x(ii),...
                            throats_conc_pore.y(ii),throats_conc_pore.z(ii)}(1);
                        jj = jj + 1;
                        end
                    end
                end
            end
        end
    end
end
%% Characterizing pore positions
for k = 1:2*nz - 1
    for j = 1:2*ny + 3
        for i = 1:2*nx - 1
            if network_model{i,j,k}(4) == 1 && rem(i,2) == 1 && ...
                    rem(j,2) == 1 && rem(k,2) == 1 
                if i == 1
                    network_model{i,j,k}(37) = 0;
                else
                    if j == 1 || j == 2*ny + 1
                        network_model{i,j,k}(37) = network_model{i-2,j,k}(37) + ...
                            mean([minThroatL,maxThroatL]) + ...
                            0.5 * (network_model{i,j,k}(13) + network_model{i-2,j,k}(13));
%                             0.5 * (network_model{i,j,k}(13) + network_model{i-2,j,k}(13));
                    else
                    network_model{i,j,k}(37) = network_model{i-2,j,k}(37) +...
                        network_model{i-1,j,k}(14) + ...
                        (network_model{i,j,k}(13) + network_model{i-2,j,k}(13));
%                         0.5 * (network_model{i,j,k}(13) +
%                         network_model{i-2,j,k}(13));
                    end
                end
                if j == 1
                    network_model{i,j,k}(38) = 0;
                else
                    network_model{i,j,k}(38) = network_model{i,j - 2,k}(38) + ...
                        network_model{i,j-1,k}(14) + ...
                        (network_model{i,j,k}(13) + network_model{i,j-2,k}(13));
%                         0.5 * (network_model{i,j,k}(13) + network_model{i,j-2,k}(13));
                end                
                
                if k == 1
                    network_model{i,j,k}(39) = 0;
                else
                    network_model{i,j,k}(39) = network_model{i,j,k - 2}(39) +...
                        network_model{i,j,k-1}(14) + ...
                      (network_model{i,j,k}(13) + network_model{i,j,k-2}(13));  
%                     0.5 * (network_model{i,j,k}(13) + network_model{i,j,k-2}(13));
                end
                
            end
        end
    end
end

%%
p = 1;
t = 1;
pore_data = zeros(n_p,50);
throat_data = zeros(n_t,50);

for k = 1:2*nz-1
    for j = 1:2*ny+3
        for i = 1:2*nx-1
            if network_model{i,j,k}(2) == 1  
                pore_data(p,1) = network_model{i,j,k}(1); % pore index
                pore_data(p,2) = network_model{i,j,k}(2); % pore or throat type
                pore_data(p,3) = network_model{i,j,k}(3); % J layer number
                pore_data(p,4) = network_model{i,j,k}(4); % connectness
                pore_data(p,5) = network_model{i,j,k}(5); % occupied by water(0) or oil(1)
                pore_data(p,6) = network_model{i,j,k}(6); % G factor
                pore_data(p,7) = network_model{i,j,k}(7); % Geometry,cir(1),tri(2),squ(3)
                pore_data(p,8) = network_model{i,j,k}(8); % half angle
                pore_data(p,12) = network_model{i,j,k}(12); % Inscribed raduis
                pore_data(p,13) = network_model{i,j,k}(13); % Length
                pore_data(p,14) = network_model{i,j,k}(14); % Cross section area
                pore_data(p,15) = network_model{i,j,k}(15); % Volume
                pore_data(p,16) = network_model{i,j,k}(16); % Receeding angle
                pore_data(p,17) = network_model{i,j,k}(17); % Advancing angle
                pore_data(p,23) = network_model{i,j,k}(23); % water volume in the pore
                pore_data(p,25) = network_model{i,j,k}(25); % single phase conductance
                pore_data(p,26) = network_model{i,j,k}(26); % water conductance
                pore_data(p,37) = network_model{i,j,k}(37); % X pore position
                pore_data(p,38) = network_model{i,j,k}(38); % Y pore position
                pore_data(p,39) = network_model{i,j,k}(39); % Z pore position
                pore_data(p,41) = network_model{i,j,k}(41); % here to the end "throat
                pore_data(p,42) = network_model{i,j,k}(42); % indexes connected to the pore"
                pore_data(p,43) = network_model{i,j,k}(43);
                pore_data(p,44) = network_model{i,j,k}(44);
                pore_data(p,45) = network_model{i,j,k}(45);
                pore_data(p,46) = network_model{i,j,k}(46);
                p = p + 1;
                
            elseif network_model{i,j,k}(2) == 2 
                throat_data(t,1) = network_model{i,j,k}(1); % throat index
                throat_data(t,2) = network_model{i,j,k}(2); % pore or throat type
                throat_data(t,4) = network_model{i,j,k}(4); % connectness
                throat_data(t,5) = network_model{i,j,k}(5); % occupied by water(0) or oil(1) 
                throat_data(t,6) = network_model{i,j,k}(6); % G factor
                throat_data(t,7) = network_model{i,j,k}(7); % Geometry,cir(1),tri(2),squ(3)
                throat_data(t,8) = network_model{i,j,k}(8); % half angle
                throat_data(t,13) = network_model{i,j,k}(13); % Inscribed radius
                throat_data(t,14) = network_model{i,j,k}(14); % Length
                throat_data(t,15) = network_model{i,j,k}(15); % cross section area
                throat_data(t,16) = network_model{i,j,k}(16); % volume
                throat_data(t,18) = network_model{i,j,k}(18); % Receeding angle
                throat_data(t,19) = network_model{i,j,k}(19); % Advancing angle
                throat_data(t,20) = network_model{i,j,k}(20); % throat total length
                throat_data(t,21) = network_model{i,j,k}(21); % pore 1 index
                throat_data(t,22) = network_model{i,j,k}(22); % pore 2 index
                throat_data(t,23) = network_model{i,j,k}(23); % length pore 1
                throat_data(t,24) = network_model{i,j,k}(24); % length pore 2
                throat_data(t,29) = network_model{i,j,k}(29); % water volume
                throat_data(t,31) = network_model{i,j,k}(31); % single phase conductance
                throat_data(t,32) = network_model{i,j,k}(32); % water conductance
                t = t + 1;                
            end
        end
    end
end