function pore_data = poreLocation(pore_data)

global nx ny nz n_p
lenX = max(pore_data(:,37));
lenY = max(pore_data(:,38));
lenZ = max(pore_data(:,39));

for ii = 1:n_p
    if rem(ii , nx) == 0
        pore_data(ii,37) = lenX - lenX / nx;
    else
        pore_data(ii,37) = (rem(ii , nx) - 1)*lenX / nx;
    end
end
for jj = 1:n_p
    if pore_data(jj,3) == -1
        pore_data(jj,38) = 0;
    elseif pore_data(jj,3) == -2
        pore_data(jj,38) = (ny + 1)*lenY / ny; %lenY - lenY/ny;
    else
        pore_data(jj,38) = pore_data(jj,3)*lenY / ny;
    end
end
for kk = 1:n_p
    A = floor(kk / (nx*ny));
    pore_data(kk,39) = A*lenZ / nz;
end