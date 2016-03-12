function [Neigh_El] = Neigh_Whole(x,y,z,nx,ny,nz)

if z == 1
%% z=1
    if y == 1
        if x == 1
            Neigh_El.x = [x+1 , x , x];
            Neigh_El.y = [y , y+1 , y];
            Neigh_El.z = [z , z , z+1];
        elseif x == 2*nx-1
            Neigh_El.x = [x-1 , x , x];
            Neigh_El.y = [y , y+1 , y];
            Neigh_El.z = [z , z , z+1];
        else
            Neigh_El.x = [x-1 , x+1 , x , x];
            Neigh_El.y = [y , y , y+1 , y];
            Neigh_El.z = [z , z , z , z+1];
        end
    elseif y == 2*ny+3
        
        if x == 1
            Neigh_El.x = [x+1 , x , x];
            Neigh_El.y = [y , y-1 , y];
            Neigh_El.z = [z , z , z+1];
        elseif x == 2*nx-1
            Neigh_El.x = [x-1 , x , x];
            Neigh_El.y = [y , y-1 , y];
            Neigh_El.z = [z , z , z+1];
        else
            Neigh_El.x = [x-1 , x+1 , x , x];
            Neigh_El.y = [y , y , y-1 , y];
            Neigh_El.z = [z , z , z , z+1];
        end
    else % 2<y<2n-2
        
        if x == 1
            Neigh_El.x = [x+1 , x , x ,x];
            Neigh_El.y = [y , y-1 , y+1 , y];
            Neigh_El.z = [z , z , z , z+1];
        elseif x == 2*nx-1
            Neigh_El.x = [x-1 , x , x ,x];
            Neigh_El.y = [y , y-1 , y+1 , y];
            Neigh_El.z = [z , z , z , z+1];
        else
            Neigh_El.x = [x-1 , x+1 , x , x , x];
            Neigh_El.y = [y , y , y-1 , y+1 , y];
            Neigh_El.z = [z , z , z, z , z+1];
        end
    end
 
%%
elseif z == 2*nz-1
%% z=2n-1
    if y == 1
        if x == 1
            Neigh_El.x = [x+1 , x , x];
            Neigh_El.y = [y , y+1 , y];
            Neigh_El.z = [z , z , z-1];
        elseif x == 2*nx-1
            Neigh_El.x = [x-1 , x , x];
            Neigh_El.y = [y , y+1 , y];
            Neigh_El.z = [z , z , z-1];
        else
            Neigh_El.x = [x-1 , x+1 , x , x];
            Neigh_El.y = [y , y , y+1 , y];
            Neigh_El.z = [z , z , z , z-1];
        end
    elseif y == 2*ny+3
        
        if x == 1
            Neigh_El.x = [x+1 , x , x];
            Neigh_El.y = [y , y-1 , y];
            Neigh_El.z = [z , z , z-1];
        elseif x == 2*nx-1
            Neigh_El.x = [x-1 , x , x];
            Neigh_El.y = [y , y-1 , y];
            Neigh_El.z = [z , z , z-1];
        else
            Neigh_El.x = [x-1 , x+1 , x , x];
            Neigh_El.y = [y , y , y-1 , y];
            Neigh_El.z = [z , z , z , z-1];
        end
    else % 2<y<2n-2
        
        if x == 1
            Neigh_El.x = [x+1 , x , x ,x];
            Neigh_El.y = [y , y-1 , y+1 , y];
            Neigh_El.z = [z , z , z , z-1];
        elseif x == 2*nx-1
            Neigh_El.x = [x-1 , x , x ,x];
            Neigh_El.y = [y , y-1 , y+1 , y];
            Neigh_El.z = [z , z , z , z-1];
        else
            Neigh_El.x = [x-1 , x+1 , x , x , x];
            Neigh_El.y = [y , y , y-1 , y+1 , y];
            Neigh_El.z = [z , z , z, z , z-1];
        end
    end
%%    
else 
    %% 2<z<2n-2
    if y == 1
        if x == 1
            Neigh_El.x = [x+1 , x , x ,x];
            Neigh_El.y = [y , y+1 , y ,y];
            Neigh_El.z = [z , z , z-1 ,z+1];
        elseif x == 2*nx-1
            Neigh_El.x = [x-1 , x , x, x];
            Neigh_El.y = [y , y+1 , y, y];
            Neigh_El.z = [z , z , z-1, z+1];
        else
            Neigh_El.x = [x-1 , x+1 , x , x, x];
            Neigh_El.y = [y , y , y+1 , y, y];
            Neigh_El.z = [z , z , z , z-1, z+1];
        end
    elseif y == 2*ny+3
        
        if x == 1
            Neigh_El.x = [x+1 , x , x, x];
            Neigh_El.y = [y , y-1 , y , y];
            Neigh_El.z = [z , z , z-1 , z+1];
        elseif x == 2*nx-1
            Neigh_El.x = [x-1 , x , x, x];
            Neigh_El.y = [y , y-1 , y , y];
            Neigh_El.z = [z , z , z-1 , z+1];
        else
            Neigh_El.x = [x-1 , x+1 , x , x , x];
            Neigh_El.y = [y , y , y-1 , y , y];
            Neigh_El.z = [z , z , z , z-1 , z+1];
        end
    else % 2<y<2n-2
        
        if x == 1
            Neigh_El.x = [x+1 , x , x ,x, x];
            Neigh_El.y = [y , y-1 , y+1 , y ,y];
            Neigh_El.z = [z , z , z , z-1, z+1];
        elseif x == 2*nx-1
            Neigh_El.x = [x-1 , x , x ,x,x];
            Neigh_El.y = [y , y-1 , y+1 , y , y];
            Neigh_El.z = [z , z , z , z-1 , z+1];
        else
            Neigh_El.x = [x-1 , x+1 , x , x , x ,x];
            Neigh_El.y = [y , y , y-1 , y+1 , y ,y];
            Neigh_El.z = [z , z , z, z , z-1 ,z+1];
        end
    end
end