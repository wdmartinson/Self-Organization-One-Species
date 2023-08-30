function [uinit] = init_condition(s,r,k,L)%arguments are positionx, positiony, refinement level and domain size
ms=L/(2^(k+1)+1);%Mesh size
xc=max(abs(s),abs(r)); % Find the maximum position of the cells in the x / y direction
uinit =0.; 
    if(xc<0.5)
        uinit =400;        % Create a 1mm x 1mm box of 400 melanophore cells centred at the origin
    end
end
