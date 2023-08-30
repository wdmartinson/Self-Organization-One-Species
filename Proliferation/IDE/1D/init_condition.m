function [uinit] = init_condition(s,k,L)%arguments are position, refinement level and domain size
ms=L/(2^(k+1)+1);%Mesh size
xc=s;%centre domain at 0.
uinit =0.; 
    if((abs(xc-0.)<ms))
        uinit =1.0/(ms);        
    end
    %uinit=20.*exp(-(xc*xc)/0.01);
end
