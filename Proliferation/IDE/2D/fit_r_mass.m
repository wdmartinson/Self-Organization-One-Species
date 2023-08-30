function  [NLR_fit] =fit_r_mass(r_growth,c_max)
% Load ABM data file 
load('melBoxBirth.mat','Nvals','densityM','meanNumber1DM','numMelBornc','numMelc','spreadAbs','step','stepBnd','stepBndY');
NLR_fit=0;
LENGTH=3.0;% Size of the (domain)^0.5
kmax=5.0;%spatial mesh size
T=150.0;%End time
dt=1.;%timestep
R=0.075; %Radius of ball for computing cell birth/death
mcpd=Nvals;%Number of darts thrown (equivalent to Nbin in the main text)
sxmm=size(densityM{1});
timeStep = ceil(T/sxmm(3));                                 % calculated each timeStep steps
numTsteps = ceil(T/timeStep);                % total number of times evaluated
im_write=ceil(timeStep/dt);%number of ts per save of data for heatmap
c_min=1.0; %Number of cells need in ball for birth to occur
for max_cpd=1:1:length(Nvals)
    L=2^kmax;
    xinterf=linspace(-LENGTH/2,LENGTH/2,2*L+1);
    for i=1:1:length(xinterf)-1
        xmid(i)=(xinterf(i)+xinterf(i+1))/2;%midpoints of the grid x values
        ymid(i)=(xinterf(i)+xinterf(i+1))/2;%midpoints of the grid y values
    end
    for i=1:1:length(xinterf)-1
        for j=1:1:length(xinterf)-1
            rho(i,j)=init_condition_box(xmid(i),ymid(j),kmax,LENGTH);
        end
    end
    Lrho=length(rho);
    dx=xmid(2)-xmid(1);%Mesh size
    dy=ymid(2)-ymid(1);%Mesh size
    t=0;
    k=1;
    tspan=0:im_write*dt:T;
    Nt=length(tspan);
    t_rho=zeros(Lrho(1),Lrho(1),Nt);
    my_i=1;
    lmr=ceil(T/dt);
    mass_rho{max_cpd} = zeros(1,lmr+1);          %Number of cells at each TS
    cells_per_day{max_cpd} = zeros(1,lmr);       %Number of cells born at each TS
    spread{max_cpd} = zeros(1,lmr);
    my_k=1;
        mass_rho{max_cpd}(my_i)=sum(sum(rho))*dx^2;
    for i=1:1:Lrho(1)
        for j=1:1:Lrho(1)
            t_rho(i,j,my_k)=rho(i,j);            %Stores denisty at each TS
        end
    end
        mass_rho{max_cpd}(my_i)=sum(sum(rho))*dx^2; %Uniform mesh so this is \int_\Omega rho
%%%%%%%%% Assembling matrices that store integration weights in Ball
%A{i}{j}(l,k)= dx^2 if (|x_i-x_l|^2+ |x_j-x_k|^2)^0.5<R else A{i}{j}(l,k)=0
for i=1:Lrho(1)
   for j=1:Lrho(1)
       A{i}{j}=zeros(Lrho(1),Lrho(1));
   end
end
L=LENGTH;
for i=1:Lrho(1)
   for j=1:Lrho(1)
       for k=1:Lrho(1)
           for l=1:Lrho(1)
               distx=abs(i-l)*dx;
               if(distx>L/2)
                   distx=(L-distx);
               end
               disty=abs(j-k)*dx;
               if(disty>L/2)
                   disty=(L-disty);
               end
               dist=(disty^2+distx^2)^0.5;
               if(dist<R)
                   A{i}{j}(l,k)=dx^2;
               end
           end
       end
   end
end
%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for t=0:dt:T
        l(my_i)=my_i;
        my_i=my_i+1;
        rhs = rhs_fun(rho,A,c_min,c_max,mcpd(max_cpd),r_growth);
        spread{max_cpd}(my_i-1)=0.;
        for i=1:1:Lrho(1)
            for j=1:1:Lrho(1)
                rho(i,j)=rho(i,j)+dt*rhs(i,j);
                if(rho(i,j)>0.)
                    sijold=spread{max_cpd}(my_i-1);
                    sij=(xmid(i)^2+ymid(j)^2)^0.5;
                    spread{max_cpd}(my_i-1)=max(sijold,sij);
                end
            end
        end
        mass_rho{max_cpd}(my_i)=sum(sum(rho))*dx^2;
        cells_per_day{max_cpd}(my_i-1)=(mass_rho{max_cpd}(my_i)-mass_rho{max_cpd}(my_i-1))/dt;
         if((mod(my_i,im_write)>(im_write-2)))
           my_k=my_k+1;
            for i=1:Lrho(1)
                for j=1:Lrho(1)
                 t_rho(i,j,my_k)=rho(i,j);
                end
            end
         end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function ode_rhs = rhs_fun(rho,A,c_min,c_max,max_cpd,r_growth)%arguments are position, refinement level and domain size
        ode_rhs=zeros(size(rho));
        rhs=zeros(size(rho));
        N=length(rho);
        for i=1:1:N(1)
            for j=1:1:N(1)
                rhs(i,j)=sum(sum(A{i}{j}.*rho));%rhs(i,j)=\int_B_R(x_i,y_j)rho dx
                ode_rhs(i,j)=r_growth*max_cpd*Indicator(c_min,rhs(i,j))*Indicator(rhs(i,j),c_max);
            end
        end
    end
a=size(numMelc{1});
for j=1:1:a(2)
    for max_cpd=1:1:length(Nvals)
        NLR_fit=NLR_fit+(mass_rho{max_cpd}(j+1)-numMelc{max_cpd}(j))^2*step;
    end
end
save('r_share_ODE_mass.mat','r_growth','c_max','mass_rho','spread','cells_per_day');
end

