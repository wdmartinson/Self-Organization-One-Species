function NLR_fit = fit_r_full_mel(r_growth,c_max)
% Load ABM data file
load('melanophores1D_ABM.mat','meanNumMel','meanNumMelBorn','spreadAbs','timeStep','melCellsPer_mm')
NLR_fit=0;
LENGTH=3.0;% Size of the domain
kmax=6.0;%spatial mesh size
T=2000.0;%End time
dt=1.;%timestep
max_cpd=1; %Maximum number of cells nucleated per day
R=0.075; %Radius of ball for computing cell birth/death
c_min=1.0; %Number of cells need in ball for birth to occur
im_write=1000;%number of ts per save of data
for max_cpd=1:1:10
L=2^kmax;
xinterf=linspace(-LENGTH/2,LENGTH/2,2*L+1);
for i=1:1:length(xinterf)-1
xmid(i)=(xinterf(i)+xinterf(i+1))/2;
rho(i)=init_condition(xmid(i),kmax,LENGTH);
W11(i)=Indicator(xmid(i),R);
end
N=length(rho);
dx=xmid(2)-xmid(1);
t=0;
k=1;
tspan=0:im_write*dt:T;
N=length(rho);
Nt=length(tspan);
t_rho=zeros(N,Nt);
my_i=1;
mass_rho(max_cpd,my_i)=sum(rho)*dx;
for i=1:1:length(rho)
    t_rho(i,k)=rho(i);
end
mass_rho(max_cpd,my_i)=sum(rho)*dx;
for t=0:dt:T
    l(my_i)=my_i;
    my_i=my_i+1;
    mass_rho(max_cpd,my_i)=sum(rho)*dx;
    cells_per_day(max_cpd,my_i-1)=(mass_rho(max_cpd,my_i)-mass_rho(max_cpd,my_i-1))/dt;
    rhs = rhs_fun(t,rho,xmid,R,dx,LENGTH,c_min,c_max,max_cpd,r_growth);
    for i=1:1:length(rho)
        rho(i)=rho(i)+dt*rhs(i);
        if(rho(i)>0.)
            spread(max_cpd,my_i)=xmid(i);
            time(my_i)=t;
        end
    end
    if((mod(my_i,im_write)>(im_write-2)))
        k=k+1;
            for i=1:1:length(rho)
            t_rho(i,k)=rho(i);
            end
    end
end
cells_per_day(max_cpd,my_i)=cells_per_day(max_cpd,my_i-1);    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ode_rhs = rhs_fun(t,rho,xmid,R,dx,LENGTH,c_min,c_max,max_cpd,r_growth)%arguments are position, refinement level and domain size
ode_rhs=zeros(size(rho));
rhs=zeros(size(rho));
N=length(rho);
for i=1:1:N
    rhs(i)=0.;%This stores mass in ball
    for j=1:1:N
    rhs(i)=rhs(i)+rho(j)*Indicator(abs(xmid(i)-xmid(j)),R)*dx;
    if((LENGTH/2-xmid(i)<R))%We are near right boundary (add on contribution from left boundary)
            rhs(i)=rhs(i)+rho(j)*Indicator(abs((xmid(i)-LENGTH)-xmid(j)),R)*dx;
    end
    if((-LENGTH/2-xmid(i)>-R))%We are near left boundary (add on contribution from right boundary)
       rhs(i)=rhs(i)+rho(j)*Indicator(abs((xmid(i)+LENGTH)-xmid(j)),R)*dx;
    end
    end
       ode_rhs(i)=r_growth*max_cpd*Indicator(c_min,rhs(i))*Indicator(rhs(i),c_max);
end
end
for i=1:10
    for j=1:1:2001
        NLR_fit=NLR_fit+(mass_rho(i,j+1)-meanNumMel{i}(j))^2;
    end
end
end
