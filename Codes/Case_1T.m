% Stainless Steel Type 304L @ T=T_avg
data = table2array(readtable("Thermal_Properties.xlsx","Range","C6:F34"));
L=0.08;
H=0.005;
T_F=200;
T_C=(T_F-32)*(5/9);
T_iso=32;
T_avg=(T_C+T_iso)/2;
k=conductivity(data,T_avg); %Interpolated from the datatable
alpha=diffusivity(data,T_avg); %Interpolated from the datatable
h=10;
T_amb=-10;
l=0.0002;

beta=4+(2*h*l/k);
dt=l^2/(alpha*beta);
tau=alpha*dt/l^2;
M=round(L/l);
N=round(H/l);
i=0;
T=zeros(M,N,75000);
error=10;
tolerance=0.01;

T(:,:,1)=T_C;
T(:,1,:)=32;
while(error>tolerance)
    i=i+1;
    for m=1:M
        for n=1:N
            q_gen=-550*m*n*l*l+log(32*(n*l)^(0.38)+50*m*l);
            X=q_gen*l^2/k;
            Y=h*l*T_amb/k;
            if(n==1) %Bottom Wall
                T(m,1,i+1)=32;
                continue;
            elseif(m==1) 
                if(n==N) %Node at Top-Left Corner
                    T(1,N,i+1)=2*tau*(T(2,N,i)+T(1,N-1,i))+T(1,N,i)*(1-(beta*tau))+tau*(X+2*Y);  
                    continue;
                else %Left Wall - Exposed to Ambient Convection
                    T(m,n,i+1)=tau*(T(m,n+1,i)+T(m,n-1,i)+2*T(m+1,n,i))+(1-(beta*tau))*T(m,n,i)+((X+2*Y)*tau);
                    continue;
                end
            elseif(m==M)
                q0=3.26*(n*l)^2;
                Z=2*q0*l/k;
                if(n==N) %Node at Top-Right Corner
                    T(m,n,i+1)=2*tau*(T(m-1,n,i)+T(m,n-1,i))+(1-4*tau)*T(m,n,i)+(Z+X)*tau;
                    continue;
                else %Right Boundary - Constant heat flux into the wall
                    T(m,n,i+1)=tau*(T(m,n+1,i)+T(m,n-1,i)+2*T(m-1,n,i))+(1-4*tau)*T(m,n,i)+(Z+X)*tau;

                end
            elseif(n==N) %Top Boundary - Perfectly Insulated Wall
                T(m,n,i+1)=tau*(T(m+1,n,i)+T(m-1,n,i)+2*T(m,n-1,i))+(1-4*tau)*T(m,n,i)+X*tau;
                continue;
            else %Internal Nodes
                T(m,n,i+1)=tau*(T(m+1,n,i)+T(m-1,n,i)+T(m,n+1,i)+T(m,n-1,i))+(1-4*tau)*T(m,n,i)+X*tau;
                continue;
            end
        end
    end
    if(Error(T(:,:,i),T(:,:,i+1),tolerance,M,N)==1)
        break;
    end
end
T(:,:,i+1:75000)=[];
T_plot=zeros(N,M,i-1);
for temp=1:i
    T_plot(:,:,temp)=T(:,:,temp)';
end
total_time=dt*i;
T_min=min(T,[],'all');
T_max=max(T,[],'all');
fprintf("The plate reaches steady state at: %.5f seconds", total_time);

%% Plotting Graphs
%Temperature Contour Plot 
time_gap=30;
counter=1;
figure("name",'Heat Transfer Visualization');
x = linspace(0,L,M);
y = linspace(0,H,N);
t1=tiledlayout(1,1);
nexttile
for j=1:time_gap:i
    [X,Y] = meshgrid(x,y);
    Z = T_plot(:,:,j);
    contourf(X,Y,Z,50);
    cb=colorbar;
    title(sprintf('Temperatures at time : %.5f seconds ',(dt*time_gap*counter)))
    caxis([T_min T_max]);
    xlabel('Length (in m)');
    ylabel('Height (in m)');
    pause(2)
    counter=counter+1;
end

%Distribution plots at four different times with Equal Intervals
p=round(linspace(1,i,5));
figure('name','Temperature Distribution plots at 4 time instances');
t2=tiledlayout(2,2);
x = linspace(0,L,M);
y = linspace(0,H,N);
title(t2,"Temperature Distribution Plots");
for o=1:4
    nexttile
    [X,Y]=meshgrid(x,y);
    Z = T_plot(:,:,p(o+1));
    contourf(X,Y,Z,50);
    cb1=colorbar;
    title(sprintf('At time t= %.5f seconds ',(dt*p(o+1))));
    caxis([T_min T_max]);
    xlabel('Length (in m)');
    ylabel('Height (in m)');
end