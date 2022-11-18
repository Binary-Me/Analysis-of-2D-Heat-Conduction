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
tolerance=0.0001;

T(:,:,1)=T_avg; %the 1st guess value of temperature.
T(:,1,:)=32;
while(error>tolerance)
    i=i+1;
    for m=1:M
        for n=1:N
            k=conductivity(data,T(m,n,i));
            alpha=diffusivity(data,T(m,n,i));
            q_gen=-550*m*n*l*l+log(32*(n*l)^(0.38)+50*m*l);
            X=q_gen*l^2/k;
            Y=h*l*T_amb/k;
            if(n==1) %Bottom Wall
                T(m,1,i+1)=32;
                continue;
            elseif(m==1) 
                if(n==N) %Node at Top-Left Corner
                    T(1,N,i+1)=(1/beta)*(2*T(2,N,i)+2*T(1,N-1,i)+X+2*Y);  
                    continue;
                else %Left Wall - Exposed to Ambient Convection
                    T(m,n,i+1)=(1/beta)*(2*T(2,n,i)+T(1,n-1,i)+T(1,n+1,i)+X+2*Y);
                    continue;
                end
            elseif(m==M)
                q0=3.26*(n*l)^2;
                Z=2*q0*l/k;
                if(n==N) %Node at Top-Right Corner
                    T(m,n,i+1)=(1/beta)*(2*T(m-1,n,i)+2*T(m,n-1,i)+Z+X);
                    continue;
                else %Right Boundary - Constant heat flux into the wall
                    T(m,n,i+1)=0.25*(T(m,n+1,i)+T(m,n-1,i)+2*T(m-1,n,i)+Z+X);

                end
            elseif(n==N) %Top Boundary - Perfectly Insulated Wall
                T(m,n,i+1)=0.25*(T(m+1,n,i)+T(m-1,n,i)+2*T(m,n-1,i)+X);
                continue;
            else %Internal Nodes
                T(m,n,i+1)=0.25*(T(m+1,n,i)+T(m-1,n,i)+T(m,n+1,i)+T(m,n-1,i)+X);
                continue;
            end
        end
    end
    if(Error(T(:,:,i),T(:,:,i+1),tolerance,M,N)==1)
        disp("Values Converged")
        break;
    end
end
T(:,:,i+1:75000)=[];
T_min=min(T,[],'all');
T_max=max(T,[],'all');
T_plot=T(:,:,i)';
disp(T_plot(round(N/2),round(M/2)));

%% Plotting Graphs
x = linspace(0,L,M);
y = linspace(0,H,N);
 [X,Y]=meshgrid(x,y);
 contourf(X,Y,T_plot,50);
 cb1=colorbar;
 title(sprintf('Steady State plot'));
 caxis([T_min T_max]);
 xlabel('Length (in m)');
 ylabel('Height (in m)');
 
 figure("name","Temperature Variations along Height and Width of the plate in Steady State");
 t3=tiledlayout(1,2);
title(t3,"Temperature Profiles at Steady State");
nexttile
plot(T_plot(:,round(M/4)),y,'b',T_plot(:,round(M/2)),y,'m');
ylabel("Slab height y (in m)");
xlabel("T (in °C)");
legend('x=L/4','x=L/2');
nexttile
plot(x,T_plot(round(N/4),:),'k',x,T_plot(round(N/2),:),'r');
ylabel("Slab Length x (in m)");
xlabel("T (in °C)");
legend('Y=H/4','Y=H/2');