function alpha=diffusivity(data,T)
    T_min=0;
    T_max=0;
    alpha_min=0;
    alpha_max=0;
    for i =2:29
        if(T<data(i,2))
            T_min=data(i-1,2);
            T_max=data(i,2);
            alpha_min=data(i-1,4);
            alpha_max=data(i,4);
            break;
        end
    end
    alpha=alpha_min+((alpha_max-alpha_min)/(T_max-T_min))*(T-T_min);
    alpha=alpha*0.01; %Converting to m^2/sec
end
