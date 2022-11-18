function k=conductivity(data,T)
    T_min=0;
    T_max=0;
    k_max=0;
    k_min=0;
    for i =2:29
        if(T<data(i,2))
            T_min=data(i-1,2);
            T_max=data(i,2);
            k_max=data(i,3);
            k_min=data(i-1,3);
            break;
        end
    end
    k=k_min+((k_max-k_min)/(T_max-T_min))*(T-T_min);
    k=k*100; %Converting to W/m-K
end
