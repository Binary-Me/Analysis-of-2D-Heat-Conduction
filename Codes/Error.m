function check=Error(T1,T2,tol,M,N)
    check=0;
    outer_check=0;
    for i=2:M
        for j=1:N
            if(abs(T1(i,j)-T2(i,j))>tol)
                check=0;
                outer_check=1;
                break;
            else
                check=1;
            end
        end
        if(outer_check==1)
            break;
        end
    end       
end

