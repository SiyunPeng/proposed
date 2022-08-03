function [u,v1]=initialize_fcm(m,n,c,v1,data,mc)



d=zeros(m,n,c);
u=zeros(m,n,c);


t=0;
while t<100
    v=v1;
    for k=1:c
        for i=1:m
            for j=1:n
                d(i,j,k)=(data(i,j)-v(k))^2+0.0001;
            end
        end
    end
    for i=1:m
        for j=1:n
            tp1=0.0;
            for k=1:c
                tp1=tp1+d(i,j,k)^(-1/(mc-1));
            end
            for k=1:c
                u(i,j,k)=d(i,j,k)^(-1/(mc-1))/tp1;
            end
        end
    end
    for k=1:c
        tp1=0.0;
        tp2=0.0;
        for i=1:m
            for j=1:n
                tp1=tp1+u(i,j,k)^mc*data(i,j);
                tp2=tp2+u(i,j,k)^mc;
            end
        end
        v1(k)=tp1/(tp2+0.0001);
    end
    temp=0.0;
    for k=1:c
        temp=temp+(v(k)-v1(k))^2;
    end
    if temp<0.0001
        break;
    end
    t=t+1;
end
end