function[t,image2]=FLICM(m,n,c,v1,data,mc)
G=zeros(m,n,c);
e=0.0001;
t=0;
err=0.1;
tt=clock;

[U,v1]=initialize_fcm(m,n,c,v1,data,mc);
while err>e && t<200
    v=v1;
    %计算邻域G
    for k=1:c
        for x=2:m-1
            for y=2:n-1
                tp1=0.0;
                for x1=-1:1
                    for y1=-1:1
                        if x1~=0 || y1~=0
                            tp1=tp1+((1-U(x+x1,y+y1,k))^mc)*((data(x+x1,y+y1)-v1(k))^2)/((x1^2+y1^2)^0.5+1);
                        end
                    end
                end
                G(x,y,k)=tp1;
            end
        end
    end
    %计算隶属度u
    for x=1:m
        for y=1:n
            tepp1=0.0;
            for k=1:c
                tepp1=tepp1+((data(x,y)-v(k))^2+G(x,y,k)+0.0001)^(-1/(mc-1));
            end
            for k=1:c
                U(x,y,k)=((data(x,y)-v(k))^2+G(x,y,k)+0.0001)^(-1/(mc-1))/tepp1;
            end
        end
    end
    %计算聚类中心v
    for k=1:c
        tp1=0.0;
        tp2=0.0;
        for x=1:m
            for y=1:n
                tp1=tp1+(U(x,y,k)^mc)*data(x,y);
                tp2=tp2+(U(x,y,k)^mc);
            end
        end
        v1(k)=tp1/tp2;
    end
    err=0;
    for k=1:c
        err=err+(v1(k)-v(k))^2;
    end
    t=t+1;
end
disp(['运行时间:  ',num2str(etime(clock,tt))]);
%%
image2=zeros(m,n);
for x=1:m
    for y=1:n
        if c==2
            if U(x,y,1)>U(x,y,2)
                image2(x,y)=0;
            else
                image2(x,y)=255;
            end
        end
        if c==3
            if U(x,y,1)>U(x,y,2) && U(x,y,1)>U(x,y,3)
                image2(x,y)=0;
            end
            if U(x,y,2)>U(x,y,1)  && U(x,y,2)>U(x,y,3)
                image2(x,y)=125;
            end
            if U(x,y,3)>U(x,y,1)&& U(x,y,3)>U(x,y,2)
                image2(x,y)=255;
            end
        end
        if c==4
            if U(x,y,1)>U(x,y,2) && U(x,y,1)>U(x,y,3)&&U(x,y,1)>U(x,y,4)
                image2(x,y)=0;
            end
            if U(x,y,2)>U(x,y,1) && U(x,y,2)>U(x,y,3)&&U(x,y,2)>U(x,y,4)
                image2(x,y)=64;
            end
            if U(x,y,3)>U(x,y,1)&& U(x,y,3)>U(x,y,2)&& U(x,y,3)>U(x,y,4)
                image2(x,y)=150;
            end
            if U(x,y,4)>U(x,y,1)&& U(x,y,4)>U(x,y,2)&& U(x,y,4)>U(x,y,3)
                image2(x,y)=255;
            end
        end
    end
end


