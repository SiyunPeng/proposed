function[t,image2]=AIT2FCM(m,n,c,v1,data,mc)
s_up=zeros(m,n,c);
s_down=zeros(m,n,c);
u1=zeros(m,n,c);
u2_u=zeros(m,n,c);
u2_d=zeros(m,n,c);
u_up=zeros(m,n,c);
u_down=zeros(m,n,c);
[U,v1]=initialize_fcm(m,n,c,v1,data,mc);
e=0.0001;
t=0;
err=0.1;
tt=clock;
while err>e && t<200
    v=v1;
    for x=1:m
        for y=1:n
            for k=1:c
                dt(x,y,k)=(data(x,y)-v(k))^2+0.0001;
            end
            tp1=0.0;
            for k=1:c
                tp1=tp1+(dt(x,y,k))^(-1/(mc-1));
            end
            for k=1:c
                u1(x,y,k)=(dt(x,y,k))^(-1/(mc-1))/tp1;
            end
        end
    end
    
    for k=1:c
        for  x=3:m-2
            for y=3:n-2
                Dx=zeros(8,1);
                for x1=-1:1
                    for y1=-1:1
                        if x1~=0 || y1~=0
                            st1=0;
                            Dr=0;
                            for x2=-1:1
                                for y2=-1:1
                                    if x2~=0 || y2~=0
                                        Dr=Dr+sqrt((x2)^2+(y2)^2);
                                        st1=st1+sqrt((x2)^2+(y2)^2)*abs(data(x+x1+x2,y+y1+y2)-v(k));
                                    end
                                end
                            end
                            sx=3*(x1+1)+(y1+2);
                            if (sx>=6)
                                sx=sx-1;
                            end
                            Dx(sx)=st1/Dr;
                        end
                    end
                end
                s_up(x,y,k)=max(reshape(Dx,numel(Dx),1));
                s_down(x,y,k)=min(reshape(Dx,numel(Dx),1));
                
            end
        end
    end
    
    
    
    % 算隶属度
    for x=3:m-2
        for y=3:n-2
            for k=1:c
                dt(x,y,k)=(data(x,y)-v(k))^2+0.0001;
            end
            tp1=0.0;
            for k=1:c
                tp1=tp1+(dt(x,y,k))^(-1/(mc-1));
            end
            for k=1:c
                u2_u(x,y,k)=s_up(x,y,k)^(-2/(mc-1))/tp1;
                u2_d(x,y,k)=s_down(x,y,k)^(-2/(mc-1))/tp1;
            end
        end
    end
    
    
    for x=1:m
        for y=1:n
            for k=1:c
                u_up(x,y,k)=max(u1(x,y,k),u2_u(x,y,k));
                u_down(x,y,k)=min(u1(x,y,k),u2_d(x,y,k));
            end
        end
    end
    
    for x=1:m
        for y=1:n
            sum0=0;
            sum1=0;
            for k=1:c
                sum0=sum0+u_up(x,y,k);
                sum1=sum1+u_down(x,y,k);
            end
            for k=1:c
                u_up(x,y,k)=u_up(x,y,k)/(sum0+0.00001);
                u_down(x,y,k)=u_down(x,y,k)/(sum1+0.00001);
            end
        end
    end
    
    for k=1:c
        tepp1=0;
        tepp2=0;
        for x=1:m
            for y=1:n
                u(x,y,k)=(u_down(x,y,k)+u_up(x,y,k))/2;
                tepp1=tepp1+u(x,y,k)*data(x,y);
                tepp2=tepp2+u(x,y,k);
            end
        end
        v1(k)=tepp1/(tepp2+0.0001);
    end
    
    
    %%降维
    xx=zeros(c,2);
    yy=zeros(c,2);
    u_new_R=zeros(m,n,c);
    u_new_L=zeros(m,n,c);
    v_R=zeros(c,1);
    v_L=zeros(c,1);
    v11=v1;
    v22=v1;
    
    data1=sort(data(:));
    data1=reshape(data1,m,n);
    comparision=0;
    while (comparision==0)
        for k=1:c
            for y=1:n
                for x=1:m-1
                    if v11(k)>data1(x,y) && v11(k)<data1(x+1,y)
                        xx(k,1)=x;
                        xx(k,2)=y;
                    end
                end
            end
        end
        
        for k=1:c
            ctt=0;
            for x=1:m
                for y=1:n
                    b=(xx(k,2)-1)*m+xx(k,1);
                    if ctt<=b
                        u_new_R(x,y,k)=u_down(x,y,k);
                    else
                        u_new_R(x,y,k)=u_up(x,y,k);
                    end
                    ctt=ctt+1;
                end
            end
        end
        
        
        for k=1:c
            tpe1=0;
            tpe2=0;
            for x=1:m
                for y=1:n
                    tpe1=tpe1+u_new_R(x,y,k)*data(x,y);
                    tpe2=tpe2+u_new_R(x,y,k);
                end
            end
            v_R(k)=tpe1/(tpe2+0.0001);
        end
        for k=1:c
            if (floor(v_R(k)) == floor(v11(k)))
                comparision=1;
            else
                v11(k)=v_R(k);
            end
        end
    end
    
    comparision1=0;
    while (comparision1==0)
        for k=1:c
            for y=1:n
                for x=1:m-1
                    if v22(k)>data1(x,y) && v22(k)<data1(x+1,y)
                        yy(k,1)=x;
                        yy(k,2)=y;
                    end
                end
            end
        end
        
        
        for k=1:c
            ctt1=0;
            for x=1:m
                for y=1:n
                    b1=(yy(k,2)-1)*m+yy(k,1);
                    if ctt1<b1
                        u_new_L(x,y,k)=u_up(x,y,k);
                    else
                        u_new_L(x,y,k)=u_down(x,y,k);
                    end
                    ctt1=ctt1+1;
                end
            end
        end
        
        
        for k=1:c
            tpe1=0;
            tpe2=0;
            for x=1:m
                for y=1:n
                    tpe1=tpe1+u_new_L(x,y,k)*data(x,y);
                    tpe2=tpe2+u_new_L(x,y,k);
                end
            end
            v_L(k)=tpe1/(tpe2+0.0001);
        end
        
        for k=1:c
            if (floor(v_L(k)) == floor(v22(k)))
                comparision1=1;
            else
                v22(k)=v_L(k);
            end
        end
    end
    
    for k=1:c
        for x=1:m
            for y=1:n
                U(x,y,k)=(u_new_R(x,y,k)+u_new_L(x,y,k))/2;
            end
        end
        v1(k)=(v_L(k)+v_R(k))/2;
    end
    
    err=0;
    for k=1:c
        err=err+(v1(k)-v(k))^2;
    end
    t=t+1;
end
disp(['运行时间:  ',num2str(etime(clock,tt))]);

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
        if c==5
            if U(x,y,1)>U(x,y,2) && U(x,y,1)>U(x,y,3)&&U(x,y,1)>U(x,y,4)&&U(x,y,1)>U(x,y,5)
                image2(x,y)=0;
            end
            if U(x,y,2)>U(x,y,1) && U(x,y,2)>U(x,y,3)&&U(x,y,2)>U(x,y,4)&&U(x,y,2)>U(x,y,5)
                image2(x,y)=65;
            end
            if U(x,y,3)>U(x,y,1)&& U(x,y,3)>U(x,y,2)&& U(x,y,3)>U(x,y,4)&& U(x,y,3)>U(x,y,5)
                image2(x,y)=130;
            end
            if U(x,y,4)>U(x,y,1)&& U(x,y,4)>U(x,y,2)&& U(x,y,4)>U(x,y,3)&& U(x,y,4)>U(x,y,5)
                image2(x,y)=200;
            end
            if U(x,y,5)>U(x,y,1)&& U(x,y,5)>U(x,y,2)&& U(x,y,5)>U(x,y,3)&& U(x,y,5)>U(x,y,4)
                image2(x,y)=255;
            end
        end
    end
end


