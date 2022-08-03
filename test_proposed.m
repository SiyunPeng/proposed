clc;
clear 
close all
%% 
data=imread('RI3.png');
data=double(data);
c=3;
data=padarray(data,[2,2],'replicate');
[m,n]=size(data);
%% 初始化参数
v=zeros(c,1);
mc=2;
p=2;
if c==2
    v1(1)=0;
    v1(2)=250;
elseif c==3
    v1(1)=0;
    v1(2)=125;
    v1(3)=250;
elseif c==4
    v1(1)=0;
    v1(2)=100;
    v1(3)=200;
    v1(4)=250;
elseif c==5
    v1(1)=0;
    v1(2)=70;
    v1(3)=150;
    v1(4)=200;
    v1(5)=250;
end

s_up=zeros(m,n,c);
s_down=zeros(m,n,c);
u1_u=zeros(m,n,c);
u1_d=zeros(m,n,c);
u2_u=zeros(m,n,c);
u2_d=zeros(m,n,c);
U=zeros(m,n,c);
G1=zeros(m,n,c);
G2=zeros(m,n,c);
K=zeros(m,n,c);

%% 计算核距离带宽sigma
%算xbar
tp1=0.0;
for x=3:m-2
    for y=3:n-2
        tp1=tp1+data(x,y);
    end
end
xbar=tp1/((m-4)*(n-4));

%算dbar
for x=3:m-2
    for y=3:n-2
        da(x,y)=abs(data(x,y)-xbar)^2;
    end
end
tp2=0.0;
for x=3:m-2
    for y=3:n-2
        tp2=tp2+da(x,y);
    end
end
dbar=tp2/(m*n);
%算sigma
tp3=0.0;
for x=3:m-2
    for y=3:n-2
        tp3=tp3+(da(x,y)-dbar)^2;
    end
end
sigma=sqrt(tp3/((m-4)*(n-4)-1));

%% 计算权衡加权模糊因子Wij
C=rand(m,n);
Cbar=rand(m,n);
eta=rand(m,n);
%计算Cj
for x=3:m-2
    for y=3:n-2
        tep1=0.0;
        daa=zeros(1,9);
        r=1;
        for x1=-1:1
            for y1=-1:1
                tep1=tep1+data(x+x1,y+y1);
                daa(1,r)=data(x+x1,y+y1);
                r=r+1;
            end
        end
        ave=tep1/9;
        C(x,y)=var(daa)/(ave^2+0.0001);
    end
end
%计算Cbar
for x=3:m-2
    for y=3:n-2
        tep2=0.0;
        for x1=-1:1
            for y1=-1:1
                tep2=tep2+C(x+x1,y+y1);
            end
        end
        Cbar(x,y)=tep2/9;
    end
end
%计算eta
for x=3:m-2
    for y=3:n-2
        tep3=0.0;
        for x1=-1:1
            for y1=-1:1
                tep3=tep3+exp(-(C(x+x1,y+y1)-Cbar(x+x1,y+y1)));
            end
        end
        eta(x,y)=exp(-(C(x,y)-Cbar(x,y)))/tep3;
    end
end
%计算Wgc
for x=3:m-2
    for y=3:n-2
        for x1=-1:1
            for y1=-1:1
                if C(x+x1,y+y1)<Cbar(x+x1,y+1)
                    Wgc(x+x1,y+y1)=2+eta(x+x1,y+y1);
                else
                    Wgc(x+x1,y+y1)=2-eta(x+x1,y+y1);
                end
            end
        end
    end
end

ssim=rand(m,n,9);
for x=3:m-2
    for y=3:n-2
        tep1=0.0;
        daa1=zeros(1,9);
        r2=1;
        for x1=-1:1
            for y1=-1:1
                tep1=tep1+data(x+x1,y+y1);
                daa1(1,r2)=data(x+x1,y+y1);
                r2=r2+1;
            end
        end
        miu1=tep1/9;
        r=1;
        for x1=-1:1
            for y1=-1:1
                tep2=0.0;
                daa2=zeros(1,9);
                r1=1;
                for x2=-1:1
                    for y2=-1:1
                        tep2=tep2+data(x+x1+x2,y+y1+y2);
                        daa2(1,r1)=data(x+x1+x2,y+y1+y2);
                        r1=r1+1;
                    end
                end
                miu2=tep2/9;
                tepp=0.0;
                for i=1:9
                    tepp=tepp+(daa1(i)-miu1)*(daa2(i)-miu2);
                end
                tepp=tepp/8;
                ssim(x,y,r)=((2*miu1*miu2+(0.01*255)^2)*(2*tepp+(0.03*255)^2))/((miu1^2+miu2^2+(0.01*255)^2)*(var(daa1)+var(daa2)+(0.03*255)^2));
                r=r+1;
            end
        end
    end
end

%% 求Wij
for x=3:m-2
    for y=3:n-2
        r=1;
        for x1=-1:1
            for y1=-1:1
                Wij(x+x1,y+y1)=(1/(sqrt(x1^2+y1^2)+1+(0.2*(1-ssim(x,y,r)))))*Wgc(x+x1,y+y1);
                r=r+1;
            end
        end
    end
end


Xq_new=zeros(m,n,3,3,9);
for x=3:m-2
    for  y=3:n-2
        count=0;
        Xq=zeros(3,3,9);
        for i1=-1:1
            for j1=-1:1
                count=count+1;
                for i2=-1:1
                    for j2=-1:1
                        %转换
                        if i2==-1
                            xx=1;
                        end
                        if j2==-1
                            yy=1;
                        end
                        if i2==0
                            xx=2;
                        end
                        if j2==0
                            yy=2;
                        end
                        if i2==1
                            xx=3;
                        end
                        if  j2==1
                            yy=3;
                        end
                        Xq(xx,yy,count)=data(x+i1+i2,y+j1+j2);
                    end
                end
                
            end
        end
        Xq_new(x,y,:,:,:)=Xq;
    end
end

%% 开始循环计时
[u1,v1]=initialize_fcm(m,n,c,v1,data,mc);
[t]=initialize_fcm(m,n,c,v1,data,mc);
e=0.0001;
ct=0;
err=0.1;
T1=t;
T2=t;
u_down=u1;
u_up=u1;
yt_up=ones(c,1);
yt_down=ones(c,1);
tt=clock;
d=zeros(m,n,c);
Vn=zeros(3,3,c);
while err>e && ct<200
    v=v1;
    for k=1:c
        temp3=zeros(3,3);
        temp4=zeros(3,3);
        for x=3:m-2
            for y=3:n-2
                temp3=temp3+ (u1(x,y,k).^mc).*(  reshape(Xq_new(x,y,:,:,5),3,3) ) ;
                temp4=temp4+ (u1(x,y,k).^mc);
            end
        end
        Vn(:,:,k)=temp3./temp4;
    end
    
    for  x=3:m-2
        for y=3:n-2
            XX=reshape(Xq_new(x,y,:,:,5),3,3);
            for k=1:c
                count2=0;
                for i1=1:3
                    for j1=1:3
                        if i1~=2 && j1~=2
                            count2 = count2 + ( XX(i1,j1) - Vn(i1,j1,k) ).^2;
                        end
                    end
                end
                d(x,y,k)=count2./9+0.25*((data(x,y)-v(k))^2);
            end
        end
    end
    
    %%
    %计算核距离K
    for k=1:c
        for x=3:m-2
            for y=3:n-2
                K(x,y,k)=exp(-(d(x,y,k))/sigma);
            end
        end
    end

    %%
    %计算邻域G
    for k=1:c
        for x=3:m-2
            for y=3:n-2
                tp1=0.0;
                tp2=0.0;
                for x1=-1:1
                    for y1=-1:1
                        if x1~=0 || y1~=0
                            tp1=tp1+Wij(x+1,y+y1)*((1-u_up(x+x1,y+y1,k))^mc*(1-T1(x+x1,y+y1,k))^p)*(1-K(x+x1,y+y1,k));
                            tp2=tp2+Wij(x+x1,y+y1)*((1-u_down(x+x1,y+y1,k))^mc*(1-T2(x+x1,y+y1,k))^p)*(1-K(x+x1,y+y1,k));
                        end
                    end
                end
                G1(x,y,k)=tp1;
                G2(x,y,k)=tp2;
            end
        end
    end
    
    
    %计算Sik
    for x=3:m-2
        for y=3:n-2
            for k=1:c
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
                                        st1=st1+sqrt((x2)^2+(y2)^2)*((1-K(x+x1+x2,y+y1+y2,k))^0.5);
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
    
    
  %%  
    for x=3:m-2
        for y=3:n-2
            for k=1:c
                d1(x,y,k)=T1(x,y,k)*(1-K(x,y,k))+G1(x,y,k)+0.0001;
                d2(x,y,k)=T2(x,y,k)*(1-K(x,y,k))+G2(x,y,k)+0.0001;
                d3(x,y,k)=T1(x,y,k)*s_up(x,y,k)^2+G1(x,y,k)+0.0001;
                d4(x,y,k)=T2(x,y,k)*s_down(x,y,k)^2+G2(x,y,k)+0.0001;
            end
            tep1=0.0;
            tep2=0.0;
            for k=1:c
                tep1=tep1+d1(x,y,k)^(-1/(mc-1));
                tep2=tep2+d2(x,y,k)^(-1/(mc-1));
            end
            for k=1:c
                u1_u(x,y,k)=d1(x,y,k)^(-1/(mc-1))/tep1;
                u1_d(x,y,k)=d2(x,y,k)^(-1/(mc-1))/tep2;
                u2_u(x,y,k)=d3(x,y,k)^(-1/(mc-1))/tep1;
                u2_d(x,y,k)=d4(x,y,k)^(-1/(mc-1))/tep2;
            end
        end
    end
    %% 
    %计算T(x,y,k)
    for k=1:c
        tp1=0;
        tp2=0;
        tp3=0;
        tp4=0;
        for x=3:m-2
            for y=3:n-2
                tp1=tp1+T1(x,y,k)*u_up(x,y,k)^mc*(1-K(x,y,k));
                tp2=tp2+T1(x,y,k)*u_up(x,y,k)^mc;
                tp3=tp3+T2(x,y,k)*u_down(x,y,k)^mc*(1-K(x,y,k));
                tp4=tp4+T2(x,y,k)*u_down(x,y,k)^mc;               
            end
        end
        yt_up(k)=tp1/(tp2+0.0001);
        yt_down(k)=tp3/(tp4+0.0001);
    end
    
    for x=3:m-2
        for y=3:n-2
            for k=1:c
                T1(x,y,k)=1/((1+2*((1-K(x,y,k)+G1(x,y,k))/(yt_up(k)+0.0001)))^(1/(p-1)));
                T2(x,y,k)=1/((1+2*((1-K(x,y,k)+G2(x,y,k))/(yt_down(k)+0.0001)))^(1/(p-1)));
            end
        end
    end
    
    %%
    for x=3:m-2
        for y=3:n-2
            for k=1:c
                u_up(x,y,k)=max(u1_u(x,y,k),u2_u(x,y,k));
                u_down(x,y,k)=min(u1_d(x,y,k),u2_d(x,y,k));
                T1(x,y,k)=max(T1(x,y,k),T2(x,y,k));
                T2(x,y,k)=min(T1(x,y,k),T2(x,y,k));
            end
        end
    end
    
%% 
    for k=1:c
        tepp1=0;
        tepp2=0;
        for x=3:m-2
            for y=3:n-2
                u(x,y,k)=(u_down(x,y,k)+u_up(x,y,k))/2;
                T(x,y,k)=(T1(x,y,k)+T2(x,y,k))/2;
                tepp1=tepp1+(u(x,y,k)^mc*T(x,y,k)^p)*data(x,y)*K(x,y,k);
                tepp2=tepp2+(u(x,y,k)^mc*T(x,y,k)^p)*K(x,y,k);
            end
        end
        v1(k)=tepp1/(tepp2+0.0001);
    end
    
    %降维
    xx=zeros(c,2);
    yy=zeros(c,2);
    u_new_R=zeros(m,n,c);
    u_new_L=zeros(m,n,c);
    T1_new=zeros(m,n,c);
    T2_new=zeros(m,n,c);
    v_R=zeros(c,1);
    v_L=zeros(c,1);
    v11=v1;
    v22=v1;
    
    data1=sort(data(:));
    data1=reshape(data1,m,n);
    comparision=0;
    while (comparision==0)
        for k=1:c
            for y=3:n-2
                for x=3:m-3
                    if v11(k)>data1(x,y) && v11(k)<data1(x+1,y)
                        xx(k,1)=x;
                        xx(k,2)=y;
                    end
                end
            end
        end
        
        for k=1:c
            ctt=0;
            for x=3:m-2
                for y=3:n-2
                    b=(xx(k,2)-1)*m+xx(k,1);
                    if ctt<=b
                        u_new_R(x,y,k)=u_down(x,y,k);
                        T1_new(x,y,k)=T2(x,y,k);
                    else
                        u_new_R(x,y,k)=u_up(x,y,k);
                        T1_new(x,y,k)=T1(x,y,k);
                    end
                    ctt=ctt+1;
                end
            end
        end
        
        
        for k=1:c
            tpe1=0;
            tpe2=0;
            for x=3:m-2
                for y=3:n-2
                    tpe1=tpe1+(u_new_R(x,y,k)*T1_new(x,y,k))*data(x,y)*K(x,y,k);
                    tpe2=tpe2+(u_new_R(x,y,k)*T1_new(x,y,k))*K(x,y,k);
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
            for y=3:n-2
                for x=3:m-3
                    if v22(k)>data1(x,y) && v22(k)<data1(x+1,y)
                        yy(k,1)=x;
                        yy(k,2)=y;
                    end
                end
            end
        end
        
        
        for k=1:c
            ctt1=0;
            for x=3:m-2
                for y=3:n-2
                    b1=(xx(k,2)-1)*m+xx(k,1);
                    if ctt1<=b1
                       u_new_L(x,y,k)=u_up(x,y,k);
                        T2_new(x,y,k)=T1(x,y,k);
                    else
                        u_new_L(x,y,k)=u_down(x,y,k);
                        T2_new(x,y,k)=T2(x,y,k);
                    end
                    ctt1=ctt1+1;
                end
            end
        end
        
        
        for k=1:c
            tpe1=0;
            tpe2=0;
            for x=3:m-2
                for y=3:n-2
                    tpe1=tpe1+(u_new_L(x,y,k)*T2_new(x,y,k))*data(x,y)*K(x,y,k);
                    tpe2=tpe2+(u_new_L(x,y,k)*T2_new(x,y,k))*K(x,y,k);
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
        for x=3:m-2
            for y=3:n-2
                U(x,y,k)=(u_new_R(x,y,k)+u_new_L(x,y,k))/2;
            end
        end
        v1(k)=(v_L(k)+v_R(k))/2;
    end
    err=0;
    for k=1:c
        err=err+(v1(k)-v(k))^2;
    end
    fprintf('迭代次数: %d ; Du: %f\n',ct,err);
    ct=ct+1;
    u1=U;
end
disp(['运行时间:  ',num2str(etime(clock,tt))]);

%%输出图片
image2=zeros(m,n);
for x=3:m-2
    for y=3:n-2
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
image2=image2(3:m-2,3:n-2);
figure(1);imshow(uint8(image2));

