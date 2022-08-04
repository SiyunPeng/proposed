clc;
clear
close all
data=imread('.png');
data=double(data);
figure(1);subplot(1,2,1);imshow(uint8(data));title('噪声图');

[m,n]=size(data);
%% 初始化参数
c=3;

v=zeros(c,1);
mc=1.5;
if c==2
    v1(1)=0;
    v1(2)=250;
end
if c==3
    v1(1)=0;
    v1(2)=125;
    v1(3)=250;
end
if c==4
    v1(1)=0;
    v1(2)=100;
    v1(3)=200;
    v1(4)=250;
end
d = zeros(m,n,c);
% 计算邻域信息
alpha=0.05;
gamma=0.499;
XK=zeros(m,n,3,3);
SS=zeros(m,n,3,3);
Xq_new=zeros(m,n,3,3,9);
for i=3:m-2
  for  j=3:n-2
     count=0; 
     Xq=zeros(3,3,9); 
     Dp1=zeros(1,9);
     Wp1=zeros(1,9);
     S_old1=zeros(1,9);
     S=zeros(3,3);
     
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
                    Xq(xx,yy,count)=data(i+i1+i2,j+j1+j2);
                end
            end

        end
    end
    Xq_new(i,j,:,:,:)=Xq;
    for kk=1:9     
      Dp1(kk)=round(sum(sum( abs( Xq(:,:,kk)-Xq(:,:,5) ) ) )./9);
%       Dp1(kk)=sum(sum( abs( Xq(:,:,kk)-Xq(:,:,5) ) ) )./9;
    end
    Dp=reshape(Dp1,3,3)';

    for k1=1:9
        FM=0;
        for ii=1:9
            FM=FM+exp(-alpha*Dp1(ii));
        end
        Wp1(k1)=exp(-alpha*Dp1(k1))/FM;
    end
    Wp=reshape(Wp1,3,3)';

    for k2=1:9
        S_old1(k2)=sum( sum( Wp.*abs( Xq(:,:,k2)-Xq(:,:,5)) ) )/9;
    end
    S_old=reshape(S_old1,3,3)';
    
    for i4=1:3
        for j4=1:3
            S(i4,j4)=exp(-gamma*S_old(i4,j4));
        end
    end
    SS(i,j,:,:)=S;
  end
end

[u0,v1]=fcm_11(m,n,c,v1,data,mc);

ee = 0.01;it = 1;
t1=clock;
V=zeros(c,1);
Vn=zeros(3,3,c);
u=u0;
while ee>0.001 && it<10     %终止条件
    
%   更新聚类中心          
     for k=1:c
        temp1=0;
        temp2=0;
        temp3=zeros(3,3);
        temp4=zeros(3,3);
        for i=3:m-2
            for j=3:n-2  
                M=reshape(Xq_new(i,j,:,:,5),3,3);
                N=reshape(SS(i,j,:,:),3,3);  
               temp1=temp1+ (u0(i,j,k).^mc).*( data(i,j) );
               temp2=temp2+( u0(i,j,k).^mc);    
               temp3=temp3+ (u0(i,j,k).^mc).*(  reshape(Xq_new(i,j,:,:,5),3,3) ) ;
               temp4=temp4+ (u0(i,j,k).^mc);  
            end
        end  
        V(k)=temp1/temp2;
        Vn(:,:,k)=temp3./temp4;
    end

% 算距离:样本data(i,j)到第k类的距离
   for  i=3:m-2
    for j=3:n-2
        XX=reshape(Xq_new(i,j,:,:,5),3,3);
        for k=1:c
            
            count2=0;
            for i1=1:3
                for j1=1:3
                    if i1~=2 && j1~=2
                          count2 = count2 + ( XX(i1,j1) - Vn(i1,j1,k) ).^2;
                    end
                end
            end
            d(i,j,k)=count2./9 +( (data(i,j)-V(k))^2 );    
        end 
    end
   end

%算Gij
G=zeros(m,n,c);
for k=1:c
    for i=3:m-2
        for j=3:n-2
            temp=0;
            for i1=-1:1
                for j1=-1:1
                    temp=temp+M(i1+2,j1+2).*( (1-u0(i+i1,j+j1,k)).^mc ).*d(i+i1,j+j1,k);
                end
            end
            G(i,j,k)=temp;
        end
    end
end

%  更新隶属度
   for k=1:c
       for i=3:m-2
           for j=3:n-2
               N=reshape(SS(i,j,:,:),3,3);  
               temp=0;
               for i1=-1:1
                   for j1=-1:1
                       temp=temp+ ( d(i,j)+G(i,j,k) );
                   end
               end
               u_up=(temp).^(-1/(mc-1));
               temp2=0;
               for k1=1:c
                   temp1=0;
                   for i2=-1:1
                       for j2=-1:1
                           temp1=temp1+( d(i,j) + G(i,j,k1));
                       end
                   end
                   temp2=temp2+ (temp1.^(-1/(mc-1)) );
               end
               u_down=temp2;
               u(i,j,k)=u_up./u_down;
           end
       end
   end
   
%  终止条件
   temp=max( max( max( abs(u0-u) ) ));
   if   temp<0.0001
        ee=0.00001;
   end
   fprintf('迭代次数: %d ; Du: %f\n',it,temp);
   it=it+1;
   u0=u;
end

disp(['运行时间：',num2str(etime(clock,t1))]);

% 聚类
image2=zeros(m,n);
for i=1:m
    for j=1:n
        if c==2
            if u(i,j,1)>u(i,j,2)
                image2(i,j)=0;
            else
                image2(i,j)=255;
            end
        end
        if c==3
            if u(i,j,1)>=u(i,j,2) && u(i,j,1)>=u(i,j,3)
                image2(i,j)=0;
            end
            if u(i,j,2)>=u(i,j,1) && u(i,j,2)>=u(i,j,3)
                image2(i,j)=125;
            end
            if u(i,j,3)>=u(i,j,1) && u(i,j,3)>=u(i,j,2)
                image2(i,j)=255;
            end
        end
        if c==4
            if u(i,j,1)>u(i,j,2) && u(i,j,1)>u(i,j,3) && u(i,j,1)>u(i,j,4)
                image2(i,j)=0;
            elseif u(i,j,2)>u(i,j,1) && u(i,j,2)>u(i,j,3) && u(i,j,2)>u(i,j,4)
                image2(i,j)=64;
            elseif u(i,j,3)>u(i,j,1) && u(i,j,3)>u(i,j,2) && u(i,j,3)>u(i,j,4)
                image2(i,j)=150;
            else 
                image2(i,j)=255;
            end
        end
    end
end
figure(2);imshow(uint8(image2));