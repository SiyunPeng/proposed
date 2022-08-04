function[t,image2]=KWFLICM(m,n,c,v1,data,mc)
G=zeros(m,n,c);
K=zeros(m,n,c);

%% 计算核距离带宽sigma
%算xbar
tp1=0.0;
for x=1:m
    for y=1:n
        tp1=tp1+data(x,y);
    end
end
xbar=tp1/(m*n);

%算dbar
for x=1:m
    for y=1:n
        da(x,y)=abs(data(x,y)-xbar)^2;
    end
end
tp2=0.0;
for x=1:m
    for y=1:n
        tp2=tp2+da(x,y);
    end
end
dbar=tp2/(m*n);
%算sigma
tp3=0.0;
for x=1:m
    for y=1:n
    tp3=tp3+(da(x,y)-dbar)^2;
    end
end
sigma=sqrt(tp3/(m*n-1));

C=rand(m,n);
Cbar=rand(m,n);
eta=rand(m,n);
%计算Cj
for x=2:m-1
    for y=2:n-1
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
for x=2:m-1
    for y=2:n-1
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
for x=2:m-1
    for y=2:n-1
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
for x=2:m-1
    for y=2:n-1
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

 %求Wij
 for x=2:m-1
     for y=2:n-1
         for x1=-1:1
             for y1=-1:1
                 Wij(x+x1,y+y1)=(1/(sqrt(x1^2+y1^2)+1))*Wgc(x+x1,y+y1);
             end
         end
     end
 end
 
%% 开始循环计时
[U,v1]=initialize_fcm(m,n,c,v1,data,mc);
e=0.0001;
t=0;
err=0.1;
tt=clock;
while err>e && t<200
    v=v1;
    %计算核距离K
    for k=1:c
        for x=1:m
            for y=1:n
                K(x,y,k)=exp(-(data(x,y)-v(k))^2/sigma);
            end
        end
    end
    
   %计算邻域G
   for k=1:c
       for x=2:m-1
           for y=2:n-1
               tp1=0.0;
               for x1=-1:1
                   for y1=-1:1
                       if x1~=0 || y1~=0
                           tp1=tp1+Wij(x+x1,y+y1)*((1-U(x+x1,y+y1,k))^mc)*(1-K(x+x1,y+y1,k));
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
               tepp1=tepp1+(1-K(x,y,k)+G(x,y,k)+0.0001)^(-1/(mc-1));
           end
           for k=1:c
               U(x,y,k)=(1-K(x,y,k)+G(x,y,k)+0.0001)^(-1/(mc-1))/tepp1;
           end
       end
   end
   %计算聚类中心v
   for k=1:c
       tp1=0.0;
       tp2=0.0;
       for x=1:m
           for y=1:n
               tp1=tp1+(U(x,y,k)^mc)*K(x,y,k)*data(x,y);
               tp2=tp2+(U(x,y,k)^mc)*K(x,y,k);
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
 