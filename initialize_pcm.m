function [T]=initialize_pcm(m,n,c,v1,data,mc,u)

d=zeros(m,n,c);
r=v1;
t=0;    
while t<100
    v=v1;
   for k=1:c
      for x=1:m
          for y=1:n
              d(x,y,k)=(data(x,y)-v(k)).^2;
          end
      end
   end
  
   
  for k=1:c
       for x=1:m
           for y=1:n
               u(x,y,k)=1./(1+(d(x,y,k)./(r(k)+0.00001)).^(1/(mc-1)));
           end
       end
  end 
  
   for k=1:c
        tp1=0.0;
        tp2=0.0;
        for x=1:m
            for y=1:n
                tp1=tp1+u(x,y,k)^mc*data(x,y);
                tp2=tp2+u(x,y,k)^mc;
            end
        end
        v1(k)=tp1/(tp2+0.000001);%中心进行迭代 找到最优中心
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
T=u;
end
