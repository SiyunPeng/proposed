clc;
clear all;
close all;
%% 
data=imread('RI3.png');
%%
figure(1);subplot(1,3,1);imshow(uint8(data));title('灰度图');
c=3;
data=double(data);
figure(1);subplot(1,3,3);imshow(uint8(data));title('噪声图')

%% 初始化参数
[m,n]=size(data);
v=zeros(c,1);
mc=2;
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

%%
[t2,image2]=FLICM(m,n,c,v1,data,mc);
figure(3); imshow(uint8(image2));title('FLICM');
fprintf('FLICM:%.4f\n',t2);

%%
[t6,image6]=AIT2FCM(m,n,c,v1,data,mc);
figure(7); imshow(uint8(image6));title('AIT2FCM');
fprintf('AIT2FCM:%.4f\n',t6);

%%
[t7,image7]=KWFLICM(m,n,c,v1,data,mc);
figure(8); imshow(uint8(image7));title('KWFLICM');
fprintf('KWFLICM:%.4f\n',t7);




