%% 生成Te、Tk冬夏值，独立于其他程序。
%% 程序编写：谢志远

%% 画图程序
Tgeout=zeros(8760,1);
for i=1:8760
    Tgeout(i,1)=16.5+rand(1,1)-0.5;
    Te(i,1)=-7.818852817914556e-7*(i-4380)^2+5;
    Tk(i,1)=-4.170054836221096e-7*(i-4380)^2+46.5;
    x(i,1)=i;
    
        %供暖季part1
    for i = 1:2520
        Tk(i,1)=55;
        Te(i,1)=Tgeout(i,1)-10;
    end
    
    %制冷季
    for i = 3625:5832
        Tk(i,1)=Tgeout(i,1)+10;
        Te(i,1)=0;
    end
    
    %供暖季part2
    for i = 7632:8760
        Tk(i,1)=55;
        Te(i,1)=Tgeout(i,1)-10;
    end
end


disp('开始绘图...')
plot(x,Te(:,1),x,Tk(:,1));
title('蒸发温度与冷凝温度全年实际值');
% axis([0 1.5 0 14]);
xlabel('\fontsize{12}小时/h');
ylabel('\fontsize{12}温度/℃');