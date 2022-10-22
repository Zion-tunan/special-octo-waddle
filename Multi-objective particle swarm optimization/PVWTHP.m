%% 风力出力、光伏出力、蒸发冷凝温度、热泵出力等构建
%% 程序编写：谢志远

%% 构建函数
function [Ppv,Pwt,Qhpsh,Qhpsc,Php,Tk,Tgeout] = PVWTHP(~)
%% PV
load solar.txt
load Temperature.txt
Pstc=0.2;                                    %额定容量、输出功率
Gac=solar;                                   %实时光照强度kw/m2  【8760,1】的表格
Gstc=1;                                      %标准测试条件下的光照强度，取值1KW/m2
alpha=-0.0037;                               %功率温度系数-0.37%
Tstc=25;                                     %标准测试条件下电池温度，取值25摄氏度
NOCT=44;                                     %额定操作条件下的光伏电池温度

for i=1:8760
    Tc(i,1)=Temperature(i,1)+0.0325*Gac(i,1);
end

for i=1:8760
    Ppv(i,1)=Pstc*0.95*(Gac(i,1)/Gstc)*(1+alpha*(Tc(i,1)-Tstc));
end

Ppv(Ppv<0) = 0 ;
% Ppv(Ppv>0.2)=Pstc;

%% WT
% load SpeedOfWind.txt
% SpeedOfWind=SpeedOfWind+3;
load SpeedWind.txt
% SpeedOfWind=SpeedWind+1;
% SpeedOfWind(4000:7000,:)=SpeedOfWind(4000:7000,:)+1;
h=10;                                        %塔架高度
h1=10;                                       %所测得10米风速
vi=2.5;                                      %风力机的启动风速
vr=12;                                       %额定风速
vc=25;                                       %截止风速
Pr=30;                                       %风力机组额定功率30KW
for i=1:8760
    v(i,1)=SpeedWind(i,1)*(h/h1)^0.14;
    if v(i,1)>vc || v(i,1)<vi
        Pwt(i,1)=0;
    elseif v(i,1)>=vi && v(i,1)<=vr
        Pwt(i,1)=((v(i,1)-vi)/(vr-vi))*Pr;
    else
        Pwt(i,1)=Pr;
    end
end

%% 构建蒸发温度与冷凝温度
load CoolingFuhe.txt
load HeatingFuhe.txt

Tgeout=zeros(8760,1);
Te=zeros(8760,1);
Tk=zeros(8760,1);
x=zeros(8760,1);
a=zeros(8760,1);
b=zeros(8760,1);
c=zeros(8760,1);
d=zeros(8760,1);

for i=1:8760
    Tgeout(i,1)=16.5+rand(1,1)-0.5;
    Te(i,1)=-7.818852817914556e-7*(i-4380)^2+5;
    Tk(i,1)=-4.170054836221096e-7*(i-4380)^2+46.5;
    x(i,1)=i;
    
    a(i,1)=i;
    b(i,1)=a(i,1)/24;
    c(i,1)=fix(b(i,1));
    d(i,1)=b(i,1)-c(i,1);
end
    %供暖季part1
    %for i = 1:2520
    %    Tk(i,1)=45+unifrnd(-1,1);
    %    Te(i,1)=Tgeout(i,1);
    %end
    
    %制冷季
    %for i = 3625:5832
    %    Tk(i,1)=Tgeout(i,1)+10;
    %    Te(i,1)=6.5+rand(1,1);
    %end
    
    %供暖季part2
    %for i = 7632:8760
    %    Tk(i,1)=45+unifrnd(-1,1);
    %    Te(i,1)=Tgeout(i,1);
    %end
    
    %供暖季part1
    for i = 1:2520
        if d(i)<0.31
            Tk(i,1)=52.5+rand(1,1);
            Te(i,1)=Tgeout(i,1)-10;
        else
            Tk(i,1)=44.5+rand(1,1);
            Te(i,1)=Tgeout(i,1)-10;
        end
    end
    
    %制冷季
    for i = 3625:5832
        if d(i)<0.31
            Tk(i,1)=Tgeout(i,1)+10;
            Te(i,1)=1.5+rand(1,1);
        else
            Tk(i,1)=Tgeout(i,1)+10;
            Te(i,1)=9.5+rand(1,1);
        end
    end
    
    %供暖季part2
    for i = 7632:8760
        if d(i)<0.31
            Tk(i,1)=52.5+rand(1,1);
            Te(i,1)=Tgeout(i,1)-10;
        else
            Tk(i,1)=44.5+rand(1,1);
            Te(i,1)=Tgeout(i,1)-10;
        end
    end
 
%% 构建热泵出力

 Qhpsh=zeros(8760,1);
 Qhpsc=zeros(8760,1);
 Php=zeros(8760,1);
 COPnom=zeros(8760,1);

for i=1:8760
    
    Qhpsh(i,1)=0;
    Qhpsc(i,1)=0;
    Php(i,1)=0;
    COPnom(i,1)=0;
    
end

    %供暖季part1
    for i = 1:2520
         Qhpsh(i,1)=478.146995+16.072630*Te(i,1)-3.769829*Tk(i,1)+0.263906*Te(i,1)^2-0.019371*Tk(i,1)*Te(i,1)-0.004210*Tk(i,1)^2+0.001502*Te(i,1)^3-0.001204*Tk(i,1)*Te(i,1)^2-0.001009*Tk(i,1)^2*Te(i,1)-0.000155*Tk(i,1)^3;
         Qhpsc(i,1)=399.7203884+13.43636568*Te(i,1)-3.151494208*Tk(i,1)+0.220619915*Te(i,1)^2-0.016194115*Tk(i,1)*Te(i,1)+0.003519194*Tk(i,1)^2+0.001255542*Te(i,1)^3-0.001006578*Tk(i,1)*Te(i,1)^2-0.000843412*Tk(i,1)^2*Te(i,1)-0.000129538*Tk(i,1)^3;
         Php(i,1)=65.895380+2.458958*Te(i,1)-1.623147*Tk(i,1)+0.036755*Te(i,1)^2-0.093511*Tk(i,1)*Te(i,1)+0.047725*Tk(i,1)^2+0.000293*Te(i,1)^3-0.000764*Tk(i,1)*Te(i,1)^2+0.001074*Tk(i,1)^2*Te(i,1)-0.000229*Tk(i,1)^3;
         COPnom(i,1)=Qhpsh(i,1)/(Php(i,1)+0.0001);%HCOP=4.48,CCOP=5.14
    end
    
    %制冷季
    for i = 3625:5832
         Qhpsh(i,1)=478.146995+16.072630*Te(i,1)-3.769829*Tk(i,1)+0.263906*Te(i,1)^2-0.019371*Tk(i,1)*Te(i,1)-0.004210*Tk(i,1)^2+0.001502*Te(i,1)^3-0.001204*Tk(i,1)*Te(i,1)^2-0.001009*Tk(i,1)^2*Te(i,1)-0.000155*Tk(i,1)^3;
         Qhpsc(i,1)=399.7203884+13.43636568*Te(i,1)-3.151494208*Tk(i,1)+0.220619915*Te(i,1)^2-0.016194115*Tk(i,1)*Te(i,1)+0.003519194*Tk(i,1)^2+0.001255542*Te(i,1)^3-0.001006578*Tk(i,1)*Te(i,1)^2-0.000843412*Tk(i,1)^2*Te(i,1)-0.000129538*Tk(i,1)^3;
         Php(i,1)=76.53493799+2.855984575*Te(i,1)-1.88522225*Tk(i,1)+0.042689932*Te(i,1)^2-0.108609143*Tk(i,1)*Te(i,1)+0.05543062*Tk(i,1)^2+0.000340357*Te(i,1)^3-0.000887126*Tk(i,1)*Te(i,1)^2+0.001247734*Tk(i,1)^2*Te(i,1)-0.000265711*Tk(i,1)^3;
         COPnom(i,1)=Qhpsc(i,1)/(Php(i,1)+0.0001);%HCOP=4.48,CCOP=5.14
    end
    
    %供暖季part2
    for i = 7632:8760
         Qhpsh(i,1)=478.146995+16.072630*Te(i,1)-3.769829*Tk(i,1)+0.263906*Te(i,1)^2-0.019371*Tk(i,1)*Te(i,1)-0.004210*Tk(i,1)^2+0.001502*Te(i,1)^3-0.001204*Tk(i,1)*Te(i,1)^2-0.001009*Tk(i,1)^2*Te(i,1)-0.000155*Tk(i,1)^3;
         Qhpsc(i,1)=399.7203884+13.43636568*Te(i,1)-3.151494208*Tk(i,1)+0.220619915*Te(i,1)^2-0.016194115*Tk(i,1)*Te(i,1)+0.003519194*Tk(i,1)^2+0.001255542*Te(i,1)^3-0.001006578*Tk(i,1)*Te(i,1)^2-0.000843412*Tk(i,1)^2*Te(i,1)-0.000129538*Tk(i,1)^3;
         Php(i,1)=65.895380+2.458958*Te(i,1)-1.623147*Tk(i,1)+0.036755*Te(i,1)^2-0.093511*Tk(i,1)*Te(i,1)+0.047725*Tk(i,1)^2+0.000293*Te(i,1)^3-0.000764*Tk(i,1)*Te(i,1)^2+0.001074*Tk(i,1)^2*Te(i,1)-0.000229*Tk(i,1)^3;
         COPnom(i,1)=Qhpsh(i,1)/(Php(i,1)+0.0001);%HCOP=4.48,CCOP=5.14
    end

