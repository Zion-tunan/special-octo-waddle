%% 目标函数与约束条件
%% 程序编写：谢志远

%% 调用Test_function，输入变量位置信息
function PopObj = Fitness_fun(X,X_num,Test_function)

switch Test_function
    case 'ZDT1'
        PopObj(1,1) = X(1);
        g = 1+9*sum(X(2:end))/(X_num-1);
        PopObj(1,2) = g*(1-(PopObj(1)/g)^0.5);
    case 'ZDT2'
        PopObj(1,1) = X(1);
        g = 1+9*sum(X(2:end))/(X_num-1);
        PopObj(1,2) = g*(1-(PopObj(1)/g)^2);
    case 'ZDT3'
        PopObj(1,1) = X(1);
        g = 1+9*sum(X(2:end))/(X_num-1);
        temp = PopObj(1)/g;
        PopObj(1,2) = g*(1-(temp)^0.5-temp*sin(10*pi* PopObj(1,1)));
    case 'ZDT4'
        PopObj(1,1) = X(1);
        temp = X(2:end);
        g = 1+10*(X_num-1)+sum(temp.^2-10*cos(4*pi*temp));
        PopObj(1,2) = g*(1-(PopObj(1)/g)^0.5);
    case 'ZDT6'
        PopObj(1,1) = X(1);
        temp = X_num-1;
        g =1+9*(sum(X(2:end)/temp)^0.25);
        PopObj(1,2) = g*(1-(PopObj(1)/g)^2);
    case 'DTLZ1'
        M =3;
        [N,D]  = size(X);
        g      = 100*(D-M+1+sum((X(:,M:end)-0.5).^2-cos(20.*pi.*(X(:,M:end)-0.5)),2));
        PopObj = 0.5*repmat(1+g,1,M).*fliplr(cumprod([ones(N,1),X(:,1:M-1)],2)).*[ones(N,1),1-X(:,M-1:-1:1)];
    case 'DTLZ2'
        M      = 3;
        g      = sum((X(:,M:end)-0.5).^2,2);
        PopObj = repmat(1+g,1,M).*fliplr(cumprod([ones(size(g,1),1),cos(X(:,1:M-1)*pi/2)],2)).*[ones(size(g,1),1),sin(X(:,M-1:-1:1)*pi/2)];
    case 'DTLZ3'
        [N,D]  = size(X);
        M      = 3;
        g      = 100*(D-M+1+sum((X(:,M:end)-0.5).^2-cos(20.*pi.*(X(:,M:end)-0.5)),2));
        PopObj = repmat(1+g,1,M).*fliplr(cumprod([ones(N,1),cos(X(:,1:M-1)*pi/2)],2)).*[ones(N,1),sin(X(:,M-1:-1:1)*pi/2)];
    case 'DTLZ4'
        M      = 3;
        X(:,1:M-1) = X(:,1:M-1).^100;
        g      = sum((X(:,M:end)-0.5).^2,2);
        PopObj = repmat(1+g,1,M).*fliplr(cumprod([ones(size(g,1),1),cos(X(:,1:M-1)*pi/2)],2)).*[ones(size(g,1),1),sin(X(:,M-1:-1:1)*pi/2)];
    case 'DTLZ5'
        M      = 3;
        g      = sum((X(:,M:end)-0.5).^2,2);
        Temp   = repmat(g,1,M-2);
        X(:,2:M-1) = (1+2*Temp.*X(:,2:M-1))./(2+2*Temp);
        PopObj = repmat(1+g,1,M).*fliplr(cumprod([ones(size(g,1),1),cos(X(:,1:M-1)*pi/2)],2)).*[ones(size(g,1),1),sin(X(:,M-1:-1:1)*pi/2)];
    case 'DTLZ6'
        M      = 3;
        g      = sum(X(:,M:end).^0.1,2);
        Temp   = repmat(g,1,M-2);
        X(:,2:M-1) = (1+2*Temp.*X(:,2:M-1))./(2+2*Temp);
        PopObj = repmat(1+g,1,M).*fliplr(cumprod([ones(size(g,1),1),cos(X(:,1:M-1)*pi/2)],2)).*[ones(size(g,1),1),sin(X(:,M-1:-1:1)*pi/2)];
    case 'DTLZ7'
        M               = 3;
        PopObj          = zeros(size(X,1),M);
        g               = 1+9*mean(X(:,M:end),2);
        PopObj(:,1:M-1) = X(:,1:M-1);
        PopObj(:,M)     = (1+g).*(M-sum(PopObj(:,1:M-1)./(1+repmat(g,1,M-1)).*(1+sin(3*pi.*PopObj(:,1:M-1))),2));
    case 'UF1'
        D = size(X,2);
        J1 = 3:2:D;
        J2 = 2:2:D;
        Y  = X - sin(6*pi*repmat(X(:,1),1,D)+repmat(1:D,size(X,1),1)*pi/D);
        PopObj(:,1) = X(:,1)         + 2*mean(Y(:,J1).^2,2);
        PopObj(:,2) = 1-sqrt(X(:,1)) + 2*mean(Y(:,J2).^2,2);
    case 'UF2'
        D = size(X,2);
        J1 = 3:2:D;
        J2 = 2:2:D;
        Y  = zeros(size(X));
        X1 = repmat(X(:,1),1,length(J1));
        Y(:,J1) = X(:,J1)-(0.3*X1.^2.*cos(24*pi*X1+4*repmat(J1,size(X,1),1)*pi/D)+0.6*X1).*cos(6*pi*X1+repmat(J1,size(X,1),1)*pi/D);
        X1 = repmat(X(:,1),1,length(J2));
        Y(:,J2) = X(:,J2)-(0.3*X1.^2.*cos(24*pi*X1+4*repmat(J2,size(X,1),1)*pi/D)+0.6*X1).*sin(6*pi*X1+repmat(J2,size(X,1),1)*pi/D);
        PopObj(:,1) = X(:,1)         + 2*mean(Y(:,J1).^2,2);
        PopObj(:,2) = 1-sqrt(X(:,1)) + 2*mean(Y(:,J2).^2,2);
    case 'UF3'
        D = size(X,2);
        J1 = 3:2:D;
        J2 = 2:2:D;
        Y  = X-repmat(X(:,1),1,D).^(0.5*(1+3*(repmat(1:D,size(X,1),1)-2)/(D-2)));
        PopObj(:,1) = X(:,1)         + 2/length(J1)*(4*sum(Y(:,J1).^2,2)-2*prod(cos(20*Y(:,J1)*pi./sqrt(repmat(J1,size(X,1),1))),2)+2);
        PopObj(:,2) = 1-sqrt(X(:,1)) + 2/length(J2)*(4*sum(Y(:,J2).^2,2)-2*prod(cos(20*Y(:,J2)*pi./sqrt(repmat(J2,size(X,1),1))),2)+2);
    case 'UF4'
        D = size(X,2);
        J1 = 3:2:D;
        J2 = 2:2:D;
        Y  = X - sin(6*pi*repmat(X(:,1),1,D)+repmat(1:D,size(X,1),1)*pi/D);
        hY = abs(Y)./(1+exp(2*abs(Y)));
        PopObj(:,1) = X(:,1)      + 2*mean(hY(:,J1),2);
        PopObj(:,2) = 1-X(:,1).^2 + 2*mean(hY(:,J2),2);
    case 'UF5'
        D = size(X,2);
        J1 = 3:2:D;
        J2 = 2:2:D;
        Y  = X - sin(6*pi*repmat(X(:,1),1,D)+repmat(1:D,size(X,1),1)*pi/D);
        hY = 2*Y.^2-cos(4*pi*Y)+1;
        PopObj(:,1) = X(:,1)   + (1/20+0.1)*abs(sin(20*pi*X(:,1)))+2*mean(hY(:,J1),2);
        PopObj(:,2) = 1-X(:,1) + (1/20+0.1)*abs(sin(20*pi*X(:,1)))+2*mean(hY(:,J2),2);
    case 'UF6'
        D = size(X,2);
        J1 = 3:2:D;
        J2 = 2:2:D;
        Y  = X - sin(6*pi*repmat(X(:,1),1,D)+repmat(1:D,size(X,1),1)*pi/D);
        PopObj(:,1) = X(:,1)   + max(0,2*(1/4+0.1)*sin(4*pi*X(:,1)))+2/length(J1)*(4*sum(Y(:,J1).^2,2)-2*prod(cos(20*Y(:,J1)*pi./sqrt(repmat(J1,size(X,1),1))),2)+2);
        PopObj(:,2) = 1-X(:,1) + max(0,2*(1/4+0.1)*sin(4*pi*X(:,1)))+2/length(J2)*(4*sum(Y(:,J2).^2,2)-2*prod(cos(20*Y(:,J2)*pi./sqrt(repmat(J2,size(X,1),1))),2)+2);
    case 'UF7'
        D = size(X,2);
        J1 = 3:2:D;
        J2 = 2:2:D;
        Y  = X - sin(6*pi*repmat(X(:,1),1,D)+repmat(1:D,size(X,1),1)*pi/D);
        PopObj(:,1) = X(:,1).^0.2   + 2*mean(Y(:,J1).^2,2);
        PopObj(:,2) = 1-X(:,1).^0.2 + 2*mean(Y(:,J2).^2,2);
    case 'UF8'
        D = size(X,2);
        J1 = 4:3:D;
        J2 = 5:3:D;
        J3 = 3:3:D;
        Y  = X-2*repmat(X(:,2),1,D).*sin(2*pi*repmat(X(:,1),1,D)+repmat(1:D,size(X,1),1)*pi/D);
        PopObj(:,1) = cos(0.5*X(:,1)*pi).*cos(0.5*X(:,2)*pi) + 2*mean(Y(:,J1).^2,2);
        PopObj(:,2) = cos(0.5*X(:,1)*pi).*sin(0.5*X(:,2)*pi) + 2*mean(Y(:,J2).^2,2);
        PopObj(:,3) = sin(0.5*X(:,1)*pi)                     + 2*mean(Y(:,J3).^2,2);
    case 'UF9'
        D = size(X,2);
        J1 = 4:3:D;
        J2 = 5:3:D;
        J3 = 3:3:D;
        Y  = X-2*repmat(X(:,2),1,D).*sin(2*pi*repmat(X(:,1),1,D)+repmat(1:D,size(X,1),1)*pi/D);
        PopObj(:,1) = 0.5*(max(0,1.1*(1-4*(2*X(:,1)-1).^2))+2*X(:,1)).*X(:,2)   + 2*mean(Y(:,J1).^2,2);
        PopObj(:,2) = 0.5*(max(0,1.1*(1-4*(2*X(:,1)-1).^2))-2*X(:,1)+2).*X(:,2) + 2*mean(Y(:,J2).^2,2);
        PopObj(:,3) = 1-X(:,2)                                                  + 2*mean(Y(:,J3).^2,2);
        
%% 三目标目标函数
    case 'UF10'
%% 基本参数输入

%引入其他函数
global Ppv Pwt Qhpsh Qhpsc Php

%加载负荷文件
load CoolingFuhe.txt
load HeatingFuhe.txt
load Fuhe.txt
load Temperature.txt
load DianJia.txt

%输入设备基本参数
SOCmax=0.9;
SOCmin=0.1;
Eb=50;                                                       %蓄电池额定容量
X=fix(X);
SOC(1,1)=0.9;
etaeh=0.95;
Qge=6;
Cse=0;

%全寿命周期经济成本（60年，GHIES，yuan）
Cpv=600*2.4+2.54*60;                                           %光伏电池成本  购买成本+运行维护成本
Cwt=12720*30*3+65*30*60;                                       %风力发电成本
Cbss=100000*5+15.9*50*60;                                      %蓄电池成本
Chp=240000*4+1200*60;                                          %热泵机组成本
Cge=20000;                                                     %地热打井成本 
Ctes=2360*3;                                                     %储冷/热水箱成本

%全寿命碳排放成本（60年，kg）
COpv=1212.90;                                               %光伏电池成本  购买成本+运行维护碳排放成本
COwt=28695.55;                                              %风力发电碳排放成本
CObss=51450;                                                %蓄电池碳排放成本
COhp=355439.5;                                              %热泵机组碳排放成本
COge=7785.5;                                                %地热打井碳排放成本 
COtes=2738.65;                                               %储热水箱碳排放


Ttessdh=53;                                          %水箱供暖时启停温度
Ttessdhh=50;
Ttessdc=2;
Ttessdcc=5;
Vtes=2.36;                                              %水箱参数
Ates=10.99;
Ktes=0.0012;
mtes=Vtes*997;
Cw=4.18;
Tbui=45;
Tbuc=10;
 
%预分配
LCC1=zeros(200,1);
Cs=zeros(200,1);
LCC=zeros(200,1);
CO2ot=zeros(200,1);
CO2op=zeros(200,1);
CO2=zeros(200,1);
EER=zeros(200,1);
Qloss=zeros(200,1);
PopObj=zeros(200,4);
Xhp=zeros(8760,1);
Xhp01=zeros(8760,1);
Xhp02=zeros(8760,1);
Phppl=zeros(8760,1);
PL=zeros(8760,1);
elecnum=zeros(8760,1);
Pbss=zeros(8760,1);
Xge=zeros(8760,1);
Xge01=zeros(8760,1);
Qdeloss=zeros(8760,1);
%Qdelossc=zeros(8760,1);
Celebuy=zeros(8760,1);
Celesell=zeros(8760,1);
Qhptes=zeros(8760,1);
Qhpbui=zeros(8760,1);
Qtes=zeros(8760,1);
COPpl=zeros(8760,1);
COPnom=zeros(8760,1);
Xloss=zeros(8760,1);
Phptes=zeros(8760,1);
Ttes=zeros(8761,1);
%Ttes(1,1)=50;
%Ttes(3625,1)=0;
%Ttes(7632,1)=50;
More=zeros(8760,1);
%Emore=zeros(8760,1); 
Qhpscpl=zeros(8760,1);
Qhpshpl=zeros(8760,1);


%% 计算解集主循环

for j= 1:200
    
    %%--单行粒子解集计算--%%
    
    %% 供暖季part01, 1-2520h
    
    for i=1:2520
        
        a(i,1)=i;
        b(i,1)=a(i,1)/24;
        c(i,1)=fix(b(i,1));
        d(i,1)=b(i,1)-c(i,1);
        
        %%%---23.00-07.00,峰谷时间，热泵向水箱供热---%%%
        if d(i,1)<0.31
            
            if Ttes(i,1)<Ttessdhh                                                                                            %水箱温度控制启停判断，水箱温度<53
                Qhptes(i,1)=((Ttessdh-Ttes(i,1))*Cw*X(j,5)*mtes+X(j,5)*Ktes*Ates*(Ttes(i,1)-Temperature(i,1)))/3600;                %热泵峰谷时期供热热负荷
                
                %%---热泵部分负荷及缺失计算---%%
                
                %第一次计算热泵数量
                Xhp01(i,1)=Qhptes(i,1)/Qhpsh(i,1);
                
                if Xhp01(i,1)<1
                    Xhp02(i,1)=0;
                else
                    Xhp02(i,1)=fix(Xhp01(i,1))+1;
                end
                
                Xhp(i,1)=Xhp02(i,1);
                
                %第二次计算热泵数量，根据边界修正,并计算缺失量与部分负荷下的热泵功耗、COP
                if Xhp(i,1)<=X(j,4) && Xhp(i,1)>0
                    PL(i,1)=Qhptes(i,1)/(Qhpsh(i,1)*Xhp02(i,1));                          %情况1，存在部分符合下的热泵
                    COPnom(i,1)=Qhpsh(i,1)/Php(i,1);
                    COPpl(i,1)=COPnom(i,1)*(-0.00006*PL(i,1)^6+0.0017*PL(i,1)^5-0.0181*PL(i,1)^4+0.096*PL(i,1)^3-0.2697*PL(i,1)^2+0.4276*PL(i,1)+0.7626);
                    Phppl(i,1)=Qhptes(i,1)/COPpl(i,1);                                    %计算部分负荷下该热泵的功耗
                    Qhpscpl(i,1)=Qhptes(i,1)-Phppl(i,1);
                    Ttes(i+1,1)=Ttessdh;
                elseif Xhp(i,1)>X(j,4)
                    Xhp02(i,1)=fix(X(j,4));
                    %Qdeloss(i,1)=Qhptes(i,1)-Qhpsh(i,1)*Xhp02(i,1);                          %Qdeloss代表缺失的热量
                    %Xloss(i,1)=Qdeloss(i,1)/Qhptes(i,1);
                    Phppl(i,1)=Php(i,1)*Xhp02(i,1);
                    Qhpscpl(i,1)=Xhp02(i,1)*Qhpsc(i,1);
                    Ttes(i+1,1)=(Xhp02(i,1)*Qhpsh(i,1)*3600-X(j,5)*Ktes*Ates*(Ttes(i,1)-Temperature(i,1)))/(mtes*Cw*X(j,5))+Ttes(i,1);
                else
                    Phppl(i,1)=0;
                    Qhpscpl(i,1)=0;
                    Ttes(i+1,1)=Ktes*Ates*(Temperature(i,1)-Ttes(i,1))/(mtes*Cw)+Ttes(i,1);
                end                                                                      %Xhp表示最终的热泵数量，带小数点，情况2 X(j,4)表示数量，即边界值
                
                Phptes(i,1)=Phppl(i,1);
                
                %计算地埋管数量
                Xge01(i,1)=Qhpscpl(i,1)/(6*0.95);                                        %地埋管传热效率设置为95%
                Xge(i,1)=fix(Xge01(i,1));
                if Xge01(i,1)-Xge(i,1) > 0                                                 %对地埋管数量取整，向上+1
                    Xge(i,1)=fix(Xge01(i,1))+1;
                else
                    Xge(i,1)=Xge01(i,1);
                end
                
            else
                Ttes(i+1,1)=Ktes*Ates*(Temperature(i,1)-Ttes(i,1))/(mtes*Cw)+Ttes(i,1);
            end
            
            
            %%%---电负荷计算---%%%
            Pbss(i,1)=Fuhe(i,1)-Ppv(i,1)*X(j,1)-Pwt(i,1)*X(j,2);
            
            %---蓄电池充放电状态判断---%
            if Pbss(i,1) < 0
                SOC(i+1,1)=(1-0.03)*SOC(i,1)-Pbss(i,1)*0.9/(X(j,3)*Eb);
                More(i,1)=abs(Pbss(i,1));
                elecnum(i,1)=0;
            elseif   Pbss(i,1) > 0
                SOC(i+1,1)=(1-0.03)*SOC(i,1)-Pbss(i,1)/(X(j,3)*0.9*Eb);
                More(i,1)=0;
                elecnum(i,1)=abs(Pbss(i,1));
            else
                SOC(i+1,1)=SOC(i,1);
                More(i,1)=0;
                elecnum(i,1)=0;
            end
            %---蓄电池充放电状态判断结束---%
            
            %---多余电量判断---%
            if  SOC(i+1,1)>=SOCmax
                SOC(i+1,1)=SOCmax;
                %先计算即将溢出时的蓄电池真实充电功率
                if  SOC(i,1)~=SOCmax
                    Pbss(i,1)=(SOC(i+1,1)-(1-0.03)*SOC(i,1))*Eb*X(j,3)/(-0.9);
                else
                    Pbss(i,1)=0;
                end
                %记录多余电量
                %在100%的情况下为弃光与弃风
                if More(i,1)>=abs(Pbss(i,1))
                    More(i,1)=More(i,1)-abs(Pbss(i,1));
                else
                    More(i,1)=0;
                end
                %记录多余电量
            else
                More(i,1)=0;
            end
            %---多余电量判断结束---%
            
            %---缺电判断---%
            if  SOC(i+1,1)<=SOCmin
                SOC(i+1,1)=SOCmin;
                %计算即将见底时的蓄电池真实放电功率
                if  SOC(i,1)~=SOCmin
                    Pbss(i,1)=(SOC(i+1,1)-(1-0.03)*SOC(i,1))*Eb*(-0.9)*X(j,3);
                else
                    Pbss(i,1)=0;
                end
                %在见底时的真实缺电率
                if elecnum(i,1)>=abs(Pbss(i,1))
                    elecnum(i,1)=elecnum(i,1)-abs(Pbss(i,1));
                else
                    elecnum(i,1)=0;
                end
            else
                elecnum(i,1)=0;
            end
            
            %外购电耗费
            Celebuy(i,1)=Phptes(i,1)*DianJia(i,1);
            
            %售电盈利
            Celesell(i,1)=More(i,1)*Cse;
            %Emore(i,1)=More(i,1)/(Ppv(i,1)*X(j,1)+Pwt(i,1)*X(j,2)+0.001);
            
            %---缺电判断结束---%
            %%%---电模块结束---%%%
        end
        
        %%%---23.00-07.00计算结束---%%%
        
        %%%---07.00-08.00计算，非上班时间，非峰谷时间---%%%
        if d(i,1)>0.31 && d(i,1)<0.36
            
            Ttes(i+1,1)=Ktes*Ates*(Temperature(i,1)-Ttes(i,1))/(mtes*Cw)+Ttes(i,1);
            
            %%%---电负荷计算---%%%
            Pbss(i,1)=Fuhe(i,1)-Ppv(i,1)*X(j,1)-Pwt(i,1)*X(j,2);
            
            %---蓄电池充放电状态判断---%
            if Pbss(i,1) < 0
                SOC(i+1,1)=(1-0.03)*SOC(i,1)-Pbss(i,1)*0.9/(X(j,3)*Eb);
                More(i,1)=abs(Pbss(i,1));
                elecnum(i,1)=0;
            elseif   Pbss(i,1) > 0
                SOC(i+1,1)=(1-0.03)*SOC(i,1)-Pbss(i,1)/(X(j,3)*0.9*Eb);
                More(i,1)=0;
                elecnum(i,1)=abs(Pbss(i,1));
            else
                SOC(i+1,1)=SOC(i,1);
                More(i,1)=0;
                elecnum(i,1)=0;
            end
            %---蓄电池充放电状态判断结束---%
            
            %---多余电量判断---%
            if  SOC(i+1,1)>=SOCmax
                SOC(i+1,1)=SOCmax;
                %先计算即将溢出时的蓄电池真实充电功率
                if  SOC(i,1)~=SOCmax
                    Pbss(i,1)=(SOC(i+1,1)-(1-0.03)*SOC(i,1))*Eb*X(j,3)/(-0.9);
                else
                    Pbss(i,1)=0;
                end
                %记录多余电量
                %在100%的情况下为弃光与弃风
                if More(i,1)>=abs(Pbss(i,1))
                    More(i,1)=More(i,1)-abs(Pbss(i,1));
                else
                    More(i,1)=0;
                end
                %记录多余电量
            else
                More(i,1)=0;
            end
            %---多余电量判断结束---%
            
            %---缺电判断---%
            if  SOC(i+1,1)<=SOCmin
                SOC(i+1,1)=SOCmin;
                %计算即将见底时的蓄电池真实放电功率
                if  SOC(i,1)~=SOCmin
                    Pbss(i,1)=(SOC(i+1,1)-(1-0.03)*SOC(i,1))*Eb*(-0.9)*X(j,3);
                else
                    Pbss(i,1)=0;
                end
                %在见底时的真实缺电率
                if elecnum(i,1)>=abs(Pbss(i,1))
                    elecnum(i,1)=elecnum(i,1)-abs(Pbss(i,1));
                else
                    elecnum(i,1)=0;
                end
            else
                elecnum(i,1)=0;
            end
            
            %外购电耗费
            Celebuy(i,1)=0;
            
            %售电盈利
            Celesell(i,1)=More(i,1)*Cse;
            %Emore(i,1)=More(i,1)/(Ppv(i,1)*X(j,1)+Pwt(i,1)*X(j,2)+0.001);
            
        end
        %%%---07.00-08.00计算结束---%%%
        
        %%%---8.00-9.00,工作时间，平常电价---%%%
        if d(i,1)>0.36 && d(i,1)<0.39
            
            %%--水箱温度计算--%%
            Ttes(i+1,1)=Ktes*Ates*(Temperature(i,1)-Ttes(i,1))/(mtes*Cw)+Ttes(i,1);
            %%--水箱温度计算结束--%%
            
            %%--热泵计算--%%
            %第一次计算热泵数量
            Xhp01(i,1)=HeatingFuhe(i,1)/Qhpsh(i,1);
            
            if Xhp01(i,1)<1
                Xhp02(i,1)=0;
            else
                Xhp02(i,1)=fix(Xhp01(i,1))+1;
            end
            
            Xhp(i,1)=Xhp02(i,1);
                
            %第二次计算热泵数量，根据边界修正,并计算缺失量与部分负荷下的热泵功耗、COP
            if Xhp(i,1)<=X(j,4) && Xhp(i,1)>0
                PL(i,1)=HeatingFuhe(i,1)/(Qhpsh(i,1)*Xhp02(i,1));                          %情况1，存在部分符合下的热泵
                COPnom(i,1)=Qhpsh(i,1)/Php(i,1);
                COPpl(i,1)=COPnom(i,1)*(-0.00006*PL(i,1)^6+0.0017*PL(i,1)^5-0.0181*PL(i,1)^4+0.096*PL(i,1)^3-0.2697*PL(i,1)^2+0.4276*PL(i,1)+0.7626);
                Phppl(i,1)=HeatingFuhe(i,1)/COPpl(i,1);
                Qhpscpl(i,1)=HeatingFuhe(i,1)-Phppl(i,1);
            elseif Xhp(i,1)>X(j,4)
                Xhp02(i,1)=fix(X(j,4));
                Qdeloss(i,1)=HeatingFuhe(i,1)-Qhpsh(i,1)*Xhp02(i,1);                          %Qdeloss代表缺失的热量
                Xloss(i,1)=Qdeloss(i,1)/HeatingFuhe(i,1);
                Phppl(i,1)=Php(i,1)*Xhp02(i,1);
                Qhpscpl(i,1)=Xhp02(i,1)*Qhpsc(i,1);
            else
                Qdeloss(i,1)=HeatingFuhe(i,1);                          %Qdeloss代表缺失的热量
                Xloss(i,1)=1;
                Phppl(i,1)=0;
                Qhpscpl(i,1)=0;
            end                                                                      %Xhp表示最终的热泵数量，带小数点，情况2 X(j,4)表示数量，即边界值
            
            %计算地埋管的数量
            Xge01(i,1)=Qhpscpl(i,1)/(6*0.95);
            Xge(i,1)=fix(Xge(i,1));
            
            if Xge01(i,1)-Xge(i,1) > 0                                                %对地埋管数量取整，向上+1
                Xge(i,1)=fix(Xge01(i,1))+1;
            else
                Xge(i,1)=Xge01(i,1);
            end
            %%--热泵计算结束--%%
            
            %%--电负荷计算--%%
            Pbss(i,1)=Fuhe(i,1)+Phppl(i,1)-Ppv(i,1)*X(j,1)-Pwt(i,1)*X(j,2);
            
            %蓄电池充放电状态判断
            if Pbss(i,1) < 0                                                         %pbss为负，代表要充电，此时的pbss代表不考虑上下限的功率
                SOC(i+1,1)=(1-0.03)*SOC(i,1)-Pbss(i,1)*0.9/(X(j,3)*Eb);
                More(i,1)=abs(Pbss(i,1));                                            %记录多余的电量
                elecnum(i,1)=0;
            elseif Pbss(i,1) > 0                                                     %pbss为正，代表要放电
                SOC(i+1,1)=(1-0.03)*SOC(i,1)-Pbss(i,1)/(X(j,3)*0.9*Eb);
                More(i,1)=0;                                                         %记录多余的浪费电量
                elecnum(i,1)=abs(Pbss(i,1));
            else
                SOC(i+1,1)=SOC(i,1);
                More(i,1)=0;                                                         %记录多余的浪费电量
                elecnum(i,1)=0;
            end
            
            %多余电量判断
            if SOC(i+1,1)>=SOCmax
                SOC(i+1,1)=SOCmax;
                %先计算即将溢出时的蓄电池真实充电功率
                if SOC(i,1)~=SOCmax                                                      %~=是不等于的意思
                    Pbss(i,1)=(SOC(i+1,1)-(1-0.03)*SOC(i,1))*Eb*X(j,3)/(-0.9);
                else
                    Pbss(i,1)=0;
                end
                %在100%的情况下为弃光与弃风
                if More(i,1)>=abs(Pbss(i,1))
                    More(i,1)=More(i,1)-abs(Pbss(i,1));
                else
                    More(i,1)=0;
                end
            else
                More(i,1)=0;
            end
            
            %缺电判断
            if SOC(i+1,1)<=SOCmin
                SOC(i+1,1)=SOCmin;
                %计算即将见底时的蓄电池真实放电功率
                if SOC(i,1)~=SOCmin
                    Pbss(i,1)=(SOC(i+1,1)-(1-0.03)*SOC(i,1))*Eb*(-0.9)*X(j,3);
                else
                    Pbss(i,1)=0;
                end
                %真实缺电率
                if elecnum(i,1)>=abs(Pbss(i,1))
                    elecnum(i,1)=elecnum(i,1)-abs(Pbss(i,1));
                else
                    elecnum(i,1)=0;
                end
            else
                elecnum(i,1)=0;
            end
            
            %外购电耗费
            Celebuy(i,1)=elecnum(i,1)*DianJia(i,1);
            
            %售电盈利
            Celesell(i,1)=More(i,1)*Cse;
            %Emore(i,1)=More(i,1)/(Ppv(i,1)*X(j,1)+Pwt(i,1)*X(j,2)+0.001);
            %%--电负荷计算结束--%%
            
        end
        %%%---8.00-9.00计算结束---%%%
        
        %%%---9.00-12.00,工作时间，峰值电价，开启水箱---%%%
        if d(i,1)>0.39 && d(i,1)<0.51
            
            %%--热泵计算--%%
            Qtes(i,1)=((Ttes(i,1)-Tbui)*Cw*X(j,5)*mtes-X(j,5)*Ktes*Ates*(Ttes(i,1)-Temperature(i,1)))/3600;
            
            %防止后续增加热泵负荷
            if Qtes(i,1)<= 0
                Qtes(i,1)=0;
                
                %%--水箱温度计算--%%
                Ttes(i+1,1)=Ktes*Ates*(Temperature(i,1)-Ttes(i,1))/(mtes*Cw)+Ttes(i,1);
                %%--水箱温度计算结束--%%
            
                %%--热泵计算--%%
                %第一次计算热泵数量
                Xhp01(i,1)=HeatingFuhe(i,1)/Qhpsh(i,1);
                
                if Xhp01(i,1)<1
                    Xhp02(i,1)=0;
                else
                    Xhp02(i,1)=fix(Xhp01(i,1))+1;
                end
                
                Xhp(i,1)=Xhp02(i,1);
                
                %第二次计算热泵数量，根据边界修正,并计算缺失量与部分负荷下的热泵功耗、COP
                if Xhp(i,1)<=X(j,4) && Xhp(i,1)>0
                    PL(i,1)=HeatingFuhe(i,1)/(Qhpsh(i,1)*Xhp02(i,1));                          %情况1，存在部分符合下的热泵
                    COPnom(i,1)=Qhpsh(i,1)/Php(i,1);
                    COPpl(i,1)=COPnom(i,1)*(-0.00006*PL(i,1)^6+0.0017*PL(i,1)^5-0.0181*PL(i,1)^4+0.096*PL(i,1)^3-0.2697*PL(i,1)^2+0.4276*PL(i,1)+0.7626);
                    Phppl(i,1)=HeatingFuhe(i,1)/COPpl(i,1);
                    Qhpscpl(i,1)=HeatingFuhe(i,1)-Phppl(i,1);
                elseif Xhp(i,1)>X(j,4)
                    Xhp02(i,1)=fix(X(j,4));
                    Qdeloss(i,1)=HeatingFuhe(i,1)-Qhpsh(i,1)*Xhp02(i,1);                          %Qdeloss代表缺失的热量
                    Xloss(i,1)=Qdeloss(i,1)/HeatingFuhe(i,1);
                    Phppl(i,1)=Php(i,1)*Xhp02(i,1);
                    Qhpscpl(i,1)=Xhp02(i,1)*Qhpsc(i,1);
                else
                    Qdeloss(i,1)=HeatingFuhe(i,1);                          %Qdeloss代表缺失的热量
                    Xloss(i,1)=1;
                    Phppl(i,1)=0;
                    Qhpscpl(i,1)=0;
                end                                                                      %Xhp表示最终的热泵数量，带小数点，情况2 X(j,4)表示数量，即边界值
            
                %计算地埋管的数量
                Xge01(i,1)=Qhpscpl(i,1)/(6*0.95);
                Xge(i,1)=fix(Xge(i,1));
            
                if Xge01(i,1)-Xge(i,1) > 0                                                %对地埋管数量取整，向上+1
                    Xge(i,1)=fix(Xge01(i,1))+1;
                else
                    Xge(i,1)=Xge01(i,1);
                end
                %%--热泵计算结束--%%
            end
            
            %情况1/2/3
            if Qtes(i,1) >= HeatingFuhe(i,1)
                Ttes(i+1,1)=(X(j,5)*Ktes*Ates*(Temperature(i,1)-Ttes(i,1))-HeatingFuhe(i,1)*3600)/(X(j,5)*mtes*Cw)+Ttes(i,1);
                Xhp02(i,1)=0;
            elseif Qtes(i,1)<HeatingFuhe(i,1) && Qtes(i,1)>0
                Qhpbui(i,1)=HeatingFuhe(i,1)-Qtes(i,1);
                Ttes(i+1,1)=Tbui;
                
                %情况2/3
                %第一次计算热泵数量
                Xhp01(i,1)=Qhpbui(i,1)/Qhpsh(i,1);
                
                if Xhp01(i,1)<1
                    Xhp02(i,1)=0;
                else
                    Xhp02(i,1)=fix(Xhp01(i,1))+1;
                end
                
                Xhp(i,1)=Xhp02(i,1);
                
                %第二次计算热泵数量，根据边界修正,并计算缺失量与部分负荷下的热泵功耗、COP
                if Xhp(i,1) <= X(j,4) && Xhp(i,1)>0
                    PL(i,1)=Qhpbui(i,1)/(Qhpsh(i,1)*Xhp02(i,1));                          %情况1，存在部分符合下的热泵
                    COPnom(i,1)=Qhpsh(i,1)/Php(i,1);
                    COPpl(i,1)=COPnom(i,1)*(-0.00006*PL(i,1)^6+0.0017*PL(i,1)^5-0.0181*PL(i,1)^4+0.096*PL(i,1)^3-0.2697*PL(i,1)^2+0.4276*PL(i,1)+0.7626);
                    Phppl(i,1)=Qhpbui(i,1)/COPpl(i,1);
                    Qhpscpl(i,1)=Qhpbui(i,1)-Phppl(i,1);
                elseif Xhp(i,1)> X(j,4)
                    Xhp02(i,1)=fix(X(j,4));
                    Qdeloss(i,1)=Qhpbui(i,1)-Qhpsh(i,1)*Xhp02(i,1);                          %Qdeloss代表缺失的热量
                    Xloss(i,1)=Qdeloss(i,1)/Qhpbui(i,1);
                    Phppl(i,1)=Php(i,1)*Xhp02(i,1);
                    Qhpscpl(i,1)=Xhp02(i,1)*Qhpsc(i,1);
                else
                    Qdeloss(i,1)=Qhpbui(i,1);                          %Qdeloss代表缺失的热量
                    Xloss(i,1)=1;
                    Phppl(i,1)=0;
                    Qhpscpl(i,1)=0;
                end                                                                      %Xhp表示最终的热泵数量，带小数点，情况2 X(j,4)表示数量，即边界值
                
                %计算地埋管数量
                Xge01(i,1)=Qhpscpl(i,1)/(6*0.95);
                Xge(i,1)=fix(Xge(i,1));
                
                if Xge01(i,1)-Xge(i,1) > 0                                               %对地埋管数量取整，向上+1
                    Xge(i,1)=fix(Xge01(i,1))+1;
                else
                    Xge(i,1)=Xge01(i,1);
                end
                
            end
            %%--热泵计算结束--%%
            
            %%--电负荷计算--%%
            Pbss(i,1)=Fuhe(i,1)+Phppl(i,1)-Ppv(i,1)*X(j,1)-Pwt(i,1)*X(j,2);
            
            %蓄电池充放电状态判断
            if Pbss(i,1) < 0                                                            %pbss为负，代表要充电，此时的pbss代表不考虑上下限的功率
                SOC(i+1,1)=(1-0.03)*SOC(i,1)-Pbss(i,1)*0.9/(X(j,3)*Eb);
                More(i,1)=abs(Pbss(i,1));                                               %记录多余的电量
                elecnum(i,1)=0;
            elseif Pbss(i,1) > 0                                                        %pbss为正，代表要放电
                SOC(i+1,1)=(1-0.03)*SOC(i,1)-Pbss(i,1)/(X(j,3)*0.9*Eb);
                More(i,1)=0;                                                            %记录多余的浪费电量
                elecnum(i,1)=abs(Pbss(i,1));
            else
                SOC(i+1,1)=SOC(i,1);
                More(i,1)=0;                                                            %记录多余的浪费电量
                elecnum(i,1)=0;
            end
            
            %多余电量判断
            if SOC(i+1,1)>=SOCmax
                SOC(i+1,1)=SOCmax;
                %先计算即将溢出时的蓄电池真实充电功率
                if SOC(i,1)~=SOCmax                                                         %~=是不等于的意思
                    Pbss(i,1)=(SOC(i+1,1)-(1-0.03)*SOC(i,1))*Eb*X(j,3)/(-0.9);
                else
                    Pbss(i,1)=0;
                end
                %在100%的情况下为弃光与弃风
                if More(i,1)>=abs(Pbss(i,1))
                    More(i,1)=More(i,1)-abs(Pbss(i,1));
                else
                    More(i,1)=0;
                end
                %记录多余电量
            else
                More(i,1)=0;
            end
            
            %缺电判断
            if SOC(i+1,1)<=SOCmin
                SOC(i+1,1)=SOCmin;
                %计算即将见底时的蓄电池真实放电功率
                if SOC(i,1)~=SOCmin
                    Pbss(i,1)=(SOC(i+1,1)-(1-0.03)*SOC(i,1))*Eb*(-0.9)*X(j,3);
                else
                    Pbss(i,1)=0;
                end
                %真实缺电率
                if elecnum(i,1)>=abs(Pbss(i,1))
                    elecnum(i,1)=elecnum(i,1)-abs(Pbss(i,1));
                else
                    elecnum(i,1)=0;
                end
            else
                elecnum(i,1)=0;
            end
            
            %外购电耗费
            Celebuy(i,1)=elecnum(i,1)*DianJia(i,1);
            
            %售电盈利
            Celesell(i,1)=More(i,1)*Cse;
            %Emore(i,1)=More(i,1)/(Ppv(i,1)*X(j,1)+Pwt(i,1)*X(j,2)+0.001);
            %%--电负荷计算结束--%%
            
        end
        %%%---9.00-12.00计算结束---%%%
        
        %%%---12.00-16.00，工作时间，平常电价---%%%
        if d(i,1)>0.51 && d(i,1)<0.69
            
            %%--水箱温度计算--%%
            Ttes(i+1,1)=Ktes*Ates*(Temperature(i,1)-Ttes(i,1))/(mtes*Cw)+Ttes(i,1);
            %%--水箱温度计算结束--%%
            
            %%--热泵计算--%%
            %第一次计算热泵数量
            Xhp01(i,1)=HeatingFuhe(i,1)/Qhpsh(i,1);
            
            if Xhp01(i,1)<1
                Xhp02(i,1)=0;
            else
                Xhp02(i,1)=fix(Xhp01(i,1))+1;
            end
            
            Xhp(i,1)=Xhp02(i,1);
                
            %第二次计算热泵数量，根据边界修正,并计算缺失量与部分负荷下的热泵功耗、COP
            if Xhp(i,1)<=X(j,4) && Xhp(i,1)>0
                PL(i,1)=HeatingFuhe(i,1)/(Qhpsh(i,1)*Xhp02(i,1));                          %情况1，存在部分符合下的热泵
                COPnom(i,1)=Qhpsh(i,1)/Php(i,1);
                COPpl(i,1)=COPnom(i,1)*(-0.00006*PL(i,1)^6+0.0017*PL(i,1)^5-0.0181*PL(i,1)^4+0.096*PL(i,1)^3-0.2697*PL(i,1)^2+0.4276*PL(i,1)+0.7626);
                Phppl(i,1)=HeatingFuhe(i,1)/COPpl(i,1);
                Qhpscpl(i,1)=HeatingFuhe(i,1)-Phppl(i,1);
            elseif Xhp(i,1)>X(j,4)
                Xhp02(i,1)=fix(X(j,4));
                Qdeloss(i,1)=HeatingFuhe(i,1)-Qhpsh(i,1)*Xhp02(i,1);                          %Qdeloss代表缺失的热量
                Xloss(i,1)=Qdeloss(i,1)/HeatingFuhe(i,1);
                Phppl(i,1)=Php(i,1)*Xhp02(i,1);
                Qhpscpl(i,1)=Xhp02(i,1)*Qhpsc(i,1);
            else
                Qdeloss(i,1)=HeatingFuhe(i,1);                          %Qdeloss代表缺失的热量
                Xloss(i,1)=1;
                Phppl(i,1)=0;
                Qhpscpl(i,1)=0;
            end                                                                      %Xhp表示最终的热泵数量，带小数点，情况2 X(j,4)表示数量，即边界值
            
            %计算地埋管的数量
            Xge01(i,1)=Qhpscpl(i,1)/(6*0.95);
            Xge(i,1)=fix(Xge(i,1));
            
            if Xge01(i,1)-Xge(i,1) > 0                                                %对地埋管数量取整，向上+1
                Xge(i,1)=fix(Xge01(i,1))+1;
            else
                Xge(i,1)=Xge01(i,1);
            end
            %%--热泵计算结束--%%
            
            %%--电负荷计算--%%
            Pbss(i,1)=Fuhe(i,1)+Phppl(i,1)-Ppv(i,1)*X(j,1)-Pwt(i,1)*X(j,2);
            
            %蓄电池充放电状态判断
            if Pbss(i,1) < 0                                                         %pbss为负，代表要充电，此时的pbss代表不考虑上下限的功率
                SOC(i+1,1)=(1-0.03)*SOC(i,1)-Pbss(i,1)*0.9/(X(j,3)*Eb);
                More(i,1)=abs(Pbss(i,1));                                            %记录多余的电量
                elecnum(i,1)=0;
            elseif Pbss(i,1) > 0                                                     %pbss为正，代表要放电
                SOC(i+1,1)=(1-0.03)*SOC(i,1)-Pbss(i,1)/(X(j,3)*0.9*Eb);
                More(i,1)=0;                                                         %记录多余的浪费电量
                elecnum(i,1)=abs(Pbss(i,1));
            else
                SOC(i+1,1)=SOC(i,1);
                More(i,1)=0;                                                         %记录多余的浪费电量
                elecnum(i,1)=0;
            end
            
            %多余电量判断
            if SOC(i+1,1)>=SOCmax
                SOC(i+1,1)=SOCmax;
                %先计算即将溢出时的蓄电池真实充电功率
                if SOC(i,1)~=SOCmax                                                      %~=是不等于的意思
                    Pbss(i,1)=(SOC(i+1,1)-(1-0.03)*SOC(i,1))*Eb*X(j,3)/(-0.9);
                else
                    Pbss(i,1)=0;
                end
                %在100%的情况下为弃光与弃风
                if More(i,1)>=abs(Pbss(i,1))
                    More(i,1)=More(i,1)-abs(Pbss(i,1));
                else
                    More(i,1)=0;
                end
                %记录多余电量
            else
                More(i,1)=0;
            end
            
            %缺电判断
            if SOC(i+1,1)<=SOCmin
                SOC(i+1,1)=SOCmin;
                %计算即将见底时的蓄电池真实放电功率
                if SOC(i,1)~=SOCmin
                    Pbss(i,1)=(SOC(i+1,1)-(1-0.03)*SOC(i,1))*Eb*(-0.9)*X(j,3);
                else
                    Pbss(i,1)=0;
                end
                %真实缺电率
                if elecnum(i,1)>=abs(Pbss(i,1))
                    elecnum(i,1)=elecnum(i,1)-abs(Pbss(i,1));
                else
                    elecnum(i,1)=0;
                end
            else
                elecnum(i,1)=0;
            end
            
            %外购电耗费
            Celebuy(i,1)=elecnum(i,1)*DianJia(i,1);
            
            %售电盈利
            Celesell(i,1)=More(i,1)*Cse;
            %Emore(i,1)=More(i,1)/(Ppv(i,1)*X(j,1)+Pwt(i,1)*X(j,2)+0.001);
            %%--电负荷计算结束--%%
            
        end
        %%%---12.00-16.00计算结束---%%%
        
        %%%---16.00-18.00,工作时间，峰值电价，开启水箱---%%%
        if d(i,1)>0.69 && d(i,1)<0.78
            
            %%--热泵计算--%%
            Qtes(i,1)=((Ttes(i,1)-Tbui)*Cw*X(j,5)*mtes-X(j,5)*Ktes*Ates*(Ttes(i,1)-Temperature(i,1)))/3600;
            
            %防止后续增加热泵负荷
            if Qtes(i,1)<0
                Qtes(i,1)=0;
                
                %%--水箱温度计算--%%
                Ttes(i+1,1)=Ktes*Ates*(Temperature(i,1)-Ttes(i,1))/(mtes*Cw)+Ttes(i,1);
                %%--水箱温度计算结束--%%
            
                %%--热泵计算--%%
                %第一次计算热泵数量
                Xhp01(i,1)=HeatingFuhe(i,1)/Qhpsh(i,1);
                
                if Xhp01(i,1)<1
                    Xhp02(i,1)=0;
                else
                    Xhp02(i,1)=fix(Xhp01(i,1))+1;
                end
                
                Xhp(i,1)=Xhp02(i,1);
                
                %第二次计算热泵数量，根据边界修正,并计算缺失量与部分负荷下的热泵功耗、COP
                if Xhp(i,1)<=X(j,4) && Xhp(i,1)>0
                    PL(i,1)=HeatingFuhe(i,1)/(Qhpsh(i,1)*Xhp02(i,1));                          %情况1，存在部分符合下的热泵
                    COPnom(i,1)=Qhpsh(i,1)/Php(i,1);
                    COPpl(i,1)=COPnom(i,1)*(-0.00006*PL(i,1)^6+0.0017*PL(i,1)^5-0.0181*PL(i,1)^4+0.096*PL(i,1)^3-0.2697*PL(i,1)^2+0.4276*PL(i,1)+0.7626);
                    Phppl(i,1)=HeatingFuhe(i,1)/COPpl(i,1);
                    Qhpscpl(i,1)=HeatingFuhe(i,1)-Phppl(i,1);
                elseif Xhp(i,1)>X(j,4)
                    Xhp02(i,1)=fix(X(j,4));
                    Qdeloss(i,1)=HeatingFuhe(i,1)-Qhpsh(i,1)*Xhp02(i,1);                          %Qdeloss代表缺失的热量
                    Xloss(i,1)=Qdeloss(i,1)/HeatingFuhe(i,1);
                    Phppl(i,1)=Php(i,1)*Xhp02(i,1);
                    Qhpscpl(i,1)=Xhp02(i,1)*Qhpsc(i,1);
                else
                    Qdeloss(i,1)=HeatingFuhe(i,1);                          %Qdeloss代表缺失的热量
                    Xloss(i,1)=1;
                    Phppl(i,1)=0;
                    Qhpscpl(i,1)=0;
                end                                                                      %Xhp表示最终的热泵数量，带小数点，情况2 X(j,4)表示数量，即边界值
            
                %计算地埋管的数量
                Xge01(i,1)=Qhpscpl(i,1)/(6*0.95);
                Xge(i,1)=fix(Xge(i,1));
            
                if Xge01(i,1)-Xge(i,1) > 0                                                %对地埋管数量取整，向上+1
                    Xge(i,1)=fix(Xge01(i,1))+1;
                else
                    Xge(i,1)=Xge01(i,1);
                end
                %%--热泵计算结束--%%
            end
            
            %情况1/2/3
            if Qtes(i,1)>=HeatingFuhe(i,1)
                Ttes(i+1,1)=(X(j,5)*Ktes*Ates*(Temperature(i,1)-Ttes(i,1))-HeatingFuhe(i,1)*3600)/(X(j,5)*mtes*Cw)+Ttes(i,1);
                Xhp02(i,1)=0;
            elseif Qtes(i,1)<HeatingFuhe(i,1) && Qtes(i,1)>0
                Qhpbui(i,1)=HeatingFuhe(i,1)-Qtes(i,1);
                Ttes(i+1,1)=Tbui;
                
                %情况2/3
                %第一次计算热泵数量
                Xhp01(i,1)=Qhpbui(i,1)/Qhpsh(i,1);
                
                if Xhp01(i,1)<1
                    Xhp02(i,1)=0;
                else
                    Xhp02(i,1)=fix(Xhp01(i,1))+1;
                end
                
                Xhp(i,1)=Xhp02(i,1);
                
                %第二次计算热泵数量，根据边界修正,并计算缺失量与部分负荷下的热泵功耗、COP
                if Xhp(i,1) <= X(j,4) && Xhp(i,1)>0
                    PL(i,1)=Qhpbui(i,1)/(Qhpsh(i,1)*Xhp02(i,1));                          %情况1，存在部分符合下的热泵
                    COPnom(i,1)=Qhpsh(i,1)/Php(i,1);
                    COPpl(i,1)=COPnom(i,1)*(-0.00006*PL(i,1)^6+0.0017*PL(i,1)^5-0.0181*PL(i,1)^4+0.096*PL(i,1)^3-0.2697*PL(i,1)^2+0.4276*PL(i,1)+0.7626);
                    Phppl(i,1)=Qhpbui(i,1)/COPpl(i,1);
                    Qhpscpl(i,1)=Qhpbui(i,1)-Phppl(i,1);
                elseif Xhp(i,1)> X(j,4)
                    Xhp02(i,1)=fix(X(j,4));
                    Qdeloss(i,1)=Qhpbui(i,1)-Qhpsh(i,1)*Xhp02(i,1);                          %Qdeloss代表缺失的热量
                    Xloss(i,1)=Qdeloss(i,1)/Qhpbui(i,1);
                    Phppl(i,1)=Php(i,1)*Xhp02(i,1);
                    Qhpscpl(i,1)=Xhp02(i,1)*Qhpsc(i,1);
                else
                    Qdeloss(i,1)=Qhpbui(i,1);                          %Qdeloss代表缺失的热量
                    Xloss(i,1)=1;
                    Phppl(i,1)=0;
                    Qhpscpl(i,1)=0;
                end                                                                      %Xhp表示最终的热泵数量，带小数点，情况2 X(j,4)表示数量，即边界值
                
                %计算地埋管数量
                Xge01(i,1)=Qhpscpl(i,1)/(6*0.95);
                Xge(i,1)=fix(Xge(i,1));
                
                if Xge01(i,1)-Xge(i,1) > 0                                               %对地埋管数量取整，向上+1
                    Xge(i,1)=fix(Xge01(i,1))+1;
                else
                    Xge(i,1)=Xge01(i,1);
                end
                
            end
            %%--热泵计算结束--%%
            
            %%--电负荷计算--%%
            Pbss(i,1)=Fuhe(i,1)+Phppl(i,1)-Ppv(i,1)*X(j,1)-Pwt(i,1)*X(j,2);
            
            %蓄电池充放电状态判断
            if Pbss(i,1) < 0                                                            %pbss为负，代表要充电，此时的pbss代表不考虑上下限的功率
                SOC(i+1,1)=(1-0.03)*SOC(i,1)-Pbss(i,1)*0.9/(X(j,3)*Eb);
                More(i,1)=abs(Pbss(i,1));                                               %记录多余的电量
                elecnum(i,1)=0;
            elseif Pbss(i,1) > 0                                                        %pbss为正，代表要放电
                SOC(i+1,1)=(1-0.03)*SOC(i,1)-Pbss(i,1)/(X(j,3)*0.9*Eb);
                More(i,1)=0;                                                            %记录多余的浪费电量
                elecnum(i,1)=abs(Pbss(i,1));
            else
                SOC(i+1,1)=SOC(i,1);
                More(i,1)=0;                                                            %记录多余的浪费电量
                elecnum(i,1)=0;
            end
            
            %多余电量判断
            if SOC(i+1,1)>=SOCmax
                SOC(i+1,1)=SOCmax;
                %先计算即将溢出时的蓄电池真实充电功率
                if SOC(i,1)~=SOCmax                                                         %~=是不等于的意思
                    Pbss(i,1)=(SOC(i+1,1)-(1-0.03)*SOC(i,1))*Eb*X(j,3)/(-0.9);
                else
                    Pbss(i,1)=0;
                end
                %在100%的情况下为弃光与弃风
                if More(i,1)>=abs(Pbss(i,1))
                    More(i,1)=More(i,1)-abs(Pbss(i,1));
                else
                    More(i,1)=0;
                end
                %记录多余电量
            else
                More(i,1)=0;
            end
            
            %缺电判断
            if SOC(i+1,1)<=SOCmin
                SOC(i+1,1)=SOCmin;
                %计算即将见底时的蓄电池真实放电功率
                if SOC(i,1)~=SOCmin
                    Pbss(i,1)=(SOC(i+1,1)-(1-0.03)*SOC(i,1))*Eb*(-0.9)*X(j,3);
                else
                    Pbss(i,1)=0;
                end
                %真实缺电率
                if elecnum(i,1)>=abs(Pbss(i,1))
                    elecnum(i,1)=elecnum(i,1)-abs(Pbss(i,1));
                else
                    elecnum(i,1)=0;
                end
            else
                elecnum(i,1)=0;
            end
            
            %外购电耗费
            Celebuy(i,1)=elecnum(i,1)*DianJia(i,1);
            
            %售电盈利
            Celesell(i,1)=More(i,1)*Cse;
            %Emore(i,1)=More(i,1)/(Ppv(i,1)*X(j,1)+Pwt(i,1)*X(j,2)+0.001);
            %%--电负荷计算结束--%%
            
        end
        %%%---16.00-18.00计算结束---%%%
        
        
        %%%---18.00-23.00，非上班时间，非峰谷电价---%%%
        if d(i,1)>0.78 && d(i,1)<0.97
            Ttes(i+1,1)=Ktes*Ates*(Temperature(i,1)-Ttes(i,1))/(mtes*Cw)+Ttes(i,1);
            
            %%%---电负荷计算---%%%
            Pbss(i,1)=Fuhe(i,1)-Ppv(i,1)*X(j,1)-Pwt(i,1)*X(j,2);
            
            %---蓄电池充放电状态判断---%
            if Pbss(i,1) < 0                                                        %pbss为负，代表要充电，此时的pbss代表不考虑上下限的功率
                SOC(i+1,1)=(1-0.03)*SOC(i,1)-Pbss(i,1)*0.9/(X(j,3)*Eb);
                More(i,1)=abs(Pbss(i,1));                                           %记录多余的电量
                elecnum(i,1)=0;
            elseif   Pbss(i,1) > 0                                                  %pbss为正，代表要放电
                SOC(i+1,1)=(1-0.03)*SOC(i,1)-Pbss(i,1)/(X(j,3)*0.9*Eb);
                More(i,1)=0;                                                        %记录多余的浪费电量
                elecnum(i,1)=abs(Pbss(i,1));
            else
                SOC(i+1,1)=SOC(i,1);
                More(i,1)=0;                                                        %记录多余的浪费电量
                elecnum(i,1)=0;
            end
            %---蓄电池充放电状态判断结束---%
            
            %---多余电量判断---%
            if  SOC(i+1,1)>=SOCmax
                SOC(i+1,1)=SOCmax;
                %先计算即将溢出时的蓄电池真实充电功率
                if  SOC(i,1)~=SOCmax                                                %~=是不等于的意思
                    Pbss(i,1)=(SOC(i+1,1)-(1-0.03)*SOC(i,1))*Eb*X(j,3)/(-0.9);
                else
                    Pbss(i,1)=0;
                end
                %记录多余电量
                %在100%的情况下为弃光与弃风
                if More(i,1)>=abs(Pbss(i,1))
                    More(i,1)=More(i,1)-abs(Pbss(i,1));                             %这个设置的很巧妙
                else
                    More(i,1)=0;
                end
                %记录多余电量
            else
                More(i,1)=0;
            end
            %---多余电量判断结束---%
            
            %---缺电判断---%
            if  SOC(i+1,1)<=SOCmin
                SOC(i+1,1)=SOCmin;
                %计算即将见底时的蓄电池真实放电功率
                if  SOC(i,1)~=SOCmin
                    Pbss(i,1)=(SOC(i+1,1)-(1-0.03)*SOC(i,1))*Eb*(-0.9)*X(j,3);
                else
                    Pbss(i,1)=0;
                end
                %在见底时的真实缺电率
                if elecnum(i,1)>=abs(Pbss(i,1))
                    elecnum(i,1)=elecnum(i,1)-abs(Pbss(i,1));                       %这个设置的很巧妙
                else
                    elecnum(i,1)=0;
                end
            else
                elecnum(i,1)=0;
            end
            
            %外购电耗费
            Celebuy(i,1)=0;
            
            %售电盈利
            Celesell(i,1)=More(i,1)*Cse;
            %Emore(i,1)=More(i,1)/(Ppv(i,1)*X(j,1)+Pwt(i,1)*X(j,2)+0.001);
            
            %---缺电判断结束---%
            %%%---电模块结束---%%%
        end
        %%%---18.00-23.00计算结束---%%%
        
        %%%---此季节计算结束---%%%
    end
    
    %% 过渡季，2521-3624h
    for i=2521:3624
        
        %%%---电负荷计算---%%%
        Pbss(i,1)=Fuhe(i,1)-Ppv(i,1)*X(j,1)-Pwt(i,1)*X(j,2);
        
        %---蓄电池充放电状态判断---%
        if Pbss(i,1) < 0                                                        %pbss为负，代表要充电，此时的pbss代表不考虑上下限的功率
            SOC(i+1,1)=(1-0.03)*SOC(i,1)-Pbss(i,1)*0.9/(X(j,3)*Eb);
            More(i,1)=abs(Pbss(i,1));                                           %记录多余的电量
            elecnum(i,1)=0;
        elseif   Pbss(i,1) > 0                                                  %pbss为正，代表要放电
            SOC(i+1,1)=(1-0.03)*SOC(i,1)-Pbss(i,1)/(X(j,3)*0.9*Eb);
            More(i,1)=0;                                                        %记录多余的浪费电量
            elecnum(i,1)=abs(Pbss(i,1));
        else
            SOC(i+1,1)=SOC(i,1);
            More(i,1)=0;                                                        %记录多余的浪费电量
            elecnum(i,1)=0;
        end
        %---蓄电池充放电状态判断结束---%
        
        %---多余电量判断---%
        if  SOC(i+1,1)>=SOCmax
            SOC(i+1,1)=SOCmax;
            %先计算即将溢出时的蓄电池真实充电功率
            if  SOC(i,1)~=SOCmax                                                %~=是不等于的意思
                Pbss(i,1)=(SOC(i+1,1)-(1-0.03)*SOC(i,1))*Eb*X(j,3)/(-0.9);
            else
                Pbss(i,1)=0;
            end
            %记录多余电量
            %在100%的情况下为弃光与弃风
            if More(i,1)>=abs(Pbss(i,1))
                More(i,1)=More(i,1)-abs(Pbss(i,1));                             %这个设置的很巧妙
            else
                More(i,1)=0;
            end
            %记录多余电量
        else
            More(i,1)=0;
        end
        %---多余电量判断结束---%
        
        %---缺电判断---%
        if  SOC(i+1,1)<=SOCmin
            SOC(i+1,1)=SOCmin;
            %计算即将见底时的蓄电池真实放电功率
            if  SOC(i,1)~=SOCmin
                Pbss(i,1)=(SOC(i+1,1)-(1-0.03)*SOC(i,1))*Eb*(-0.9)*X(j,3);
            else
                Pbss(i,1)=0;
            end
            %在见底时的真实缺电率
            if elecnum(i,1)>=abs(Pbss(i,1))
                elecnum(i,1)=elecnum(i,1)-abs(Pbss(i,1));                       %这个设置的很巧妙
            else
                elecnum(i,1)=0;
            end
        else
            elecnum(i,1)=0;
        end
        
        %外购电耗费
        Celebuy(i,1)=elecnum(i,1)*DianJia(i,1);
        
        %售电盈利
        Celesell(i,1)=More(i,1)*Cse;
        %Emore(i,1)=More(i,1)/(Ppv(i,1)*X(j,1)+Pwt(i,1)*X(j,2)+0.001);
        
        %---缺电判断结束---%
        %%%---电模块结束---%%%
        
    end
    
    %% 制冷季，3625h-5832h
    for i=3625:5832
        
        a(i,1)=i;
        b(i,1)=a(i,1)/24;
        c(i,1)=fix(b(i,1));
        d(i,1)=b(i,1)-c(i,1);
        
        %%%---23.00-07.00，峰谷时间，热泵向水箱制冷---%%%
        if d(i,1)<0.31
            
            if Ttes(i,1)>Ttessdcc                                                                                       %水箱温度控制启停判断
                Qhptes(i,1)=abs((Ttessdc-Ttes(i,1))*Cw*X(j,5)*mtes+X(j,5)*Ktes*Ates*(Ttes(i,1)-Temperature(i,1)))/3600;      %热泵峰谷制冷负荷
                
                %%---热泵部分负荷及缺失计算---%%
                
                %第一次计算热泵数量
                Xhp01(i,1)=Qhptes(i,1)/Qhpsc(i,1);
                
                if Xhp01(i,1)<1
                    Xhp02(i,1)=0;
                else
                    Xhp02(i,1)=fix(Xhp01(i,1))+1;
                end
                
                Xhp(i,1)=Xhp02(i,1);
                
                %第二次计算热泵数量，根据边界修正,并计算缺失量与部分负荷下的热泵功耗、COP
                if Xhp(i,1)<=X(j,4) && Xhp(i,1)>0
                    PL(i,1)=Qhptes(i,1)/(Qhpsc(i,1)*Xhp02(i,1));                          %情况1，存在部分符合下的热泵
                    COPnom(i,1)=Qhpsc(i,1)/Php(i,1);
                    COPpl(i,1)=COPnom(i,1)*(-0.00006*PL(i,1)^6+0.0017*PL(i,1)^5-0.0181*PL(i,1)^4+0.096*PL(i,1)^3-0.2697*PL(i,1)^2+0.4276*PL(i,1)+0.7626);
                    Phppl(i,1)=Qhptes(i,1)/COPpl(i,1);                                    %计算部分负荷下该热泵的功耗
                    Qhpshpl(i,1)=Qhptes(i,1)+Phppl(i,1);
                    Ttes(i+1,1)=Ttessdc;
                elseif Xhp(i,1)>X(j,4)
                    Xhp02(i,1)=fix(X(j,4));
                    %Qdeloss(i,1)=Qhptes(i,1)-Qhpsh(i,1)*Xhp02(i,1);                          %Qdeloss代表缺失的热量
                    %Xloss(i,1)=Qdeloss(i,1)/Qhptes(i,1);
                    Phppl(i,1)=Php(i,1)*Xhp02(i,1);
                    Qhpshpl(i,1)=Xhp02(i,1)*Qhpsh(i,1);
                    Ttes(i+1,1)=(-Xhp02(i,1)*Qhpsc(i,1)*3600-X(j,5)*Ktes*Ates*(Ttes(i,1)-Temperature(i,1)))/(mtes*Cw*X(j,5))+Ttes(i,1);
                else
                    Phppl(i,1)=0;
                    Qhpshpl(i,1)=0;
                    Ttes(i+1,1)=Ktes*Ates*(Temperature(i,1)-Ttes(i,1))/(mtes*Cw)+Ttes(i,1);
                end                                                                      %Xhp表示最终的热泵数量，带小数点，情况2 X(j,4)表示数量，即边界值
                
                Phptes(i,1)=Phppl(i,1);
                
                %计算地埋管数量
                Xge01(i,1)=Qhpshpl(i,1)/(6*0.95);
                Xge(i,1)=fix(Xge01(i,1));
                if Xge01(i,1)-Xge(i,1) > 0                                                %对地埋管数量取整，向上+1
                    Xge(i,1)=fix(Xge01(i,1))+1;
                else
                    Xge(i,1)=Xge01(i,1);
                end
                
            else
                Ttes(i+1,1)=Ktes*Ates*(Temperature(i,1)-Ttes(i,1))/(mtes*Cw)+Ttes(i,1);
            end
            
            %%---热泵部分负荷及缺失计算结束---%%
            
            %%%---电负荷计算---%%%
            Pbss(i,1)=Fuhe(i,1)-Ppv(i,1)*X(j,1)-Pwt(i,1)*X(j,2);
            
            %---蓄电池充放电状态判断---%
            if Pbss(i,1) < 0
                SOC(i+1,1)=(1-0.03)*SOC(i,1)-Pbss(i,1)*0.9/(X(j,3)*Eb);
                More(i,1)=abs(Pbss(i,1));
                elecnum(i,1)=0;
            elseif   Pbss(i,1) > 0
                SOC(i+1,1)=(1-0.03)*SOC(i,1)-Pbss(i,1)/(X(j,3)*0.9*Eb);
                More(i,1)=0;
                elecnum(i,1)=abs(Pbss(i,1));
            else
                SOC(i+1,1)=SOC(i,1);
                More(i,1)=0;
                elecnum(i,1)=0;
            end
            %---蓄电池充放电状态判断结束---%
            
            %---多余电量判断---%
            if  SOC(i+1,1)>=SOCmax
                SOC(i+1,1)=SOCmax;
                %先计算即将溢出时的蓄电池真实充电功率
                if  SOC(i,1)~=SOCmax
                    Pbss(i,1)=(SOC(i+1,1)-(1-0.03)*SOC(i,1))*Eb*X(j,3)/(-0.9);
                else
                    Pbss(i,1)=0;
                end
                %记录多余电量
                %在100%的情况下为弃光与弃风
                if More(i,1)>=abs(Pbss(i,1))
                    More(i,1)=More(i,1)-abs(Pbss(i,1));
                else
                    More(i,1)=0;
                end
                %记录多余电量
            else
                More(i,1)=0;
            end
            %---多余电量判断结束---%
            
            %---缺电判断---%
            if  SOC(i+1,1)<=SOCmin
                SOC(i+1,1)=SOCmin;
                %计算即将见底时的蓄电池真实放电功率
                if  SOC(i,1)~=SOCmin
                    Pbss(i,1)=(SOC(i+1,1)-(1-0.03)*SOC(i,1))*Eb*(-0.9)*X(j,3);
                else
                    Pbss(i,1)=0;
                end
                %在见底时的真实缺电率
                if elecnum(i,1)>=abs(Pbss(i,1))
                    elecnum(i,1)=elecnum(i,1)-abs(Pbss(i,1));
                else
                    elecnum(i,1)=0;
                end
            else
                elecnum(i,1)=0;
            end
            
            %外购电耗费
            Celebuy(i,1)=Phptes(i,1)*DianJia(i,1);
            
            %售电盈利
            Celesell(i,1)=More(i,1)*Cse;
            %Emore(i,1)=More(i,1)/(Ppv(i,1)*X(j,1)+Pwt(i,1)*X(j,2)+0.001);
            
        end
        
        %%%---23.00-07.00计算结束---%%%
        
        %%%---07.00-08.00，非工作时间，非峰谷电价---%%%
        if d(i,1)>0.31 && d(i,1)<0.36
            Ttes(i+1,1)=Ktes*Ates*(Temperature(i,1)-Ttes(i,1))/(mtes*Cw)+Ttes(i,1);
            
            %%%---电负荷计算---%%%
            Pbss(i,1)=Fuhe(i,1)-Ppv(i,1)*X(j,1)-Pwt(i,1)*X(j,2);
            
            %---蓄电池充放电状态判断---%
            if Pbss(i,1) < 0
                SOC(i+1,1)=(1-0.03)*SOC(i,1)-Pbss(i,1)*0.9/(X(j,3)*Eb);
                More(i,1)=abs(Pbss(i,1));
                elecnum(i,1)=0;
            elseif   Pbss(i,1) > 0
                SOC(i+1,1)=(1-0.03)*SOC(i,1)-Pbss(i,1)/(X(j,3)*0.9*Eb);
                More(i,1)=0;
                elecnum(i,1)=abs(Pbss(i,1));
            else
                SOC(i+1,1)=SOC(i,1);
                More(i,1)=0;
                elecnum(i,1)=0;
            end
            %---蓄电池充放电状态判断结束---%
            
            %---多余电量判断---%
            if  SOC(i+1,1)>=SOCmax
                SOC(i+1,1)=SOCmax;
                %先计算即将溢出时的蓄电池真实充电功率
                if  SOC(i,1)~=SOCmax
                    Pbss(i,1)=(SOC(i+1,1)-(1-0.03)*SOC(i,1))*Eb*X(j,3)/(-0.9);
                else
                    Pbss(i,1)=0;
                end
                %记录多余电量
                %在100%的情况下为弃光与弃风
                if More(i,1)>=abs(Pbss(i,1))
                    More(i,1)=More(i,1)-abs(Pbss(i,1));
                else
                    More(i,1)=0;
                end
                %记录多余电量
            else
                More(i,1)=0;
            end
            %---多余电量判断结束---%
            
            %---缺电判断---%
            if  SOC(i+1,1)<=SOCmin
                SOC(i+1,1)=SOCmin;
                %计算即将见底时的蓄电池真实放电功率
                if  SOC(i,1)~=SOCmin
                    Pbss(i,1)=(SOC(i+1,1)-(1-0.03)*SOC(i,1))*Eb*(-0.9)*X(j,3);
                else
                    Pbss(i,1)=0;
                end
                %在见底时的真实缺电率
                if elecnum(i,1)>=abs(Pbss(i,1))
                    elecnum(i,1)=elecnum(i,1)-abs(Pbss(i,1));
                else
                    elecnum(i,1)=0;
                end
            else
                elecnum(i,1)=0;
            end
            
            %外购电耗费
            Celebuy(i,1)=0;
            
            %售电盈利
            Celesell(i,1)=More(i,1)*Cse;
            %Emore(i,1)=More(i,1)/(Ppv(i,1)*X(j,1)+Pwt(i,1)*X(j,2)+0.001);
            
        end
        %%%---07.00-08.00计算结束---%%%
        
        %%%---8.00-9.00,工作时间，平常电价---%%%
        if d(i,1)>0.36 && d(i,1)<0.39
            
            %%--水箱温度计算--%%
            Ttes(i+1,1)=Ktes*Ates*(Temperature(i,1)-Ttes(i,1))/(mtes*Cw)+Ttes(i,1);
            %%--水箱温度计算结束--%%
            
            %%--热泵计算--%%
            
            %第一次计算热泵数量
            Xhp01(i,1)=CoolingFuhe(i,1)/Qhpsc(i,1);
            
            if Xhp01(i,1)<1
                Xhp02(i,1)=0;
            else
                Xhp02(i,1)=fix(Xhp01(i,1))+1;
            end
            
            Xhp(i,1)=Xhp02(i,1);
                
            %第二次计算热泵数量，根据边界修正,并计算缺失量与部分负荷下的热泵功耗、COP
            if Xhp(i,1)<=X(j,4) && Xhp(i,1)>0
                PL(i,1)=CoolingFuhe(i,1)/(Qhpsc(i,1)*Xhp02(i,1));                          %情况1，存在部分符合下的热泵
                COPnom(i,1)=Qhpsc(i,1)/Php(i,1);
                COPpl(i,1)=COPnom(i,1)*(-0.00006*PL(i,1)^6+0.0017*PL(i,1)^5-0.0181*PL(i,1)^4+0.096*PL(i,1)^3-0.2697*PL(i,1)^2+0.4276*PL(i,1)+0.7626);
                Phppl(i,1)=CoolingFuhe(i,1)/COPpl(i,1);
                Qhpshpl(i,1)=CoolingFuhe(i,1)+Phppl(i,1);
            elseif Xhp(i,1)>X(j,4)
                Xhp02(i,1)=fix(X(j,4));
                Qdeloss(i,1)=CoolingFuhe(i,1)-Qhpsc(i,1)*Xhp02(i,1);                          %Qdeloss代表缺失的热量
                Xloss(i,1)=Qdeloss(i,1)/CoolingFuhe(i,1);
                Phppl(i,1)=Php(i,1)*Xhp02(i,1);
                Qhpshpl(i,1)=Xhp02(i,1)*Qhpsh(i,1);
            else
                Qdeloss(i,1)=CoolingFuhe(i,1);                          %Qdeloss代表缺失的热量
                Xloss(i,1)=1;
                Phppl(i,1)=0;
                Qhpshpl(i,1)=0;
            end                                                                      %Xhp表示最终的热泵数量，带小数点，情况2 X(j,4)表示数量，即边界值
            
            %计算地埋管数量
            Xge01(i,1)=Qhpshpl(i,1)/(6*0.95);
            Xge(i,1)=fix(Xge01(i,1));
            
            if Xge01(i,1)-Xge(i,1) > 0                                                %对地埋管数量取整，向上+1
                Xge(i,1)=fix(Xge01(i,1))+1;
            else
                Xge(i,1)=Xge01(i,1);
            end
            %%---热泵计算结束---%%
            
            %%--电负荷计算--%%
            Pbss(i,1)=Fuhe(i,1)+Phppl(i,1)-Ppv(i,1)*X(j,1)-Pwt(i,1)*X(j,2);
            
            %蓄电池充放电状态判断
            if Pbss(i,1) < 0                                                         %pbss为负，代表要充电，此时的pbss代表不考虑上下限的功率
                SOC(i+1,1)=(1-0.03)*SOC(i,1)-Pbss(i,1)*0.9/(X(j,3)*Eb);
                More(i,1)=abs(Pbss(i,1));                                            %记录多余的电量
                elecnum(i,1)=0;
            elseif Pbss(i,1) > 0                                                     %pbss为正，代表要放电
                SOC(i+1,1)=(1-0.03)*SOC(i,1)-Pbss(i,1)/(X(j,3)*0.9*Eb);
                More(i,1)=0;                                                         %记录多余的浪费电量
                elecnum(i,1)=abs(Pbss(i,1));
            else
                SOC(i+1,1)=SOC(i,1);
                More(i,1)=0;                                                         %记录多余的浪费电量
                elecnum(i,1)=0;
            end
            
            %多余电量判断
            if SOC(i+1,1)>=SOCmax
                SOC(i+1,1)=SOCmax;
                %先计算即将溢出时的蓄电池真实充电功率
                if SOC(i,1)~=SOCmax                                                      %~=是不等于的意思
                    Pbss(i,1)=(SOC(i+1,1)-(1-0.03)*SOC(i,1))*Eb*X(j,3)/(-0.9);
                else
                    Pbss(i,1)=0;
                end
                %在100%的情况下为弃光与弃风
                if More(i,1)>=abs(Pbss(i,1))
                    More(i,1)=More(i,1)-abs(Pbss(i,1));
                else
                    More(i,1)=0;
                end
                %记录多余电量
            else
                More(i,1)=0;
            end
            
            %缺电判断
            if SOC(i+1,1)<=SOCmin
                SOC(i+1,1)=SOCmin;
                %计算即将见底时的蓄电池真实放电功率
                if SOC(i,1)~=SOCmin
                    Pbss(i,1)=(SOC(i+1,1)-(1-0.03)*SOC(i,1))*Eb*(-0.9)*X(j,3);
                else
                    Pbss(i,1)=0;
                end
                %真实缺电率
                if elecnum(i,1)>=abs(Pbss(i,1))
                    elecnum(i,1)=elecnum(i,1)-abs(Pbss(i,1));
                else
                    elecnum(i,1)=0;
                end
            else
                elecnum(i,1)=0;
            end
            
            %外购电耗费
            Celebuy(i,1)=elecnum(i,1)*DianJia(i,1);
            
            %售电盈利
            Celesell(i,1)=More(i,1)*Cse;
            %Emore(i,1)=More(i,1)/(Ppv(i,1)*X(j,1)+Pwt(i,1)*X(j,2)+0.001);
            %%--电负荷计算结束--%%
            
        end
        %%%---8.00-9.00计算结束---%%%
        
        %%%---9.00-12.00,工作时间，峰值电价，开启水箱---%%%
        if d(i,1)>0.39 && d(i,1)<0.51
            
            %%--热泵计算--%%
            Qtes(i,1)=((Tbuc-Ttes(i,1))*Cw*X(j,5)*mtes+X(j,5)*Ktes*Ates*(Ttes(i,1)-Temperature(i,1)))/3600;
            
            if Qtes(i,1)<=0
                Qtes(i,1)=0;
                
                %%--水箱温度计算--%%
                Ttes(i+1,1)=Ktes*Ates*(Temperature(i,1)-Ttes(i,1))/(mtes*Cw)+Ttes(i,1);
                %%--水箱温度计算结束--%%
            
                %%--热泵计算--%%
            
                %第一次计算热泵数量
                Xhp01(i,1)=CoolingFuhe(i,1)/Qhpsc(i,1);
                
                if Xhp01(i,1)<1
                    Xhp02(i,1)=0;
                else
                    Xhp02(i,1)=fix(Xhp01(i,1))+1;
                end
                
                Xhp(i,1)=Xhp02(i,1);
                
                %第二次计算热泵数量，根据边界修正,并计算缺失量与部分负荷下的热泵功耗、COP
                if Xhp(i,1)<=X(j,4) && Xhp(i,1)>0
                    PL(i,1)=CoolingFuhe(i,1)/(Qhpsc(i,1)*Xhp02(i,1));                          %情况1，存在部分符合下的热泵
                    COPnom(i,1)=Qhpsc(i,1)/Php(i,1);
                    COPpl(i,1)=COPnom(i,1)*(-0.00006*PL(i,1)^6+0.0017*PL(i,1)^5-0.0181*PL(i,1)^4+0.096*PL(i,1)^3-0.2697*PL(i,1)^2+0.4276*PL(i,1)+0.7626);
                    Phppl(i,1)=CoolingFuhe(i,1)/COPpl(i,1);
                    Qhpshpl(i,1)=CoolingFuhe(i,1)+Phppl(i,1);
                elseif Xhp(i,1)>X(j,4)
                    Xhp02(i,1)=fix(X(j,4));
                    Qdeloss(i,1)=CoolingFuhe(i,1)-Qhpsc(i,1)*Xhp02(i,1);                          %Qdeloss代表缺失的热量
                    Xloss(i,1)=Qdeloss(i,1)/CoolingFuhe(i,1);
                    Phppl(i,1)=Php(i,1)*Xhp02(i,1);
                    Qhpshpl(i,1)=Xhp02(i,1)*Qhpsh(i,1);
                else
                    Qdeloss(i,1)=CoolingFuhe(i,1);                          %Qdeloss代表缺失的热量
                    Xloss(i,1)=1;
                    Phppl(i,1)=0;
                    Qhpshpl(i,1)=0;
                end                                                                      %Xhp表示最终的热泵数量，带小数点，情况2 X(j,4)表示数量，即边界值
            
                %计算地埋管数量
                Xge01(i,1)=Qhpshpl(i,1)/(6*0.95);
                Xge(i,1)=fix(Xge01(i,1));
            
                if Xge01(i,1)-Xge(i,1) > 0                                                %对地埋管数量取整，向上+1
                    Xge(i,1)=fix(Xge01(i,1))+1;
                else
                    Xge(i,1)=Xge01(i,1);
                end
                %%---热泵计算结束---%%

            end
            
            %情况1/2/3
            if Qtes(i,1)>=CoolingFuhe(i,1)
                Ttes(i+1,1)=(X(j,5)*Ktes*Ates*(Temperature(i,1)-Ttes(i,1))+CoolingFuhe(i,1)*3600)/(X(j,5)*mtes*Cw)+Ttes(i,1);
                Xhp02(i,1)=0;
            elseif Qtes(i,1)<CoolingFuhe(i,1) && Qtes(i,1)>0
                Qhpbui(i,1)=CoolingFuhe(i,1)-Qtes(i,1);
                Ttes(i+1,1)=Tbuc;
                
                %情况2/3
                %第一次计算热泵数量
                Xhp01(i,1)=Qhpbui(i,1)/Qhpsc(i,1);
                
                if Xhp01(i,1)<1
                    Xhp02(i,1)=0;
                else
                    Xhp02(i,1)=fix(Xhp01(i,1))+1;
                end
                
                Xhp(i,1)=Xhp02(i,1);
                
                %第二次计算热泵数量，根据边界修正,并计算缺失量与部分负荷下的热泵功耗、COP
                if Xhp(i,1)<=X(j,4) && Xhp(i,1)>0
                    PL(i,1)=Qhpbui(i,1)/(Qhpsc(i,1)*Xhp02(i,1));                          %情况1，存在部分符合下的热泵
                    COPnom(i,1)=Qhpsc(i,1)/Php(i,1);
                    COPpl(i,1)=COPnom(i,1)*(-0.00006*PL(i,1)^6+0.0017*PL(i,1)^5-0.0181*PL(i,1)^4+0.096*PL(i,1)^3-0.2697*PL(i,1)^2+0.4276*PL(i,1)+0.7626);
                    Phppl(i,1)=Qhpbui(i,1)/COPpl(i,1);
                    Qhpshpl(i,1)=Qhpbui(i,1)+Phppl(i,1);
                elseif Xhp(i,1)>X(j,4)
                    Xhp02(i,1)=fix(X(j,4));
                    Qdeloss(i,1)=Qhpbui(i,1)-Qhpsc(i,1)*Xhp02(i,1);                          %Qdeloss代表缺失的热量
                    Xloss(i,1)=Qdeloss(i,1)/Qhpbui(i,1);
                    Phppl(i,1)=Php(i,1)*Xhp02(i,1);
                    Qhpshpl(i,1)=Xhp02(i,1)*Qhpsh(i,1);
                else
                    Qdeloss(i,1)=Qhpbui(i,1);                          %Qdeloss代表缺失的热量
                    Xloss(i,1)=1;
                    Phppl(i,1)=0;
                    Qhpshpl(i,1)=0;
                end                                                                      %Xhp表示最终的热泵数量，带小数点，情况2 X(j,4)表示数量，即边界值
                
                %计算地埋管数量
                Xge01(i,1)=Qhpshpl(i,1)/(6*0.95);
                Xge(i,1)=fix(Xge(i,1));
                
                if Xge01(i,1)-Xge(i,1) > 0                                               %对地埋管数量取整，向上+1
                    Xge(i,1)=fix(Xge01(i,1))+1;
                else
                    Xge(i,1)=Xge01(i,1);
                end
                
            end
            %%--热泵计算结束--%%
            
            %%--电负荷计算--%%
            Pbss(i,1)=Fuhe(i,1)+Phppl(i,1)-Ppv(i,1)*X(j,1)-Pwt(i,1)*X(j,2);
            
            %蓄电池充放电状态判断
            if Pbss(i,1) < 0                                                            %pbss为负，代表要充电，此时的pbss代表不考虑上下限的功率
                SOC(i+1,1)=(1-0.03)*SOC(i,1)-Pbss(i,1)*0.9/(X(j,3)*Eb);
                More(i,1)=abs(Pbss(i,1));                                               %记录多余的电量
                elecnum(i,1)=0;
            elseif Pbss(i,1) > 0                                                        %pbss为正，代表要放电
                SOC(i+1,1)=(1-0.03)*SOC(i,1)-Pbss(i,1)/(X(j,3)*0.9*Eb);
                More(i,1)=0;                                                            %记录多余的浪费电量
                elecnum(i,1)=abs(Pbss(i,1));
            else
                SOC(i+1,1)=SOC(i,1);
                More(i,1)=0;                                                            %记录多余的浪费电量
                elecnum(i,1)=0;
            end
            
            %多余电量判断
            if SOC(i+1,1)>=SOCmax
                SOC(i+1,1)=SOCmax;
                %先计算即将溢出时的蓄电池真实充电功率
                if SOC(i,1)~=SOCmax                                                         %~=是不等于的意思
                    Pbss(i,1)=(SOC(i+1,1)-(1-0.03)*SOC(i,1))*Eb*X(j,3)/(-0.9);
                else
                    Pbss(i,1)=0;
                end
                %在100%的情况下为弃光与弃风
                if More(i,1)>=abs(Pbss(i,1))
                    More(i,1)=More(i,1)-abs(Pbss(i,1));
                else
                    More(i,1)=0;
                end
                %记录多余电量
            else
                More(i,1)=0;
            end
            
            %缺电判断
            if SOC(i+1,1)<=SOCmin
                SOC(i+1,1)=SOCmin;
                %计算即将见底时的蓄电池真实放电功率
                if SOC(i,1)~=SOCmin
                    Pbss(i,1)=(SOC(i+1,1)-(1-0.03)*SOC(i,1))*Eb*(-0.9)*X(j,3);
                else
                    Pbss(i,1)=0;
                end
                %真实缺电率
                if elecnum(i,1)>=abs(Pbss(i,1))
                    elecnum(i,1)=elecnum(i,1)-abs(Pbss(i,1));
                else
                    elecnum(i,1)=0;
                end
            else
                elecnum(i,1)=0;
            end
            
            %外购电耗费
            Celebuy(i,1)=elecnum(i,1)*DianJia(i,1);
            
            %售电盈利
            Celesell(i,1)=More(i,1)*Cse;
            %Emore(i,1)=More(i,1)/(Ppv(i,1)*X(j,1)+Pwt(i,1)*X(j,2)+0.001);
            %%--电负荷计算结束--%%
            
        end
        %%%---9.00-12.00计算结束---%%%
        
        %%%---12.00-16.00,工作时间，平常电价---%%%
        if d(i,1)>0.51 && d(i,1)<0.69
            
            %%--水箱温度计算--%%
            Ttes(i+1,1)=Ktes*Ates*(Temperature(i,1)-Ttes(i,1))/(mtes*Cw)+Ttes(i,1);
            %%--水箱温度计算结束--%%
            
            %%--热泵计算--%%
            
            %第一次计算热泵数量
            Xhp01(i,1)=CoolingFuhe(i,1)/Qhpsc(i,1);
            
            if Xhp01(i,1)<1
                Xhp02(i,1)=0;
            else
                Xhp02(i,1)=fix(Xhp01(i,1))+1;
            end
            
            Xhp(i,1)=Xhp02(i,1);
                
            %第二次计算热泵数量，根据边界修正,并计算缺失量与部分负荷下的热泵功耗、COP
            if Xhp(i,1)<=X(j,4) && Xhp(i,1)>0
                PL(i,1)=CoolingFuhe(i,1)/(Qhpsc(i,1)*Xhp02(i,1));                          %情况1，存在部分符合下的热泵
                COPnom(i,1)=Qhpsc(i,1)/Php(i,1);
                COPpl(i,1)=COPnom(i,1)*(-0.00006*PL(i,1)^6+0.0017*PL(i,1)^5-0.0181*PL(i,1)^4+0.096*PL(i,1)^3-0.2697*PL(i,1)^2+0.4276*PL(i,1)+0.7626);
                Phppl(i,1)=CoolingFuhe(i,1)/COPpl(i,1);
                Qhpshpl(i,1)=CoolingFuhe(i,1)+Phppl(i,1);
            elseif Xhp(i,1)>X(j,4)
                Xhp02(i,1)=fix(X(j,4));
                Qdeloss(i,1)=CoolingFuhe(i,1)-Qhpsc(i,1)*Xhp02(i,1);                          %Qdeloss代表缺失的热量
                Xloss(i,1)=Qdeloss(i,1)/CoolingFuhe(i,1);
                Phppl(i,1)=Php(i,1)*Xhp02(i,1);
                Qhpshpl(i,1)=Xhp02(i,1)*Qhpsh(i,1);
            else
                Qdeloss(i,1)=CoolingFuhe(i,1);                          %Qdeloss代表缺失的热量
                Xloss(i,1)=1;
                Phppl(i,1)=0;
                Qhpshpl(i,1)=0;
            end                                                                      %Xhp表示最终的热泵数量，带小数点，情况2 X(j,4)表示数量，即边界值
            
            %计算地埋管数量
            Xge01(i,1)=Qhpshpl(i,1)/(6*0.95);
            Xge(i,1)=fix(Xge01(i,1));
            
            if Xge01(i,1)-Xge(i,1) > 0                                                %对地埋管数量取整，向上+1
                Xge(i,1)=fix(Xge01(i,1))+1;
            else
                Xge(i,1)=Xge01(i,1);
            end
            %%---热泵计算结束---%%
            
            %%--电负荷计算--%%
            Pbss(i,1)=Fuhe(i,1)+Phppl(i,1)-Ppv(i,1)*X(j,1)-Pwt(i,1)*X(j,2);
            
            %蓄电池充放电状态判断
            if Pbss(i,1) < 0                                                         %pbss为负，代表要充电，此时的pbss代表不考虑上下限的功率
                SOC(i+1,1)=(1-0.03)*SOC(i,1)-Pbss(i,1)*0.9/(X(j,3)*Eb);
                More(i,1)=abs(Pbss(i,1));                                            %记录多余的电量
                elecnum(i,1)=0;
            elseif Pbss(i,1) > 0                                                     %pbss为正，代表要放电
                SOC(i+1,1)=(1-0.03)*SOC(i,1)-Pbss(i,1)/(X(j,3)*0.9*Eb);
                More(i,1)=0;                                                         %记录多余的浪费电量
                elecnum(i,1)=abs(Pbss(i,1));
            else
                SOC(i+1,1)=SOC(i,1);
                More(i,1)=0;                                                         %记录多余的浪费电量
                elecnum(i,1)=0;
            end
            
            %多余电量判断
            if SOC(i+1,1)>=SOCmax
                SOC(i+1,1)=SOCmax;
                %先计算即将溢出时的蓄电池真实充电功率
                if SOC(i,1)~=SOCmax                                                      %~=是不等于的意思
                    Pbss(i,1)=(SOC(i+1,1)-(1-0.03)*SOC(i,1))*Eb*X(j,3)/(-0.9);
                else
                    Pbss(i,1)=0;
                end
                %在100%的情况下为弃光与弃风
                if More(i,1)>=abs(Pbss(i,1))
                    More(i,1)=More(i,1)-abs(Pbss(i,1));
                else
                    More(i,1)=0;
                end
                %记录多余电量
            else
                More(i,1)=0;
            end
            
            %缺电判断
            if SOC(i+1,1)<=SOCmin
                SOC(i+1,1)=SOCmin;
                %计算即将见底时的蓄电池真实放电功率
                if SOC(i,1)~=SOCmin
                    Pbss(i,1)=(SOC(i+1,1)-(1-0.03)*SOC(i,1))*Eb*(-0.9)*X(j,3);
                else
                    Pbss(i,1)=0;
                end
                %真实缺电率
                if elecnum(i,1)>=abs(Pbss(i,1))
                    elecnum(i,1)=elecnum(i,1)-abs(Pbss(i,1));
                else
                    elecnum(i,1)=0;
                end
            else
                elecnum(i,1)=0;
            end
            
            %外购电耗费
            Celebuy(i,1)=elecnum(i,1)*DianJia(i,1);
            
            %售电盈利
            Celesell(i,1)=More(i,1)*Cse;
            %Emore(i,1)=More(i,1)/(Ppv(i,1)*X(j,1)+Pwt(i,1)*X(j,2)+0.001);
            %%--电负荷计算结束--%%
            
        end
        %%%---12.00-16.00计算结束---%%%
        
        %%%---16.00-18.00,工作时间，峰值电价，开启水箱---%%%
        if d(i,1)>0.69 && d(i,1)<0.78
            
            %%--热泵计算--%%
            Qtes(i,1)=((Tbuc-Ttes(i,1))*Cw*X(j,5)*mtes+X(j,5)*Ktes*Ates*(Ttes(i,1)-Temperature(i,1)))/3600;
            
            if Qtes(i,1)<=0
                Qtes(i,1)=0;
                
                %%--水箱温度计算--%%
                Ttes(i+1,1)=Ktes*Ates*(Temperature(i,1)-Ttes(i,1))/(mtes*Cw)+Ttes(i,1);
                %%--水箱温度计算结束--%%
            
                %%--热泵计算--%%
            
                %第一次计算热泵数量
                Xhp01(i,1)=CoolingFuhe(i,1)/Qhpsc(i,1);
                
                if Xhp01(i,1)<1
                    Xhp02(i,1)=0;
                else
                    Xhp02(i,1)=fix(Xhp01(i,1))+1;
                end
                
                Xhp(i,1)=Xhp02(i,1);
                
                %第二次计算热泵数量，根据边界修正,并计算缺失量与部分负荷下的热泵功耗、COP
                if Xhp(i,1)<=X(j,4) && Xhp(i,1)>0
                    PL(i,1)=CoolingFuhe(i,1)/(Qhpsc(i,1)*Xhp02(i,1));                          %情况1，存在部分符合下的热泵
                    COPnom(i,1)=Qhpsc(i,1)/Php(i,1);
                    COPpl(i,1)=COPnom(i,1)*(-0.00006*PL(i,1)^6+0.0017*PL(i,1)^5-0.0181*PL(i,1)^4+0.096*PL(i,1)^3-0.2697*PL(i,1)^2+0.4276*PL(i,1)+0.7626);
                    Phppl(i,1)=CoolingFuhe(i,1)/COPpl(i,1);
                    Qhpshpl(i,1)=CoolingFuhe(i,1)+Phppl(i,1);
                elseif Xhp(i,1)>X(j,4)
                    Xhp02(i,1)=fix(X(j,4));
                    Qdeloss(i,1)=CoolingFuhe(i,1)-Qhpsc(i,1)*Xhp02(i,1);                          %Qdeloss代表缺失的热量
                    Xloss(i,1)=Qdeloss(i,1)/CoolingFuhe(i,1);
                    Phppl(i,1)=Php(i,1)*Xhp02(i,1);
                    Qhpshpl(i,1)=Xhp02(i,1)*Qhpsh(i,1);
                else
                    Qdeloss(i,1)=CoolingFuhe(i,1);                          %Qdeloss代表缺失的热量
                    Xloss(i,1)=1;
                    Phppl(i,1)=0;
                    Qhpshpl(i,1)=0;
                end                                                                      %Xhp表示最终的热泵数量，带小数点，情况2 X(j,4)表示数量，即边界值
            
                %计算地埋管数量
                Xge01(i,1)=Qhpshpl(i,1)/(6*0.95);
                Xge(i,1)=fix(Xge01(i,1));
            
                if Xge01(i,1)-Xge(i,1) > 0                                                %对地埋管数量取整，向上+1
                    Xge(i,1)=fix(Xge01(i,1))+1;
                else
                    Xge(i,1)=Xge01(i,1);
                end
                %%---热泵计算结束---%%

            end
            
            %情况1/2/3
            if Qtes(i,1)>=CoolingFuhe(i,1)
                Ttes(i+1,1)=(X(j,5)*Ktes*Ates*(Temperature(i,1)-Ttes(i,1))+CoolingFuhe(i,1)*3600)/(X(j,5)*mtes*Cw)+Ttes(i,1);
                Xhp02(i,1)=0;
            elseif Qtes(i,1)<CoolingFuhe(i,1) && Qtes(i,1)>0
                Qhpbui(i,1)=CoolingFuhe(i,1)-Qtes(i,1);
                Ttes(i+1,1)=Tbuc;
                
                %情况2/3
                %第一次计算热泵数量
                Xhp01(i,1)=Qhpbui(i,1)/Qhpsc(i,1);
                
                if Xhp01(i,1)<1
                    Xhp02(i,1)=0;
                else
                    Xhp02(i,1)=fix(Xhp01(i,1))+1;
                end
                
                Xhp(i,1)=Xhp02(i,1);
                
                %第二次计算热泵数量，根据边界修正,并计算缺失量与部分负荷下的热泵功耗、COP
                if Xhp(i,1)<=X(j,4) && Xhp(i,1)>0
                    PL(i,1)=Qhpbui(i,1)/(Qhpsc(i,1)*Xhp02(i,1));                          %情况1，存在部分符合下的热泵
                    COPnom(i,1)=Qhpsc(i,1)/Php(i,1);
                    COPpl(i,1)=COPnom(i,1)*(-0.00006*PL(i,1)^6+0.0017*PL(i,1)^5-0.0181*PL(i,1)^4+0.096*PL(i,1)^3-0.2697*PL(i,1)^2+0.4276*PL(i,1)+0.7626);
                    Phppl(i,1)=Qhpbui(i,1)/COPpl(i,1);
                    Qhpshpl(i,1)=Qhpbui(i,1)+Phppl(i,1);
                elseif Xhp(i,1)>X(j,4)
                    Xhp02(i,1)=fix(X(j,4));
                    Qdeloss(i,1)=Qhpbui(i,1)-Qhpsc(i,1)*Xhp02(i,1);                          %Qdeloss代表缺失的热量
                    Xloss(i,1)=Qdeloss(i,1)/Qhpbui(i,1);
                    Phppl(i,1)=Php(i,1)*Xhp02(i,1);
                    Qhpshpl(i,1)=Xhp02(i,1)*Qhpsh(i,1);
                else
                    Qdeloss(i,1)=Qhpbui(i,1);                          %Qdeloss代表缺失的热量
                    Xloss(i,1)=1;
                    Phppl(i,1)=0;
                    Qhpshpl(i,1)=0;
                end                                                                      %Xhp表示最终的热泵数量，带小数点，情况2 X(j,4)表示数量，即边界值
                
                %计算地埋管数量
                Xge01(i,1)=Qhpshpl(i,1)/(6*0.95);
                Xge(i,1)=fix(Xge(i,1));
                
                if Xge01(i,1)-Xge(i,1) > 0                                               %对地埋管数量取整，向上+1
                    Xge(i,1)=fix(Xge01(i,1))+1;
                else
                    Xge(i,1)=Xge01(i,1);
                end
                
            end
            %%--热泵计算结束--%%
            
            %%--电负荷计算--%%
            Pbss(i,1)=Fuhe(i,1)+Phppl(i,1)-Ppv(i,1)*X(j,1)-Pwt(i,1)*X(j,2);
            
            %蓄电池充放电状态判断
            if Pbss(i,1) < 0                                                            %pbss为负，代表要充电，此时的pbss代表不考虑上下限的功率
                SOC(i+1,1)=(1-0.03)*SOC(i,1)-Pbss(i,1)*0.9/(X(j,3)*Eb);
                More(i,1)=abs(Pbss(i,1));                                               %记录多余的电量
                elecnum(i,1)=0;
            elseif Pbss(i,1) > 0                                                        %pbss为正，代表要放电
                SOC(i+1,1)=(1-0.03)*SOC(i,1)-Pbss(i,1)/(X(j,3)*0.9*Eb);
                More(i,1)=0;                                                            %记录多余的浪费电量
                elecnum(i,1)=abs(Pbss(i,1));
            else
                SOC(i+1,1)=SOC(i,1);
                More(i,1)=0;                                                            %记录多余的浪费电量
                elecnum(i,1)=0;
            end
            
            %多余电量判断
            if SOC(i+1,1)>=SOCmax
                SOC(i+1,1)=SOCmax;
                %先计算即将溢出时的蓄电池真实充电功率
                if SOC(i,1)~=SOCmax                                                         %~=是不等于的意思
                    Pbss(i,1)=(SOC(i+1,1)-(1-0.03)*SOC(i,1))*Eb*X(j,3)/(-0.9);
                else
                    Pbss(i,1)=0;
                end
                %在100%的情况下为弃光与弃风
                if More(i,1)>=abs(Pbss(i,1))
                    More(i,1)=More(i,1)-abs(Pbss(i,1));
                else
                    More(i,1)=0;
                end
                %记录多余电量
            else
                More(i,1)=0;
            end
            
            %缺电判断
            if SOC(i+1,1)<=SOCmin
                SOC(i+1,1)=SOCmin;
                %计算即将见底时的蓄电池真实放电功率
                if SOC(i,1)~=SOCmin
                    Pbss(i,1)=(SOC(i+1,1)-(1-0.03)*SOC(i,1))*Eb*(-0.9)*X(j,3);
                else
                    Pbss(i,1)=0;
                end
                %真实缺电率
                if elecnum(i,1)>=abs(Pbss(i,1))
                    elecnum(i,1)=elecnum(i,1)-abs(Pbss(i,1));
                else
                    elecnum(i,1)=0;
                end
            else
                elecnum(i,1)=0;
            end
            
            %外购电耗费
            Celebuy(i,1)=elecnum(i,1)*DianJia(i,1);
            
            %售电盈利
            Celesell(i,1)=More(i,1)*Cse;
            %Emore(i,1)=More(i,1)/(Ppv(i,1)*X(j,1)+Pwt(i,1)*X(j,2)+0.001);
            %%--电负荷计算结束--%%
            
        end
        %%%---16.00-18.00计算结束---%%%
        
        
        %%%---18.00-23.00，非上班时间，非峰值电价---%%%
        if d(i,1)>0.78 && d(i,1)<0.97
            Ttes(i+1,1)=Ktes*Ates*(Temperature(i,1)-Ttes(i,1))/(mtes*Cw)+Ttes(i,1);
            
            %%%---电负荷计算---%%%
            Pbss(i,1)=Fuhe(i,1)-Ppv(i,1)*X(j,1)-Pwt(i,1)*X(j,2);
            
            %---蓄电池充放电状态判断---%
            if Pbss(i,1) < 0                                                        %pbss为负，代表要充电，此时的pbss代表不考虑上下限的功率
                SOC(i+1,1)=(1-0.03)*SOC(i,1)-Pbss(i,1)*0.9/(X(j,3)*Eb);
                More(i,1)=abs(Pbss(i,1));                                           %记录多余的电量
                elecnum(i,1)=0;
            elseif   Pbss(i,1) > 0                                                  %pbss为正，代表要放电
                SOC(i+1,1)=(1-0.03)*SOC(i,1)-Pbss(i,1)/(X(j,3)*0.9*Eb);
                More(i,1)=0;                                                        %记录多余的浪费电量
                elecnum(i,1)=abs(Pbss(i,1));
            else
                SOC(i+1,1)=SOC(i,1);
                More(i,1)=0;                                                        %记录多余的浪费电量
                elecnum(i,1)=0;
            end
            %---蓄电池充放电状态判断结束---%
            
            %---多余电量判断---%
            if  SOC(i+1,1)>=SOCmax
                SOC(i+1,1)=SOCmax;
                %先计算即将溢出时的蓄电池真实充电功率
                if  SOC(i,1)~=SOCmax                                                %~=是不等于的意思
                    Pbss(i,1)=(SOC(i+1,1)-(1-0.03)*SOC(i,1))*Eb*X(j,3)/(-0.9);
                else
                    Pbss(i,1)=0;
                end
                %记录多余电量
                %在100%的情况下为弃光与弃风
                if More(i,1)>=abs(Pbss(i,1))
                    More(i,1)=More(i,1)-abs(Pbss(i,1));                             %这个设置的很巧妙
                else
                    More(i,1)=0;
                end
                %记录多余电量
            else
                More(i,1)=0;
            end
            %---多余电量判断结束---%
            
            %---缺电判断---%
            if  SOC(i+1,1)<=SOCmin
                SOC(i+1,1)=SOCmin;
                %计算即将见底时的蓄电池真实放电功率
                if  SOC(i,1)~=SOCmin
                    Pbss(i,1)=(SOC(i+1,1)-(1-0.03)*SOC(i,1))*Eb*(-0.9)*X(j,3);
                else
                    Pbss(i,1)=0;
                end
                %在见底时的真实缺电率
                if elecnum(i,1)>=abs(Pbss(i,1))
                    elecnum(i,1)=elecnum(i,1)-abs(Pbss(i,1));                       %这个设置的很巧妙
                else
                    elecnum(i,1)=0;
                end
            else
                elecnum(i,1)=0;
            end
            
            %外购电耗费
            Celebuy(i,1)=0;
            
            %售电盈利
            Celesell(i,1)=More(i,1)*Cse;
            %Emore(i,1)=More(i,1)/(Ppv(i,1)*X(j,1)+Pwt(i,1)*X(j,2)+0.001);
            
        end
        %%%---18.00-23.00计算结束---%%%
        
    end
    %%%---制冷季计算结束---%%%
    %% 过渡季
    for i=5833:7631
        
        %%%---电负荷计算---%%%
        Pbss(i,1)=Fuhe(i,1)-Ppv(i,1)*X(j,1)-Pwt(i,1)*X(j,2);
        
        %---蓄电池充放电状态判断---%
        if Pbss(i,1) < 0                                                        %pbss为负，代表要充电，此时的pbss代表不考虑上下限的功率
            SOC(i+1,1)=(1-0.03)*SOC(i,1)-Pbss(i,1)*0.9/(X(j,3)*Eb);
            More(i,1)=abs(Pbss(i,1));                                           %记录多余的电量
            elecnum(i,1)=0;
        elseif   Pbss(i,1) > 0                                                  %pbss为正，代表要放电
            SOC(i+1,1)=(1-0.03)*SOC(i,1)-Pbss(i,1)/(X(j,3)*0.9*Eb);
            More(i,1)=0;                                                        %记录多余的浪费电量
            elecnum(i,1)=abs(Pbss(i,1));
        else
            SOC(i+1,1)=SOC(i,1);
            More(i,1)=0;                                                        %记录多余的浪费电量
            elecnum(i,1)=0;
        end
        %---蓄电池充放电状态判断结束---%
        
        %---多余电量判断---%
        if  SOC(i+1,1)>=SOCmax
            SOC(i+1,1)=SOCmax;
            %先计算即将溢出时的蓄电池真实充电功率
            if  SOC(i,1)~=SOCmax                                                %~=是不等于的意思
                Pbss(i,1)=(SOC(i+1,1)-(1-0.03)*SOC(i,1))*Eb*X(j,3)/(-0.9);
            else
                Pbss(i,1)=0;
            end
            %记录多余电量
            %在100%的情况下为弃光与弃风
            if More(i,1)>=abs(Pbss(i,1))
                More(i,1)=More(i,1)-abs(Pbss(i,1));                             %这个设置的很巧妙
            else
                More(i,1)=0;
            end
            %记录多余电量
        else
            More(i,1)=0;
        end
        %---多余电量判断结束---%
        
        %---缺电判断---%
        if  SOC(i+1,1)<=SOCmin
            SOC(i+1,1)=SOCmin;
            %计算即将见底时的蓄电池真实放电功率
            if  SOC(i,1)~=SOCmin
                Pbss(i,1)=(SOC(i+1,1)-(1-0.03)*SOC(i,1))*Eb*(-0.9)*X(j,3);
            else
                Pbss(i,1)=0;
            end
            %在见底时的真实缺电率
            if elecnum(i,1)>=abs(Pbss(i,1))
                elecnum(i,1)=elecnum(i,1)-abs(Pbss(i,1));                       %这个设置的很巧妙
            else
                elecnum(i,1)=0;
            end
        else
            elecnum(i,1)=0;
        end
        
        %外购电耗费
        Celebuy(i,1)=elecnum(i,1)*DianJia(i,1);
        
        %售电盈利
        Celesell(i,1)=More(i,1)*Cse;
        %Emore(i,1)=More(i,1)/(Ppv(i,1)*X(j,1)+Pwt(i,1)*X(j,2)+0.001);
        
        %---缺电判断结束---%
        %%%---电模块结束---%%%
    end
    
    %% 供暖季 part02
    for i=7632:8760
        
        a(i,1)=i;
        b(i,1)=a(i,1)/24;
        c(i,1)=fix(b(i,1));
        d(i,1)=b(i,1)-c(i,1);
        
        %%%---23.00-07.00,峰谷时间，热泵向水箱供热---%%%
        if d(i,1)<0.31
            
            if Ttes(i,1)<Ttessdhh                                                                                            %水箱温度控制启停判断，水箱温度<53
                Qhptes(i,1)=((Ttessdh-Ttes(i,1))*Cw*X(j,5)*mtes+X(j,5)*Ktes*Ates*(Ttes(i,1)-Temperature(i,1)))/3600;                %热泵峰谷时期供热热负荷
                
                %%---热泵部分负荷及缺失计算---%%
                
                %第一次计算热泵数量
                Xhp01(i,1)=Qhptes(i,1)/Qhpsh(i,1);
                
                if Xhp01(i,1)<1
                    Xhp02(i,1)=0;
                else
                    Xhp02(i,1)=fix(Xhp01(i,1))+1;
                end
                
                Xhp(i,1)=Xhp02(i,1);
                
                %第二次计算热泵数量，根据边界修正,并计算缺失量与部分负荷下的热泵功耗、COP
                if Xhp(i,1)<=X(j,4) && Xhp(i,1)>0
                    PL(i,1)=Qhptes(i,1)/(Qhpsh(i,1)*Xhp02(i,1));                          %情况1，存在部分符合下的热泵
                    COPnom(i,1)=Qhpsh(i,1)/Php(i,1);
                    COPpl(i,1)=COPnom(i,1)*(-0.00006*PL(i,1)^6+0.0017*PL(i,1)^5-0.0181*PL(i,1)^4+0.096*PL(i,1)^3-0.2697*PL(i,1)^2+0.4276*PL(i,1)+0.7626);
                    Phppl(i,1)=Qhptes(i,1)/COPpl(i,1);                                    %计算部分负荷下该热泵的功耗
                    Qhpscpl(i,1)=Qhptes(i,1)-Phppl(i,1);
                    Ttes(i+1,1)=Ttessdh;
                elseif Xhp(i,1)>X(j,4)
                    Xhp02(i,1)=fix(X(j,4));
                    %Qdeloss(i,1)=Qhptes(i,1)-Qhpsh(i,1)*Xhp02(i,1);                          %Qdeloss代表缺失的热量
                    %Xloss(i,1)=Qdeloss(i,1)/Qhptes(i,1);
                    Phppl(i,1)=Php(i,1)*Xhp02(i,1);
                    Qhpscpl(i,1)=Xhp02(i,1)*Qhpsc(i,1);
                    Ttes(i+1,1)=(Xhp02(i,1)*Qhpsh(i,1)*3600-X(j,5)*Ktes*Ates*(Ttes(i,1)-Temperature(i,1)))/(mtes*Cw*X(j,5))+Ttes(i,1);
                else
                    Phppl(i,1)=0;
                    Qhpscpl(i,1)=0;
                    Ttes(i+1,1)=Ktes*Ates*(Temperature(i,1)-Ttes(i,1))/(mtes*Cw)+Ttes(i,1);
                end                                                                      %Xhp表示最终的热泵数量，带小数点，情况2 X(j,4)表示数量，即边界值
                
                Phptes(i,1)=Phppl(i,1);
                
                %计算地埋管数量
                Xge01(i,1)=Qhpscpl(i,1)/(6*0.95);                                        %地埋管传热效率设置为95%
                Xge(i,1)=fix(Xge01(i,1));
                if Xge01(i,1)-Xge(i,1) > 0                                                 %对地埋管数量取整，向上+1
                    Xge(i,1)=fix(Xge01(i,1))+1;
                else
                    Xge(i,1)=Xge01(i,1);
                end
                
            else
                Ttes(i+1,1)=Ktes*Ates*(Temperature(i,1)-Ttes(i,1))/(mtes*Cw)+Ttes(i,1);
            end
            
            
            %%%---电负荷计算---%%%
            Pbss(i,1)=Fuhe(i,1)-Ppv(i,1)*X(j,1)-Pwt(i,1)*X(j,2);
            
            %---蓄电池充放电状态判断---%
            if Pbss(i,1) < 0
                SOC(i+1,1)=(1-0.03)*SOC(i,1)-Pbss(i,1)*0.9/(X(j,3)*Eb);
                More(i,1)=abs(Pbss(i,1));
                elecnum(i,1)=0;
            elseif   Pbss(i,1) > 0
                SOC(i+1,1)=(1-0.03)*SOC(i,1)-Pbss(i,1)/(X(j,3)*0.9*Eb);
                More(i,1)=0;
                elecnum(i,1)=abs(Pbss(i,1));
            else
                SOC(i+1,1)=SOC(i,1);
                More(i,1)=0;
                elecnum(i,1)=0;
            end
            %---蓄电池充放电状态判断结束---%
            
            %---多余电量判断---%
            if  SOC(i+1,1)>=SOCmax
                SOC(i+1,1)=SOCmax;
                %先计算即将溢出时的蓄电池真实充电功率
                if  SOC(i,1)~=SOCmax
                    Pbss(i,1)=(SOC(i+1,1)-(1-0.03)*SOC(i,1))*Eb*X(j,3)/(-0.9);
                else
                    Pbss(i,1)=0;
                end
                %记录多余电量
                %在100%的情况下为弃光与弃风
                if More(i,1)>=abs(Pbss(i,1))
                    More(i,1)=More(i,1)-abs(Pbss(i,1));
                else
                    More(i,1)=0;
                end
                %记录多余电量
            else
                More(i,1)=0;
            end
            %---多余电量判断结束---%
            
            %---缺电判断---%
            if  SOC(i+1,1)<=SOCmin
                SOC(i+1,1)=SOCmin;
                %计算即将见底时的蓄电池真实放电功率
                if  SOC(i,1)~=SOCmin
                    Pbss(i,1)=(SOC(i+1,1)-(1-0.03)*SOC(i,1))*Eb*(-0.9)*X(j,3);
                else
                    Pbss(i,1)=0;
                end
                %在见底时的真实缺电率
                if elecnum(i,1)>=abs(Pbss(i,1))
                    elecnum(i,1)=elecnum(i,1)-abs(Pbss(i,1));
                else
                    elecnum(i,1)=0;
                end
            else
                elecnum(i,1)=0;
            end
            
            %外购电耗费
            Celebuy(i,1)=Phptes(i,1)*DianJia(i,1);
            
            %售电盈利
            Celesell(i,1)=More(i,1)*Cse;
            %Emore(i,1)=More(i,1)/(Ppv(i,1)*X(j,1)+Pwt(i,1)*X(j,2)+0.001);
            
            %---缺电判断结束---%
            %%%---电模块结束---%%%
        end
        
        %%%---23.00-07.00计算结束---%%%
        
        %%%---07.00-08.00计算，非上班时间，非峰谷时间---%%%
        if d(i,1)>0.31 && d(i,1)<0.36
            
            Ttes(i+1,1)=Ktes*Ates*(Temperature(i,1)-Ttes(i,1))/(mtes*Cw)+Ttes(i,1);
            
            %%%---电负荷计算---%%%
            Pbss(i,1)=Fuhe(i,1)-Ppv(i,1)*X(j,1)-Pwt(i,1)*X(j,2);
            
            %---蓄电池充放电状态判断---%
            if Pbss(i,1) < 0
                SOC(i+1,1)=(1-0.03)*SOC(i,1)-Pbss(i,1)*0.9/(X(j,3)*Eb);
                More(i,1)=abs(Pbss(i,1));
                elecnum(i,1)=0;
            elseif   Pbss(i,1) > 0
                SOC(i+1,1)=(1-0.03)*SOC(i,1)-Pbss(i,1)/(X(j,3)*0.9*Eb);
                More(i,1)=0;
                elecnum(i,1)=abs(Pbss(i,1));
            else
                SOC(i+1,1)=SOC(i,1);
                More(i,1)=0;
                elecnum(i,1)=0;
            end
            %---蓄电池充放电状态判断结束---%
            
            %---多余电量判断---%
            if  SOC(i+1,1)>=SOCmax
                SOC(i+1,1)=SOCmax;
                %先计算即将溢出时的蓄电池真实充电功率
                if  SOC(i,1)~=SOCmax
                    Pbss(i,1)=(SOC(i+1,1)-(1-0.03)*SOC(i,1))*Eb*X(j,3)/(-0.9);
                else
                    Pbss(i,1)=0;
                end
                %记录多余电量
                %在100%的情况下为弃光与弃风
                if More(i,1)>=abs(Pbss(i,1))
                    More(i,1)=More(i,1)-abs(Pbss(i,1));
                else
                    More(i,1)=0;
                end
                %记录多余电量
            else
                More(i,1)=0;
            end
            %---多余电量判断结束---%
            
            %---缺电判断---%
            if  SOC(i+1,1)<=SOCmin
                SOC(i+1,1)=SOCmin;
                %计算即将见底时的蓄电池真实放电功率
                if  SOC(i,1)~=SOCmin
                    Pbss(i,1)=(SOC(i+1,1)-(1-0.03)*SOC(i,1))*Eb*(-0.9)*X(j,3);
                else
                    Pbss(i,1)=0;
                end
                %在见底时的真实缺电率
                if elecnum(i,1)>=abs(Pbss(i,1))
                    elecnum(i,1)=elecnum(i,1)-abs(Pbss(i,1));
                else
                    elecnum(i,1)=0;
                end
            else
                elecnum(i,1)=0;
            end
            
            %外购电耗费
            Celebuy(i,1)=0;
            
            %售电盈利
            Celesell(i,1)=More(i,1)*Cse;
            %Emore(i,1)=More(i,1)/(Ppv(i,1)*X(j,1)+Pwt(i,1)*X(j,2)+0.001);
            
        end
        %%%---07.00-08.00计算结束---%%%
        
        %%%---8.00-9.00,工作时间，平常电价---%%%
        if d(i,1)>0.36 && d(i,1)<0.39
            
            %%--水箱温度计算--%%
            Ttes(i+1,1)=Ktes*Ates*(Temperature(i,1)-Ttes(i,1))/(mtes*Cw)+Ttes(i,1);
            %%--水箱温度计算结束--%%
            
            %%--热泵计算--%%
            %第一次计算热泵数量
            Xhp01(i,1)=HeatingFuhe(i,1)/Qhpsh(i,1);
            
            if Xhp01(i,1)<1
                Xhp02(i,1)=0;
            else
                Xhp02(i,1)=fix(Xhp01(i,1))+1;
            end
            
            Xhp(i,1)=Xhp02(i,1);
                
            %第二次计算热泵数量，根据边界修正,并计算缺失量与部分负荷下的热泵功耗、COP
            if Xhp(i,1)<=X(j,4) && Xhp(i,1)>0
                PL(i,1)=HeatingFuhe(i,1)/(Qhpsh(i,1)*Xhp02(i,1));                          %情况1，存在部分符合下的热泵
                COPnom(i,1)=Qhpsh(i,1)/Php(i,1);
                COPpl(i,1)=COPnom(i,1)*(-0.00006*PL(i,1)^6+0.0017*PL(i,1)^5-0.0181*PL(i,1)^4+0.096*PL(i,1)^3-0.2697*PL(i,1)^2+0.4276*PL(i,1)+0.7626);
                Phppl(i,1)=HeatingFuhe(i,1)/COPpl(i,1);
                Qhpscpl(i,1)=HeatingFuhe(i,1)-Phppl(i,1);
            elseif Xhp(i,1)>X(j,4)
                Xhp02(i,1)=fix(X(j,4));
                Qdeloss(i,1)=HeatingFuhe(i,1)-Qhpsh(i,1)*Xhp02(i,1);                          %Qdeloss代表缺失的热量
                Xloss(i,1)=Qdeloss(i,1)/HeatingFuhe(i,1);
                Phppl(i,1)=Php(i,1)*Xhp02(i,1);
                Qhpscpl(i,1)=Xhp02(i,1)*Qhpsc(i,1);
            else
                Qdeloss(i,1)=HeatingFuhe(i,1);                          %Qdeloss代表缺失的热量
                Xloss(i,1)=1;
                Phppl(i,1)=0;
                Qhpscpl(i,1)=0;
            end                                                                      %Xhp表示最终的热泵数量，带小数点，情况2 X(j,4)表示数量，即边界值
            
            %计算地埋管的数量
            Xge01(i,1)=Qhpscpl(i,1)/(6*0.95);
            Xge(i,1)=fix(Xge(i,1));
            
            if Xge01(i,1)-Xge(i,1) > 0                                                %对地埋管数量取整，向上+1
                Xge(i,1)=fix(Xge01(i,1))+1;
            else
                Xge(i,1)=Xge01(i,1);
            end
            %%--热泵计算结束--%%
            
            %%--电负荷计算--%%
            Pbss(i,1)=Fuhe(i,1)+Phppl(i,1)-Ppv(i,1)*X(j,1)-Pwt(i,1)*X(j,2);
            
            %蓄电池充放电状态判断
            if Pbss(i,1) < 0                                                         %pbss为负，代表要充电，此时的pbss代表不考虑上下限的功率
                SOC(i+1,1)=(1-0.03)*SOC(i,1)-Pbss(i,1)*0.9/(X(j,3)*Eb);
                More(i,1)=abs(Pbss(i,1));                                            %记录多余的电量
                elecnum(i,1)=0;
            elseif Pbss(i,1) > 0                                                     %pbss为正，代表要放电
                SOC(i+1,1)=(1-0.03)*SOC(i,1)-Pbss(i,1)/(X(j,3)*0.9*Eb);
                More(i,1)=0;                                                         %记录多余的浪费电量
                elecnum(i,1)=abs(Pbss(i,1));
            else
                SOC(i+1,1)=SOC(i,1);
                More(i,1)=0;                                                         %记录多余的浪费电量
                elecnum(i,1)=0;
            end
            
            %多余电量判断
            if SOC(i+1,1)>=SOCmax
                SOC(i+1,1)=SOCmax;
                %先计算即将溢出时的蓄电池真实充电功率
                if SOC(i,1)~=SOCmax                                                      %~=是不等于的意思
                    Pbss(i,1)=(SOC(i+1,1)-(1-0.03)*SOC(i,1))*Eb*X(j,3)/(-0.9);
                else
                    Pbss(i,1)=0;
                end
                %在100%的情况下为弃光与弃风
                if More(i,1)>=abs(Pbss(i,1))
                    More(i,1)=More(i,1)-abs(Pbss(i,1));
                else
                    More(i,1)=0;
                end
            else
                More(i,1)=0;
            end
            
            %缺电判断
            if SOC(i+1,1)<=SOCmin
                SOC(i+1,1)=SOCmin;
                %计算即将见底时的蓄电池真实放电功率
                if SOC(i,1)~=SOCmin
                    Pbss(i,1)=(SOC(i+1,1)-(1-0.03)*SOC(i,1))*Eb*(-0.9)*X(j,3);
                else
                    Pbss(i,1)=0;
                end
                %真实缺电率
                if elecnum(i,1)>=abs(Pbss(i,1))
                    elecnum(i,1)=elecnum(i,1)-abs(Pbss(i,1));
                else
                    elecnum(i,1)=0;
                end
            else
                elecnum(i,1)=0;
            end
            
            %外购电耗费
            Celebuy(i,1)=elecnum(i,1)*DianJia(i,1);
            
            %售电盈利
            Celesell(i,1)=More(i,1)*Cse;
            %Emore(i,1)=More(i,1)/(Ppv(i,1)*X(j,1)+Pwt(i,1)*X(j,2)+0.001);
            %%--电负荷计算结束--%%
            
        end
        %%%---8.00-9.00计算结束---%%%
        
        %%%---9.00-12.00,工作时间，峰值电价，开启水箱---%%%
        if d(i,1)>0.39 && d(i,1)<0.51
            
            %%--热泵计算--%%
            Qtes(i,1)=((Ttes(i,1)-Tbui)*Cw*X(j,5)*mtes-X(j,5)*Ktes*Ates*(Ttes(i,1)-Temperature(i,1)))/3600;
            
            %防止后续增加热泵负荷
            if Qtes(i,1)<= 0
                Qtes(i,1)=0;
                
                %%--水箱温度计算--%%
                Ttes(i+1,1)=Ktes*Ates*(Temperature(i,1)-Ttes(i,1))/(mtes*Cw)+Ttes(i,1);
                %%--水箱温度计算结束--%%
            
                %%--热泵计算--%%
                %第一次计算热泵数量
                Xhp01(i,1)=HeatingFuhe(i,1)/Qhpsh(i,1);
                
                if Xhp01(i,1)<1
                    Xhp02(i,1)=0;
                else
                    Xhp02(i,1)=fix(Xhp01(i,1))+1;
                end
                
                Xhp(i,1)=Xhp02(i,1);
                
                %第二次计算热泵数量，根据边界修正,并计算缺失量与部分负荷下的热泵功耗、COP
                if Xhp(i,1)<=X(j,4) && Xhp(i,1)>0
                    PL(i,1)=HeatingFuhe(i,1)/(Qhpsh(i,1)*Xhp02(i,1));                          %情况1，存在部分符合下的热泵
                    COPnom(i,1)=Qhpsh(i,1)/Php(i,1);
                    COPpl(i,1)=COPnom(i,1)*(-0.00006*PL(i,1)^6+0.0017*PL(i,1)^5-0.0181*PL(i,1)^4+0.096*PL(i,1)^3-0.2697*PL(i,1)^2+0.4276*PL(i,1)+0.7626);
                    Phppl(i,1)=HeatingFuhe(i,1)/COPpl(i,1);
                    Qhpscpl(i,1)=HeatingFuhe(i,1)-Phppl(i,1);
                elseif Xhp(i,1)>X(j,4)
                    Xhp02(i,1)=fix(X(j,4));
                    Qdeloss(i,1)=HeatingFuhe(i,1)-Qhpsh(i,1)*Xhp02(i,1);                          %Qdeloss代表缺失的热量
                    Xloss(i,1)=Qdeloss(i,1)/HeatingFuhe(i,1);
                    Phppl(i,1)=Php(i,1)*Xhp02(i,1);
                    Qhpscpl(i,1)=Xhp02(i,1)*Qhpsc(i,1);
                else
                    Qdeloss(i,1)=HeatingFuhe(i,1);                          %Qdeloss代表缺失的热量
                    Xloss(i,1)=1;
                    Phppl(i,1)=0;
                    Qhpscpl(i,1)=0;
                end                                                                      %Xhp表示最终的热泵数量，带小数点，情况2 X(j,4)表示数量，即边界值
            
                %计算地埋管的数量
                Xge01(i,1)=Qhpscpl(i,1)/(6*0.95);
                Xge(i,1)=fix(Xge(i,1));
            
                if Xge01(i,1)-Xge(i,1) > 0                                                %对地埋管数量取整，向上+1
                    Xge(i,1)=fix(Xge01(i,1))+1;
                else
                    Xge(i,1)=Xge01(i,1);
                end
                %%--热泵计算结束--%%
            end
            
            %情况1/2/3
            if Qtes(i,1) >= HeatingFuhe(i,1)
                Ttes(i+1,1)=(X(j,5)*Ktes*Ates*(Temperature(i,1)-Ttes(i,1))-HeatingFuhe(i,1)*3600)/(X(j,5)*mtes*Cw)+Ttes(i,1);
                Xhp02(i,1)=0;
            elseif Qtes(i,1)<HeatingFuhe(i,1) && Qtes(i,1)>0
                Qhpbui(i,1)=HeatingFuhe(i,1)-Qtes(i,1);
                Ttes(i+1,1)=Tbui;
                
                %情况2/3
                %第一次计算热泵数量
                Xhp01(i,1)=Qhpbui(i,1)/Qhpsh(i,1);
                
                if Xhp01(i,1)<1
                    Xhp02(i,1)=0;
                else
                    Xhp02(i,1)=fix(Xhp01(i,1))+1;
                end
                
                Xhp(i,1)=Xhp02(i,1);
                
                %第二次计算热泵数量，根据边界修正,并计算缺失量与部分负荷下的热泵功耗、COP
                if Xhp(i,1) <= X(j,4) && Xhp(i,1)>0
                    PL(i,1)=Qhpbui(i,1)/(Qhpsh(i,1)*Xhp02(i,1));                          %情况1，存在部分符合下的热泵
                    COPnom(i,1)=Qhpsh(i,1)/Php(i,1);
                    COPpl(i,1)=COPnom(i,1)*(-0.00006*PL(i,1)^6+0.0017*PL(i,1)^5-0.0181*PL(i,1)^4+0.096*PL(i,1)^3-0.2697*PL(i,1)^2+0.4276*PL(i,1)+0.7626);
                    Phppl(i,1)=Qhpbui(i,1)/COPpl(i,1);
                    Qhpscpl(i,1)=Qhpbui(i,1)-Phppl(i,1);
                elseif Xhp(i,1)> X(j,4)
                    Xhp02(i,1)=fix(X(j,4));
                    Qdeloss(i,1)=Qhpbui(i,1)-Qhpsh(i,1)*Xhp02(i,1);                          %Qdeloss代表缺失的热量
                    Xloss(i,1)=Qdeloss(i,1)/Qhpbui(i,1);
                    Phppl(i,1)=Php(i,1)*Xhp02(i,1);
                    Qhpscpl(i,1)=Xhp02(i,1)*Qhpsc(i,1);
                else
                    Qdeloss(i,1)=Qhpbui(i,1);                          %Qdeloss代表缺失的热量
                    Xloss(i,1)=1;
                    Phppl(i,1)=0;
                    Qhpscpl(i,1)=0;
                end                                                                      %Xhp表示最终的热泵数量，带小数点，情况2 X(j,4)表示数量，即边界值
                
                %计算地埋管数量
                Xge01(i,1)=Qhpscpl(i,1)/(6*0.95);
                Xge(i,1)=fix(Xge(i,1));
                
                if Xge01(i,1)-Xge(i,1) > 0                                               %对地埋管数量取整，向上+1
                    Xge(i,1)=fix(Xge01(i,1))+1;
                else
                    Xge(i,1)=Xge01(i,1);
                end
                
            end
            %%--热泵计算结束--%%
            
            %%--电负荷计算--%%
            Pbss(i,1)=Fuhe(i,1)+Phppl(i,1)-Ppv(i,1)*X(j,1)-Pwt(i,1)*X(j,2);
            
            %蓄电池充放电状态判断
            if Pbss(i,1) < 0                                                            %pbss为负，代表要充电，此时的pbss代表不考虑上下限的功率
                SOC(i+1,1)=(1-0.03)*SOC(i,1)-Pbss(i,1)*0.9/(X(j,3)*Eb);
                More(i,1)=abs(Pbss(i,1));                                               %记录多余的电量
                elecnum(i,1)=0;
            elseif Pbss(i,1) > 0                                                        %pbss为正，代表要放电
                SOC(i+1,1)=(1-0.03)*SOC(i,1)-Pbss(i,1)/(X(j,3)*0.9*Eb);
                More(i,1)=0;                                                            %记录多余的浪费电量
                elecnum(i,1)=abs(Pbss(i,1));
            else
                SOC(i+1,1)=SOC(i,1);
                More(i,1)=0;                                                            %记录多余的浪费电量
                elecnum(i,1)=0;
            end
            
            %多余电量判断
            if SOC(i+1,1)>=SOCmax
                SOC(i+1,1)=SOCmax;
                %先计算即将溢出时的蓄电池真实充电功率
                if SOC(i,1)~=SOCmax                                                         %~=是不等于的意思
                    Pbss(i,1)=(SOC(i+1,1)-(1-0.03)*SOC(i,1))*Eb*X(j,3)/(-0.9);
                else
                    Pbss(i,1)=0;
                end
                %在100%的情况下为弃光与弃风
                if More(i,1)>=abs(Pbss(i,1))
                    More(i,1)=More(i,1)-abs(Pbss(i,1));
                else
                    More(i,1)=0;
                end
                %记录多余电量
            else
                More(i,1)=0;
            end
            
            %缺电判断
            if SOC(i+1,1)<=SOCmin
                SOC(i+1,1)=SOCmin;
                %计算即将见底时的蓄电池真实放电功率
                if SOC(i,1)~=SOCmin
                    Pbss(i,1)=(SOC(i+1,1)-(1-0.03)*SOC(i,1))*Eb*(-0.9)*X(j,3);
                else
                    Pbss(i,1)=0;
                end
                %真实缺电率
                if elecnum(i,1)>=abs(Pbss(i,1))
                    elecnum(i,1)=elecnum(i,1)-abs(Pbss(i,1));
                else
                    elecnum(i,1)=0;
                end
            else
                elecnum(i,1)=0;
            end
            
            %外购电耗费
            Celebuy(i,1)=elecnum(i,1)*DianJia(i,1);
            
            %售电盈利
            Celesell(i,1)=More(i,1)*Cse;
            %Emore(i,1)=More(i,1)/(Ppv(i,1)*X(j,1)+Pwt(i,1)*X(j,2)+0.001);
            %%--电负荷计算结束--%%
            
        end
        %%%---9.00-12.00计算结束---%%%
        
        %%%---12.00-16.00，工作时间，平常电价---%%%
        if d(i,1)>0.51 && d(i,1)<0.69
            
            %%--水箱温度计算--%%
            Ttes(i+1,1)=Ktes*Ates*(Temperature(i,1)-Ttes(i,1))/(mtes*Cw)+Ttes(i,1);
            %%--水箱温度计算结束--%%
            
            %%--热泵计算--%%
            %第一次计算热泵数量
            Xhp01(i,1)=HeatingFuhe(i,1)/Qhpsh(i,1);
            
            if Xhp01(i,1)<1
                Xhp02(i,1)=0;
            else
                Xhp02(i,1)=fix(Xhp01(i,1))+1;
            end
            
            Xhp(i,1)=Xhp02(i,1);
                
            %第二次计算热泵数量，根据边界修正,并计算缺失量与部分负荷下的热泵功耗、COP
            if Xhp(i,1)<=X(j,4) && Xhp(i,1)>0
                PL(i,1)=HeatingFuhe(i,1)/(Qhpsh(i,1)*Xhp02(i,1));                          %情况1，存在部分符合下的热泵
                COPnom(i,1)=Qhpsh(i,1)/Php(i,1);
                COPpl(i,1)=COPnom(i,1)*(-0.00006*PL(i,1)^6+0.0017*PL(i,1)^5-0.0181*PL(i,1)^4+0.096*PL(i,1)^3-0.2697*PL(i,1)^2+0.4276*PL(i,1)+0.7626);
                Phppl(i,1)=HeatingFuhe(i,1)/COPpl(i,1);
                Qhpscpl(i,1)=HeatingFuhe(i,1)-Phppl(i,1);
            elseif Xhp(i,1)>X(j,4)
                Xhp02(i,1)=fix(X(j,4));
                Qdeloss(i,1)=HeatingFuhe(i,1)-Qhpsh(i,1)*Xhp02(i,1);                          %Qdeloss代表缺失的热量
                Xloss(i,1)=Qdeloss(i,1)/HeatingFuhe(i,1);
                Phppl(i,1)=Php(i,1)*Xhp02(i,1);
                Qhpscpl(i,1)=Xhp02(i,1)*Qhpsc(i,1);
            else
                Qdeloss(i,1)=HeatingFuhe(i,1);                          %Qdeloss代表缺失的热量
                Xloss(i,1)=1;
                Phppl(i,1)=0;
                Qhpscpl(i,1)=0;
            end                                                                      %Xhp表示最终的热泵数量，带小数点，情况2 X(j,4)表示数量，即边界值
            
            %计算地埋管的数量
            Xge01(i,1)=Qhpscpl(i,1)/(6*0.95);
            Xge(i,1)=fix(Xge(i,1));
            
            if Xge01(i,1)-Xge(i,1) > 0                                                %对地埋管数量取整，向上+1
                Xge(i,1)=fix(Xge01(i,1))+1;
            else
                Xge(i,1)=Xge01(i,1);
            end
            %%--热泵计算结束--%%
            
            %%--电负荷计算--%%
            Pbss(i,1)=Fuhe(i,1)+Phppl(i,1)-Ppv(i,1)*X(j,1)-Pwt(i,1)*X(j,2);
            
            %蓄电池充放电状态判断
            if Pbss(i,1) < 0                                                         %pbss为负，代表要充电，此时的pbss代表不考虑上下限的功率
                SOC(i+1,1)=(1-0.03)*SOC(i,1)-Pbss(i,1)*0.9/(X(j,3)*Eb);
                More(i,1)=abs(Pbss(i,1));                                            %记录多余的电量
                elecnum(i,1)=0;
            elseif Pbss(i,1) > 0                                                     %pbss为正，代表要放电
                SOC(i+1,1)=(1-0.03)*SOC(i,1)-Pbss(i,1)/(X(j,3)*0.9*Eb);
                More(i,1)=0;                                                         %记录多余的浪费电量
                elecnum(i,1)=abs(Pbss(i,1));
            else
                SOC(i+1,1)=SOC(i,1);
                More(i,1)=0;                                                         %记录多余的浪费电量
                elecnum(i,1)=0;
            end
            
            %多余电量判断
            if SOC(i+1,1)>=SOCmax
                SOC(i+1,1)=SOCmax;
                %先计算即将溢出时的蓄电池真实充电功率
                if SOC(i,1)~=SOCmax                                                      %~=是不等于的意思
                    Pbss(i,1)=(SOC(i+1,1)-(1-0.03)*SOC(i,1))*Eb*X(j,3)/(-0.9);
                else
                    Pbss(i,1)=0;
                end
                %在100%的情况下为弃光与弃风
                if More(i,1)>=abs(Pbss(i,1))
                    More(i,1)=More(i,1)-abs(Pbss(i,1));
                else
                    More(i,1)=0;
                end
                %记录多余电量
            else
                More(i,1)=0;
            end
            
            %缺电判断
            if SOC(i+1,1)<=SOCmin
                SOC(i+1,1)=SOCmin;
                %计算即将见底时的蓄电池真实放电功率
                if SOC(i,1)~=SOCmin
                    Pbss(i,1)=(SOC(i+1,1)-(1-0.03)*SOC(i,1))*Eb*(-0.9)*X(j,3);
                else
                    Pbss(i,1)=0;
                end
                %真实缺电率
                if elecnum(i,1)>=abs(Pbss(i,1))
                    elecnum(i,1)=elecnum(i,1)-abs(Pbss(i,1));
                else
                    elecnum(i,1)=0;
                end
            else
                elecnum(i,1)=0;
            end
            
            %外购电耗费
            Celebuy(i,1)=elecnum(i,1)*DianJia(i,1);
            
            %售电盈利
            Celesell(i,1)=More(i,1)*Cse;
            %Emore(i,1)=More(i,1)/(Ppv(i,1)*X(j,1)+Pwt(i,1)*X(j,2)+0.001);
            %%--电负荷计算结束--%%
            
        end
        %%%---12.00-16.00计算结束---%%%
        
        %%%---16.00-18.00,工作时间，峰值电价，开启水箱---%%%
        if d(i,1)>0.69 && d(i,1)<0.78
            
            %%--热泵计算--%%
            Qtes(i,1)=((Ttes(i,1)-Tbui)*Cw*X(j,5)*mtes-X(j,5)*Ktes*Ates*(Ttes(i,1)-Temperature(i,1)))/3600;
            
            %防止后续增加热泵负荷
            if Qtes(i,1)<0
                Qtes(i,1)=0;
                
                %%--水箱温度计算--%%
                Ttes(i+1,1)=Ktes*Ates*(Temperature(i,1)-Ttes(i,1))/(mtes*Cw)+Ttes(i,1);
                %%--水箱温度计算结束--%%
            
                %%--热泵计算--%%
                %第一次计算热泵数量
                Xhp01(i,1)=HeatingFuhe(i,1)/Qhpsh(i,1);
                
                if Xhp01(i,1)<1
                    Xhp02(i,1)=0;
                else
                    Xhp02(i,1)=fix(Xhp01(i,1))+1;
                end
                
                Xhp(i,1)=Xhp02(i,1);
                
                %第二次计算热泵数量，根据边界修正,并计算缺失量与部分负荷下的热泵功耗、COP
                if Xhp(i,1)<=X(j,4) && Xhp(i,1)>0
                    PL(i,1)=HeatingFuhe(i,1)/(Qhpsh(i,1)*Xhp02(i,1));                          %情况1，存在部分符合下的热泵
                    COPnom(i,1)=Qhpsh(i,1)/Php(i,1);
                    COPpl(i,1)=COPnom(i,1)*(-0.00006*PL(i,1)^6+0.0017*PL(i,1)^5-0.0181*PL(i,1)^4+0.096*PL(i,1)^3-0.2697*PL(i,1)^2+0.4276*PL(i,1)+0.7626);
                    Phppl(i,1)=HeatingFuhe(i,1)/COPpl(i,1);
                    Qhpscpl(i,1)=HeatingFuhe(i,1)-Phppl(i,1);
                elseif Xhp(i,1)>X(j,4)
                    Xhp02(i,1)=fix(X(j,4));
                    Qdeloss(i,1)=HeatingFuhe(i,1)-Qhpsh(i,1)*Xhp02(i,1);                          %Qdeloss代表缺失的热量
                    Xloss(i,1)=Qdeloss(i,1)/HeatingFuhe(i,1);
                    Phppl(i,1)=Php(i,1)*Xhp02(i,1);
                    Qhpscpl(i,1)=Xhp02(i,1)*Qhpsc(i,1);
                else
                    Qdeloss(i,1)=HeatingFuhe(i,1);                          %Qdeloss代表缺失的热量
                    Xloss(i,1)=1;
                    Phppl(i,1)=0;
                    Qhpscpl(i,1)=0;
                end                                                                      %Xhp表示最终的热泵数量，带小数点，情况2 X(j,4)表示数量，即边界值
            
                %计算地埋管的数量
                Xge01(i,1)=Qhpscpl(i,1)/(6*0.95);
                Xge(i,1)=fix(Xge(i,1));
            
                if Xge01(i,1)-Xge(i,1) > 0                                                %对地埋管数量取整，向上+1
                    Xge(i,1)=fix(Xge01(i,1))+1;
                else
                    Xge(i,1)=Xge01(i,1);
                end
                %%--热泵计算结束--%%
            end
            
            %情况1/2/3
            if Qtes(i,1)>=HeatingFuhe(i,1)
                Ttes(i+1,1)=(X(j,5)*Ktes*Ates*(Temperature(i,1)-Ttes(i,1))-HeatingFuhe(i,1)*3600)/(X(j,5)*mtes*Cw)+Ttes(i,1);
                Xhp02(i,1)=0;
            elseif Qtes(i,1)<HeatingFuhe(i,1) && Qtes(i,1)>0
                Qhpbui(i,1)=HeatingFuhe(i,1)-Qtes(i,1);
                Ttes(i+1,1)=Tbui;
                
                %情况2/3
                %第一次计算热泵数量
                Xhp01(i,1)=Qhpbui(i,1)/Qhpsh(i,1);
                
                if Xhp01(i,1)<1
                    Xhp02(i,1)=0;
                else
                    Xhp02(i,1)=fix(Xhp01(i,1))+1;
                end
                
                Xhp(i,1)=Xhp02(i,1);
                
                %第二次计算热泵数量，根据边界修正,并计算缺失量与部分负荷下的热泵功耗、COP
                if Xhp(i,1) <= X(j,4) && Xhp(i,1)>0
                    PL(i,1)=Qhpbui(i,1)/(Qhpsh(i,1)*Xhp02(i,1));                          %情况1，存在部分符合下的热泵
                    COPnom(i,1)=Qhpsh(i,1)/Php(i,1);
                    COPpl(i,1)=COPnom(i,1)*(-0.00006*PL(i,1)^6+0.0017*PL(i,1)^5-0.0181*PL(i,1)^4+0.096*PL(i,1)^3-0.2697*PL(i,1)^2+0.4276*PL(i,1)+0.7626);
                    Phppl(i,1)=Qhpbui(i,1)/COPpl(i,1);
                    Qhpscpl(i,1)=Qhpbui(i,1)-Phppl(i,1);
                elseif Xhp(i,1)> X(j,4)
                    Xhp02(i,1)=fix(X(j,4));
                    Qdeloss(i,1)=Qhpbui(i,1)-Qhpsh(i,1)*Xhp02(i,1);                          %Qdeloss代表缺失的热量
                    Xloss(i,1)=Qdeloss(i,1)/Qhpbui(i,1);
                    Phppl(i,1)=Php(i,1)*Xhp02(i,1);
                    Qhpscpl(i,1)=Xhp02(i,1)*Qhpsc(i,1);
                else
                    Qdeloss(i,1)=Qhpbui(i,1);                          %Qdeloss代表缺失的热量
                    Xloss(i,1)=1;
                    Phppl(i,1)=0;
                    Qhpscpl(i,1)=0;
                end                                                                      %Xhp表示最终的热泵数量，带小数点，情况2 X(j,4)表示数量，即边界值
                
                %计算地埋管数量
                Xge01(i,1)=Qhpscpl(i,1)/(6*0.95);
                Xge(i,1)=fix(Xge(i,1));
                
                if Xge01(i,1)-Xge(i,1) > 0                                               %对地埋管数量取整，向上+1
                    Xge(i,1)=fix(Xge01(i,1))+1;
                else
                    Xge(i,1)=Xge01(i,1);
                end
                
            end
            %%--热泵计算结束--%%
            
            %%--电负荷计算--%%
            Pbss(i,1)=Fuhe(i,1)+Phppl(i,1)-Ppv(i,1)*X(j,1)-Pwt(i,1)*X(j,2);
            
            %蓄电池充放电状态判断
            if Pbss(i,1) < 0                                                            %pbss为负，代表要充电，此时的pbss代表不考虑上下限的功率
                SOC(i+1,1)=(1-0.03)*SOC(i,1)-Pbss(i,1)*0.9/(X(j,3)*Eb);
                More(i,1)=abs(Pbss(i,1));                                               %记录多余的电量
                elecnum(i,1)=0;
            elseif Pbss(i,1) > 0                                                        %pbss为正，代表要放电
                SOC(i+1,1)=(1-0.03)*SOC(i,1)-Pbss(i,1)/(X(j,3)*0.9*Eb);
                More(i,1)=0;                                                            %记录多余的浪费电量
                elecnum(i,1)=abs(Pbss(i,1));
            else
                SOC(i+1,1)=SOC(i,1);
                More(i,1)=0;                                                            %记录多余的浪费电量
                elecnum(i,1)=0;
            end
            
            %多余电量判断
            if SOC(i+1,1)>=SOCmax
                SOC(i+1,1)=SOCmax;
                %先计算即将溢出时的蓄电池真实充电功率
                if SOC(i,1)~=SOCmax                                                         %~=是不等于的意思
                    Pbss(i,1)=(SOC(i+1,1)-(1-0.03)*SOC(i,1))*Eb*X(j,3)/(-0.9);
                else
                    Pbss(i,1)=0;
                end
                %在100%的情况下为弃光与弃风
                if More(i,1)>=abs(Pbss(i,1))
                    More(i,1)=More(i,1)-abs(Pbss(i,1));
                else
                    More(i,1)=0;
                end
                %记录多余电量
            else
                More(i,1)=0;
            end
            
            %缺电判断
            if SOC(i+1,1)<=SOCmin
                SOC(i+1,1)=SOCmin;
                %计算即将见底时的蓄电池真实放电功率
                if SOC(i,1)~=SOCmin
                    Pbss(i,1)=(SOC(i+1,1)-(1-0.03)*SOC(i,1))*Eb*(-0.9)*X(j,3);
                else
                    Pbss(i,1)=0;
                end
                %真实缺电率
                if elecnum(i,1)>=abs(Pbss(i,1))
                    elecnum(i,1)=elecnum(i,1)-abs(Pbss(i,1));
                else
                    elecnum(i,1)=0;
                end
            else
                elecnum(i,1)=0;
            end
            
            %外购电耗费
            Celebuy(i,1)=elecnum(i,1)*DianJia(i,1);
            
            %售电盈利
            Celesell(i,1)=More(i,1)*Cse;
            %Emore(i,1)=More(i,1)/(Ppv(i,1)*X(j,1)+Pwt(i,1)*X(j,2)+0.001);
            %%--电负荷计算结束--%%
            
        end
        %%%---16.00-18.00计算结束---%%%
        
        
        %%%---18.00-23.00，非上班时间，非峰谷电价---%%%
        if d(i,1)>0.78 && d(i,1)<0.97
            Ttes(i+1,1)=Ktes*Ates*(Temperature(i,1)-Ttes(i,1))/(mtes*Cw)+Ttes(i,1);
            
            %%%---电负荷计算---%%%
            Pbss(i,1)=Fuhe(i,1)-Ppv(i,1)*X(j,1)-Pwt(i,1)*X(j,2);
            
            %---蓄电池充放电状态判断---%
            if Pbss(i,1) < 0                                                        %pbss为负，代表要充电，此时的pbss代表不考虑上下限的功率
                SOC(i+1,1)=(1-0.03)*SOC(i,1)-Pbss(i,1)*0.9/(X(j,3)*Eb);
                More(i,1)=abs(Pbss(i,1));                                           %记录多余的电量
                elecnum(i,1)=0;
            elseif   Pbss(i,1) > 0                                                  %pbss为正，代表要放电
                SOC(i+1,1)=(1-0.03)*SOC(i,1)-Pbss(i,1)/(X(j,3)*0.9*Eb);
                More(i,1)=0;                                                        %记录多余的浪费电量
                elecnum(i,1)=abs(Pbss(i,1));
            else
                SOC(i+1,1)=SOC(i,1);
                More(i,1)=0;                                                        %记录多余的浪费电量
                elecnum(i,1)=0;
            end
            %---蓄电池充放电状态判断结束---%
            
            %---多余电量判断---%
            if  SOC(i+1,1)>=SOCmax
                SOC(i+1,1)=SOCmax;
                %先计算即将溢出时的蓄电池真实充电功率
                if  SOC(i,1)~=SOCmax                                                %~=是不等于的意思
                    Pbss(i,1)=(SOC(i+1,1)-(1-0.03)*SOC(i,1))*Eb*X(j,3)/(-0.9);
                else
                    Pbss(i,1)=0;
                end
                %记录多余电量
                %在100%的情况下为弃光与弃风
                if More(i,1)>=abs(Pbss(i,1))
                    More(i,1)=More(i,1)-abs(Pbss(i,1));                             %这个设置的很巧妙
                else
                    More(i,1)=0;
                end
                %记录多余电量
            else
                More(i,1)=0;
            end
            %---多余电量判断结束---%
            
            %---缺电判断---%
            if  SOC(i+1,1)<=SOCmin
                SOC(i+1,1)=SOCmin;
                %计算即将见底时的蓄电池真实放电功率
                if  SOC(i,1)~=SOCmin
                    Pbss(i,1)=(SOC(i+1,1)-(1-0.03)*SOC(i,1))*Eb*(-0.9)*X(j,3);
                else
                    Pbss(i,1)=0;
                end
                %在见底时的真实缺电率
                if elecnum(i,1)>=abs(Pbss(i,1))
                    elecnum(i,1)=elecnum(i,1)-abs(Pbss(i,1));                       %这个设置的很巧妙
                else
                    elecnum(i,1)=0;
                end
            else
                elecnum(i,1)=0;
            end
            
            %外购电耗费
            Celebuy(i,1)=0;
            
            %售电盈利
            Celesell(i,1)=More(i,1)*Cse;
            %Emore(i,1)=More(i,1)/(Ppv(i,1)*X(j,1)+Pwt(i,1)*X(j,2)+0.001);
            
            %---缺电判断结束---%
            %%%---电模块结束---%%%
        end
        %%%---18.00-23.00计算结束---%%%
        
        %%%---此季节计算结束---%%%
    end
    
    
    %目标函数计算
    %Chengfa(j,1)=sum(LPS)*0.5*20;                                                          %1000:惩罚系数k；   20：电费单价；
    Xgem=max(Xge);
    LCC1(j,1)=Cpv*X(j,1)+Cwt*X(j,2)+Cbss*X(j,3)+Chp*X(j,4)+Cge*Xgem+Ctes*X(j,5);            %系统购置
    Cs(j,1)=(sum(Celebuy)-sum(Celesell))*60;                                                %购电与售电
    LCC(j,1)=LCC1(j,1)+Cs(j,1);                                                             %目标函数1：全寿命周期系统经济成本
    CO2ot(j,1)=COpv*X(j,1)+COwt*X(j,2)+CObss*X(j,3)+COhp*X(j,4)+COge*Xgem+COtes*X(j,5);     %初始碳排放
    CO2op(j,1)=(sum(elecnum)+sum(Phptes))*0.89*60;                                          %运行碳排放
    CO2(j,1)=CO2ot(j,1)+CO2op(j,1);                                                         %目标函数2：全寿命周期碳排放
    EER(j,1)=(sum(Pwt)*X(j,2)+sum(Ppv)*X(j,1)-sum(More)+sum(Qhpscpl))/(sum(Pwt)*X(j,2)+sum(Ppv)*X(j,1)-sum(More)+sum(elecnum)+sum(Phptes)+sum(HeatingFuhe)+sum(CoolingFuhe));                    %目标函数3：外电比
    %Qloss(j,1)=sum(Xloss)+sum(More)/(sum(Pwt)*X(j,2)+sum(Ppv)*X(j,1));                                                       %约束条件，Xloss(j,1)=0
    %Qloss(j,1)=sum(Xloss)+sum(More);
    Qloss(j,1)=sum(Qdeloss)+sum(More);
    %Qloss(j,1)=sum(More);
    
    %LPSP(j,1)=sum(LPS)/(sum(Fuhe)+sum(Php)*X(j,4));                                         %缺电率
    %SPSP(j,1)=sum(More)/(sum(Fuhe)+sum(Php)*X(j,4));                                        %浪费电率
    %SHSP(j,1)=(sum(MoreHeating)+sum(MoreCooling))/(sum(CoolingFuhe)+sum(HeatingFuhe));      %浪费热/冷率
    %SHSP(j,1)=sum(MoreHeating)/sum(HeatingFuhe)+sum(MoreCooling)/sum(CoolingFuhe);
    %LHSP(j,1)=(sum(LossHeating)+sum(LossCooling))/(sum(CoolingFuhe)+sum(HeatingFuhe));      %缺热/冷率
    %LHSP(j,1)=sum(LossHeating)/sum(HeatingFuhe)+sum(LossCooling)/sum(CoolingFuhe);
    %LGWMR(j,1)=sum(LossGeHeating)/(sum(Qhpscc)*X(j,4));                                     %地热缺失率
    %MGWMR(j,1)=sum(MoreGeHeating)/(sum(Qhpscc)*X(j,4));                                     %地热过载率
    %PopObj(j,1) = LCC(j,1)/100000000*2/3;
    PopObj(j,1) = LCC(j,1);
    PopObj(j,2) = CO2(j,1);
    PopObj(j,3) = 1-EER(j,1);
    PopObj(j,4) = Qloss(j,1);
end
    %%--单行粒子解集计算结束--%%
    

    
end


%F2=0.45*(LPSP+SPSP)+0.2*LHSP+0.2*SHSP+0.15*LSGSP;                                 %目标函数2：系统可靠性目标
%F2=0.35*LHSP+0.35*SHSP+0.3*LSGSP;
%F2=0.4*LPSP+0.3*LHSP+0.3*LGWMR;         


end
