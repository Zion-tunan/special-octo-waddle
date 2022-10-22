%% Ŀ�꺯����Լ������
%% �����д��л־Զ

%% ����Test_function���������λ����Ϣ
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
        
%% ��Ŀ��Ŀ�꺯��
    case 'UF10'
%% ������������

%������������
global Ppv Pwt Qhpsh Qhpsc Php

%���ظ����ļ�
load CoolingFuhe.txt
load HeatingFuhe.txt
load Fuhe.txt
load Temperature.txt
load DianJia.txt

%�����豸��������
SOCmax=0.9;
SOCmin=0.1;
Eb=50;                                                       %���ض����
X=fix(X);
SOC(1,1)=0.9;
etaeh=0.95;
Qge=6;
Cse=0;

%ȫ�������ھ��óɱ���60�꣬GHIES��yuan��
Cpv=600*2.4+2.54*60;                                           %�����سɱ�  ����ɱ�+����ά���ɱ�
Cwt=12720*30*3+65*30*60;                                       %��������ɱ�
Cbss=100000*5+15.9*50*60;                                      %���سɱ�
Chp=240000*4+1200*60;                                          %�ȱû���ɱ�
Cge=20000;                                                     %���ȴ򾮳ɱ� 
Ctes=2360*3;                                                     %����/��ˮ��ɱ�

%ȫ����̼�ŷųɱ���60�꣬kg��
COpv=1212.90;                                               %�����سɱ�  ����ɱ�+����ά��̼�ŷųɱ�
COwt=28695.55;                                              %��������̼�ŷųɱ�
CObss=51450;                                                %����̼�ŷųɱ�
COhp=355439.5;                                              %�ȱû���̼�ŷųɱ�
COge=7785.5;                                                %���ȴ�̼�ŷųɱ� 
COtes=2738.65;                                               %����ˮ��̼�ŷ�


Ttessdh=53;                                          %ˮ�乩ůʱ��ͣ�¶�
Ttessdhh=50;
Ttessdc=2;
Ttessdcc=5;
Vtes=2.36;                                              %ˮ�����
Ates=10.99;
Ktes=0.0012;
mtes=Vtes*997;
Cw=4.18;
Tbui=45;
Tbuc=10;
 
%Ԥ����
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


%% ����⼯��ѭ��

for j= 1:200
    
    %%--�������ӽ⼯����--%%
    
    %% ��ů��part01, 1-2520h
    
    for i=1:2520
        
        a(i,1)=i;
        b(i,1)=a(i,1)/24;
        c(i,1)=fix(b(i,1));
        d(i,1)=b(i,1)-c(i,1);
        
        %%%---23.00-07.00,���ʱ�䣬�ȱ���ˮ�乩��---%%%
        if d(i,1)<0.31
            
            if Ttes(i,1)<Ttessdhh                                                                                            %ˮ���¶ȿ�����ͣ�жϣ�ˮ���¶�<53
                Qhptes(i,1)=((Ttessdh-Ttes(i,1))*Cw*X(j,5)*mtes+X(j,5)*Ktes*Ates*(Ttes(i,1)-Temperature(i,1)))/3600;                %�ȱ÷��ʱ�ڹ����ȸ���
                
                %%---�ȱò��ָ��ɼ�ȱʧ����---%%
                
                %��һ�μ����ȱ�����
                Xhp01(i,1)=Qhptes(i,1)/Qhpsh(i,1);
                
                if Xhp01(i,1)<1
                    Xhp02(i,1)=0;
                else
                    Xhp02(i,1)=fix(Xhp01(i,1))+1;
                end
                
                Xhp(i,1)=Xhp02(i,1);
                
                %�ڶ��μ����ȱ����������ݱ߽�����,������ȱʧ���벿�ָ����µ��ȱù��ġ�COP
                if Xhp(i,1)<=X(j,4) && Xhp(i,1)>0
                    PL(i,1)=Qhptes(i,1)/(Qhpsh(i,1)*Xhp02(i,1));                          %���1�����ڲ��ַ����µ��ȱ�
                    COPnom(i,1)=Qhpsh(i,1)/Php(i,1);
                    COPpl(i,1)=COPnom(i,1)*(-0.00006*PL(i,1)^6+0.0017*PL(i,1)^5-0.0181*PL(i,1)^4+0.096*PL(i,1)^3-0.2697*PL(i,1)^2+0.4276*PL(i,1)+0.7626);
                    Phppl(i,1)=Qhptes(i,1)/COPpl(i,1);                                    %���㲿�ָ����¸��ȱõĹ���
                    Qhpscpl(i,1)=Qhptes(i,1)-Phppl(i,1);
                    Ttes(i+1,1)=Ttessdh;
                elseif Xhp(i,1)>X(j,4)
                    Xhp02(i,1)=fix(X(j,4));
                    %Qdeloss(i,1)=Qhptes(i,1)-Qhpsh(i,1)*Xhp02(i,1);                          %Qdeloss����ȱʧ������
                    %Xloss(i,1)=Qdeloss(i,1)/Qhptes(i,1);
                    Phppl(i,1)=Php(i,1)*Xhp02(i,1);
                    Qhpscpl(i,1)=Xhp02(i,1)*Qhpsc(i,1);
                    Ttes(i+1,1)=(Xhp02(i,1)*Qhpsh(i,1)*3600-X(j,5)*Ktes*Ates*(Ttes(i,1)-Temperature(i,1)))/(mtes*Cw*X(j,5))+Ttes(i,1);
                else
                    Phppl(i,1)=0;
                    Qhpscpl(i,1)=0;
                    Ttes(i+1,1)=Ktes*Ates*(Temperature(i,1)-Ttes(i,1))/(mtes*Cw)+Ttes(i,1);
                end                                                                      %Xhp��ʾ���յ��ȱ���������С���㣬���2 X(j,4)��ʾ���������߽�ֵ
                
                Phptes(i,1)=Phppl(i,1);
                
                %������������
                Xge01(i,1)=Qhpscpl(i,1)/(6*0.95);                                        %����ܴ���Ч������Ϊ95%
                Xge(i,1)=fix(Xge01(i,1));
                if Xge01(i,1)-Xge(i,1) > 0                                                 %�Ե��������ȡ��������+1
                    Xge(i,1)=fix(Xge01(i,1))+1;
                else
                    Xge(i,1)=Xge01(i,1);
                end
                
            else
                Ttes(i+1,1)=Ktes*Ates*(Temperature(i,1)-Ttes(i,1))/(mtes*Cw)+Ttes(i,1);
            end
            
            
            %%%---�縺�ɼ���---%%%
            Pbss(i,1)=Fuhe(i,1)-Ppv(i,1)*X(j,1)-Pwt(i,1)*X(j,2);
            
            %---���س�ŵ�״̬�ж�---%
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
            %---���س�ŵ�״̬�жϽ���---%
            
            %---��������ж�---%
            if  SOC(i+1,1)>=SOCmax
                SOC(i+1,1)=SOCmax;
                %�ȼ��㼴�����ʱ��������ʵ��繦��
                if  SOC(i,1)~=SOCmax
                    Pbss(i,1)=(SOC(i+1,1)-(1-0.03)*SOC(i,1))*Eb*X(j,3)/(-0.9);
                else
                    Pbss(i,1)=0;
                end
                %��¼�������
                %��100%�������Ϊ����������
                if More(i,1)>=abs(Pbss(i,1))
                    More(i,1)=More(i,1)-abs(Pbss(i,1));
                else
                    More(i,1)=0;
                end
                %��¼�������
            else
                More(i,1)=0;
            end
            %---��������жϽ���---%
            
            %---ȱ���ж�---%
            if  SOC(i+1,1)<=SOCmin
                SOC(i+1,1)=SOCmin;
                %���㼴������ʱ��������ʵ�ŵ繦��
                if  SOC(i,1)~=SOCmin
                    Pbss(i,1)=(SOC(i+1,1)-(1-0.03)*SOC(i,1))*Eb*(-0.9)*X(j,3);
                else
                    Pbss(i,1)=0;
                end
                %�ڼ���ʱ����ʵȱ����
                if elecnum(i,1)>=abs(Pbss(i,1))
                    elecnum(i,1)=elecnum(i,1)-abs(Pbss(i,1));
                else
                    elecnum(i,1)=0;
                end
            else
                elecnum(i,1)=0;
            end
            
            %�⹺��ķ�
            Celebuy(i,1)=Phptes(i,1)*DianJia(i,1);
            
            %�۵�ӯ��
            Celesell(i,1)=More(i,1)*Cse;
            %Emore(i,1)=More(i,1)/(Ppv(i,1)*X(j,1)+Pwt(i,1)*X(j,2)+0.001);
            
            %---ȱ���жϽ���---%
            %%%---��ģ�����---%%%
        end
        
        %%%---23.00-07.00�������---%%%
        
        %%%---07.00-08.00���㣬���ϰ�ʱ�䣬�Ƿ��ʱ��---%%%
        if d(i,1)>0.31 && d(i,1)<0.36
            
            Ttes(i+1,1)=Ktes*Ates*(Temperature(i,1)-Ttes(i,1))/(mtes*Cw)+Ttes(i,1);
            
            %%%---�縺�ɼ���---%%%
            Pbss(i,1)=Fuhe(i,1)-Ppv(i,1)*X(j,1)-Pwt(i,1)*X(j,2);
            
            %---���س�ŵ�״̬�ж�---%
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
            %---���س�ŵ�״̬�жϽ���---%
            
            %---��������ж�---%
            if  SOC(i+1,1)>=SOCmax
                SOC(i+1,1)=SOCmax;
                %�ȼ��㼴�����ʱ��������ʵ��繦��
                if  SOC(i,1)~=SOCmax
                    Pbss(i,1)=(SOC(i+1,1)-(1-0.03)*SOC(i,1))*Eb*X(j,3)/(-0.9);
                else
                    Pbss(i,1)=0;
                end
                %��¼�������
                %��100%�������Ϊ����������
                if More(i,1)>=abs(Pbss(i,1))
                    More(i,1)=More(i,1)-abs(Pbss(i,1));
                else
                    More(i,1)=0;
                end
                %��¼�������
            else
                More(i,1)=0;
            end
            %---��������жϽ���---%
            
            %---ȱ���ж�---%
            if  SOC(i+1,1)<=SOCmin
                SOC(i+1,1)=SOCmin;
                %���㼴������ʱ��������ʵ�ŵ繦��
                if  SOC(i,1)~=SOCmin
                    Pbss(i,1)=(SOC(i+1,1)-(1-0.03)*SOC(i,1))*Eb*(-0.9)*X(j,3);
                else
                    Pbss(i,1)=0;
                end
                %�ڼ���ʱ����ʵȱ����
                if elecnum(i,1)>=abs(Pbss(i,1))
                    elecnum(i,1)=elecnum(i,1)-abs(Pbss(i,1));
                else
                    elecnum(i,1)=0;
                end
            else
                elecnum(i,1)=0;
            end
            
            %�⹺��ķ�
            Celebuy(i,1)=0;
            
            %�۵�ӯ��
            Celesell(i,1)=More(i,1)*Cse;
            %Emore(i,1)=More(i,1)/(Ppv(i,1)*X(j,1)+Pwt(i,1)*X(j,2)+0.001);
            
        end
        %%%---07.00-08.00�������---%%%
        
        %%%---8.00-9.00,����ʱ�䣬ƽ�����---%%%
        if d(i,1)>0.36 && d(i,1)<0.39
            
            %%--ˮ���¶ȼ���--%%
            Ttes(i+1,1)=Ktes*Ates*(Temperature(i,1)-Ttes(i,1))/(mtes*Cw)+Ttes(i,1);
            %%--ˮ���¶ȼ������--%%
            
            %%--�ȱü���--%%
            %��һ�μ����ȱ�����
            Xhp01(i,1)=HeatingFuhe(i,1)/Qhpsh(i,1);
            
            if Xhp01(i,1)<1
                Xhp02(i,1)=0;
            else
                Xhp02(i,1)=fix(Xhp01(i,1))+1;
            end
            
            Xhp(i,1)=Xhp02(i,1);
                
            %�ڶ��μ����ȱ����������ݱ߽�����,������ȱʧ���벿�ָ����µ��ȱù��ġ�COP
            if Xhp(i,1)<=X(j,4) && Xhp(i,1)>0
                PL(i,1)=HeatingFuhe(i,1)/(Qhpsh(i,1)*Xhp02(i,1));                          %���1�����ڲ��ַ����µ��ȱ�
                COPnom(i,1)=Qhpsh(i,1)/Php(i,1);
                COPpl(i,1)=COPnom(i,1)*(-0.00006*PL(i,1)^6+0.0017*PL(i,1)^5-0.0181*PL(i,1)^4+0.096*PL(i,1)^3-0.2697*PL(i,1)^2+0.4276*PL(i,1)+0.7626);
                Phppl(i,1)=HeatingFuhe(i,1)/COPpl(i,1);
                Qhpscpl(i,1)=HeatingFuhe(i,1)-Phppl(i,1);
            elseif Xhp(i,1)>X(j,4)
                Xhp02(i,1)=fix(X(j,4));
                Qdeloss(i,1)=HeatingFuhe(i,1)-Qhpsh(i,1)*Xhp02(i,1);                          %Qdeloss����ȱʧ������
                Xloss(i,1)=Qdeloss(i,1)/HeatingFuhe(i,1);
                Phppl(i,1)=Php(i,1)*Xhp02(i,1);
                Qhpscpl(i,1)=Xhp02(i,1)*Qhpsc(i,1);
            else
                Qdeloss(i,1)=HeatingFuhe(i,1);                          %Qdeloss����ȱʧ������
                Xloss(i,1)=1;
                Phppl(i,1)=0;
                Qhpscpl(i,1)=0;
            end                                                                      %Xhp��ʾ���յ��ȱ���������С���㣬���2 X(j,4)��ʾ���������߽�ֵ
            
            %�������ܵ�����
            Xge01(i,1)=Qhpscpl(i,1)/(6*0.95);
            Xge(i,1)=fix(Xge(i,1));
            
            if Xge01(i,1)-Xge(i,1) > 0                                                %�Ե��������ȡ��������+1
                Xge(i,1)=fix(Xge01(i,1))+1;
            else
                Xge(i,1)=Xge01(i,1);
            end
            %%--�ȱü������--%%
            
            %%--�縺�ɼ���--%%
            Pbss(i,1)=Fuhe(i,1)+Phppl(i,1)-Ppv(i,1)*X(j,1)-Pwt(i,1)*X(j,2);
            
            %���س�ŵ�״̬�ж�
            if Pbss(i,1) < 0                                                         %pbssΪ��������Ҫ��磬��ʱ��pbss�������������޵Ĺ���
                SOC(i+1,1)=(1-0.03)*SOC(i,1)-Pbss(i,1)*0.9/(X(j,3)*Eb);
                More(i,1)=abs(Pbss(i,1));                                            %��¼����ĵ���
                elecnum(i,1)=0;
            elseif Pbss(i,1) > 0                                                     %pbssΪ��������Ҫ�ŵ�
                SOC(i+1,1)=(1-0.03)*SOC(i,1)-Pbss(i,1)/(X(j,3)*0.9*Eb);
                More(i,1)=0;                                                         %��¼������˷ѵ���
                elecnum(i,1)=abs(Pbss(i,1));
            else
                SOC(i+1,1)=SOC(i,1);
                More(i,1)=0;                                                         %��¼������˷ѵ���
                elecnum(i,1)=0;
            end
            
            %��������ж�
            if SOC(i+1,1)>=SOCmax
                SOC(i+1,1)=SOCmax;
                %�ȼ��㼴�����ʱ��������ʵ��繦��
                if SOC(i,1)~=SOCmax                                                      %~=�ǲ����ڵ���˼
                    Pbss(i,1)=(SOC(i+1,1)-(1-0.03)*SOC(i,1))*Eb*X(j,3)/(-0.9);
                else
                    Pbss(i,1)=0;
                end
                %��100%�������Ϊ����������
                if More(i,1)>=abs(Pbss(i,1))
                    More(i,1)=More(i,1)-abs(Pbss(i,1));
                else
                    More(i,1)=0;
                end
            else
                More(i,1)=0;
            end
            
            %ȱ���ж�
            if SOC(i+1,1)<=SOCmin
                SOC(i+1,1)=SOCmin;
                %���㼴������ʱ��������ʵ�ŵ繦��
                if SOC(i,1)~=SOCmin
                    Pbss(i,1)=(SOC(i+1,1)-(1-0.03)*SOC(i,1))*Eb*(-0.9)*X(j,3);
                else
                    Pbss(i,1)=0;
                end
                %��ʵȱ����
                if elecnum(i,1)>=abs(Pbss(i,1))
                    elecnum(i,1)=elecnum(i,1)-abs(Pbss(i,1));
                else
                    elecnum(i,1)=0;
                end
            else
                elecnum(i,1)=0;
            end
            
            %�⹺��ķ�
            Celebuy(i,1)=elecnum(i,1)*DianJia(i,1);
            
            %�۵�ӯ��
            Celesell(i,1)=More(i,1)*Cse;
            %Emore(i,1)=More(i,1)/(Ppv(i,1)*X(j,1)+Pwt(i,1)*X(j,2)+0.001);
            %%--�縺�ɼ������--%%
            
        end
        %%%---8.00-9.00�������---%%%
        
        %%%---9.00-12.00,����ʱ�䣬��ֵ��ۣ�����ˮ��---%%%
        if d(i,1)>0.39 && d(i,1)<0.51
            
            %%--�ȱü���--%%
            Qtes(i,1)=((Ttes(i,1)-Tbui)*Cw*X(j,5)*mtes-X(j,5)*Ktes*Ates*(Ttes(i,1)-Temperature(i,1)))/3600;
            
            %��ֹ���������ȱø���
            if Qtes(i,1)<= 0
                Qtes(i,1)=0;
                
                %%--ˮ���¶ȼ���--%%
                Ttes(i+1,1)=Ktes*Ates*(Temperature(i,1)-Ttes(i,1))/(mtes*Cw)+Ttes(i,1);
                %%--ˮ���¶ȼ������--%%
            
                %%--�ȱü���--%%
                %��һ�μ����ȱ�����
                Xhp01(i,1)=HeatingFuhe(i,1)/Qhpsh(i,1);
                
                if Xhp01(i,1)<1
                    Xhp02(i,1)=0;
                else
                    Xhp02(i,1)=fix(Xhp01(i,1))+1;
                end
                
                Xhp(i,1)=Xhp02(i,1);
                
                %�ڶ��μ����ȱ����������ݱ߽�����,������ȱʧ���벿�ָ����µ��ȱù��ġ�COP
                if Xhp(i,1)<=X(j,4) && Xhp(i,1)>0
                    PL(i,1)=HeatingFuhe(i,1)/(Qhpsh(i,1)*Xhp02(i,1));                          %���1�����ڲ��ַ����µ��ȱ�
                    COPnom(i,1)=Qhpsh(i,1)/Php(i,1);
                    COPpl(i,1)=COPnom(i,1)*(-0.00006*PL(i,1)^6+0.0017*PL(i,1)^5-0.0181*PL(i,1)^4+0.096*PL(i,1)^3-0.2697*PL(i,1)^2+0.4276*PL(i,1)+0.7626);
                    Phppl(i,1)=HeatingFuhe(i,1)/COPpl(i,1);
                    Qhpscpl(i,1)=HeatingFuhe(i,1)-Phppl(i,1);
                elseif Xhp(i,1)>X(j,4)
                    Xhp02(i,1)=fix(X(j,4));
                    Qdeloss(i,1)=HeatingFuhe(i,1)-Qhpsh(i,1)*Xhp02(i,1);                          %Qdeloss����ȱʧ������
                    Xloss(i,1)=Qdeloss(i,1)/HeatingFuhe(i,1);
                    Phppl(i,1)=Php(i,1)*Xhp02(i,1);
                    Qhpscpl(i,1)=Xhp02(i,1)*Qhpsc(i,1);
                else
                    Qdeloss(i,1)=HeatingFuhe(i,1);                          %Qdeloss����ȱʧ������
                    Xloss(i,1)=1;
                    Phppl(i,1)=0;
                    Qhpscpl(i,1)=0;
                end                                                                      %Xhp��ʾ���յ��ȱ���������С���㣬���2 X(j,4)��ʾ���������߽�ֵ
            
                %�������ܵ�����
                Xge01(i,1)=Qhpscpl(i,1)/(6*0.95);
                Xge(i,1)=fix(Xge(i,1));
            
                if Xge01(i,1)-Xge(i,1) > 0                                                %�Ե��������ȡ��������+1
                    Xge(i,1)=fix(Xge01(i,1))+1;
                else
                    Xge(i,1)=Xge01(i,1);
                end
                %%--�ȱü������--%%
            end
            
            %���1/2/3
            if Qtes(i,1) >= HeatingFuhe(i,1)
                Ttes(i+1,1)=(X(j,5)*Ktes*Ates*(Temperature(i,1)-Ttes(i,1))-HeatingFuhe(i,1)*3600)/(X(j,5)*mtes*Cw)+Ttes(i,1);
                Xhp02(i,1)=0;
            elseif Qtes(i,1)<HeatingFuhe(i,1) && Qtes(i,1)>0
                Qhpbui(i,1)=HeatingFuhe(i,1)-Qtes(i,1);
                Ttes(i+1,1)=Tbui;
                
                %���2/3
                %��һ�μ����ȱ�����
                Xhp01(i,1)=Qhpbui(i,1)/Qhpsh(i,1);
                
                if Xhp01(i,1)<1
                    Xhp02(i,1)=0;
                else
                    Xhp02(i,1)=fix(Xhp01(i,1))+1;
                end
                
                Xhp(i,1)=Xhp02(i,1);
                
                %�ڶ��μ����ȱ����������ݱ߽�����,������ȱʧ���벿�ָ����µ��ȱù��ġ�COP
                if Xhp(i,1) <= X(j,4) && Xhp(i,1)>0
                    PL(i,1)=Qhpbui(i,1)/(Qhpsh(i,1)*Xhp02(i,1));                          %���1�����ڲ��ַ����µ��ȱ�
                    COPnom(i,1)=Qhpsh(i,1)/Php(i,1);
                    COPpl(i,1)=COPnom(i,1)*(-0.00006*PL(i,1)^6+0.0017*PL(i,1)^5-0.0181*PL(i,1)^4+0.096*PL(i,1)^3-0.2697*PL(i,1)^2+0.4276*PL(i,1)+0.7626);
                    Phppl(i,1)=Qhpbui(i,1)/COPpl(i,1);
                    Qhpscpl(i,1)=Qhpbui(i,1)-Phppl(i,1);
                elseif Xhp(i,1)> X(j,4)
                    Xhp02(i,1)=fix(X(j,4));
                    Qdeloss(i,1)=Qhpbui(i,1)-Qhpsh(i,1)*Xhp02(i,1);                          %Qdeloss����ȱʧ������
                    Xloss(i,1)=Qdeloss(i,1)/Qhpbui(i,1);
                    Phppl(i,1)=Php(i,1)*Xhp02(i,1);
                    Qhpscpl(i,1)=Xhp02(i,1)*Qhpsc(i,1);
                else
                    Qdeloss(i,1)=Qhpbui(i,1);                          %Qdeloss����ȱʧ������
                    Xloss(i,1)=1;
                    Phppl(i,1)=0;
                    Qhpscpl(i,1)=0;
                end                                                                      %Xhp��ʾ���յ��ȱ���������С���㣬���2 X(j,4)��ʾ���������߽�ֵ
                
                %������������
                Xge01(i,1)=Qhpscpl(i,1)/(6*0.95);
                Xge(i,1)=fix(Xge(i,1));
                
                if Xge01(i,1)-Xge(i,1) > 0                                               %�Ե��������ȡ��������+1
                    Xge(i,1)=fix(Xge01(i,1))+1;
                else
                    Xge(i,1)=Xge01(i,1);
                end
                
            end
            %%--�ȱü������--%%
            
            %%--�縺�ɼ���--%%
            Pbss(i,1)=Fuhe(i,1)+Phppl(i,1)-Ppv(i,1)*X(j,1)-Pwt(i,1)*X(j,2);
            
            %���س�ŵ�״̬�ж�
            if Pbss(i,1) < 0                                                            %pbssΪ��������Ҫ��磬��ʱ��pbss�������������޵Ĺ���
                SOC(i+1,1)=(1-0.03)*SOC(i,1)-Pbss(i,1)*0.9/(X(j,3)*Eb);
                More(i,1)=abs(Pbss(i,1));                                               %��¼����ĵ���
                elecnum(i,1)=0;
            elseif Pbss(i,1) > 0                                                        %pbssΪ��������Ҫ�ŵ�
                SOC(i+1,1)=(1-0.03)*SOC(i,1)-Pbss(i,1)/(X(j,3)*0.9*Eb);
                More(i,1)=0;                                                            %��¼������˷ѵ���
                elecnum(i,1)=abs(Pbss(i,1));
            else
                SOC(i+1,1)=SOC(i,1);
                More(i,1)=0;                                                            %��¼������˷ѵ���
                elecnum(i,1)=0;
            end
            
            %��������ж�
            if SOC(i+1,1)>=SOCmax
                SOC(i+1,1)=SOCmax;
                %�ȼ��㼴�����ʱ��������ʵ��繦��
                if SOC(i,1)~=SOCmax                                                         %~=�ǲ����ڵ���˼
                    Pbss(i,1)=(SOC(i+1,1)-(1-0.03)*SOC(i,1))*Eb*X(j,3)/(-0.9);
                else
                    Pbss(i,1)=0;
                end
                %��100%�������Ϊ����������
                if More(i,1)>=abs(Pbss(i,1))
                    More(i,1)=More(i,1)-abs(Pbss(i,1));
                else
                    More(i,1)=0;
                end
                %��¼�������
            else
                More(i,1)=0;
            end
            
            %ȱ���ж�
            if SOC(i+1,1)<=SOCmin
                SOC(i+1,1)=SOCmin;
                %���㼴������ʱ��������ʵ�ŵ繦��
                if SOC(i,1)~=SOCmin
                    Pbss(i,1)=(SOC(i+1,1)-(1-0.03)*SOC(i,1))*Eb*(-0.9)*X(j,3);
                else
                    Pbss(i,1)=0;
                end
                %��ʵȱ����
                if elecnum(i,1)>=abs(Pbss(i,1))
                    elecnum(i,1)=elecnum(i,1)-abs(Pbss(i,1));
                else
                    elecnum(i,1)=0;
                end
            else
                elecnum(i,1)=0;
            end
            
            %�⹺��ķ�
            Celebuy(i,1)=elecnum(i,1)*DianJia(i,1);
            
            %�۵�ӯ��
            Celesell(i,1)=More(i,1)*Cse;
            %Emore(i,1)=More(i,1)/(Ppv(i,1)*X(j,1)+Pwt(i,1)*X(j,2)+0.001);
            %%--�縺�ɼ������--%%
            
        end
        %%%---9.00-12.00�������---%%%
        
        %%%---12.00-16.00������ʱ�䣬ƽ�����---%%%
        if d(i,1)>0.51 && d(i,1)<0.69
            
            %%--ˮ���¶ȼ���--%%
            Ttes(i+1,1)=Ktes*Ates*(Temperature(i,1)-Ttes(i,1))/(mtes*Cw)+Ttes(i,1);
            %%--ˮ���¶ȼ������--%%
            
            %%--�ȱü���--%%
            %��һ�μ����ȱ�����
            Xhp01(i,1)=HeatingFuhe(i,1)/Qhpsh(i,1);
            
            if Xhp01(i,1)<1
                Xhp02(i,1)=0;
            else
                Xhp02(i,1)=fix(Xhp01(i,1))+1;
            end
            
            Xhp(i,1)=Xhp02(i,1);
                
            %�ڶ��μ����ȱ����������ݱ߽�����,������ȱʧ���벿�ָ����µ��ȱù��ġ�COP
            if Xhp(i,1)<=X(j,4) && Xhp(i,1)>0
                PL(i,1)=HeatingFuhe(i,1)/(Qhpsh(i,1)*Xhp02(i,1));                          %���1�����ڲ��ַ����µ��ȱ�
                COPnom(i,1)=Qhpsh(i,1)/Php(i,1);
                COPpl(i,1)=COPnom(i,1)*(-0.00006*PL(i,1)^6+0.0017*PL(i,1)^5-0.0181*PL(i,1)^4+0.096*PL(i,1)^3-0.2697*PL(i,1)^2+0.4276*PL(i,1)+0.7626);
                Phppl(i,1)=HeatingFuhe(i,1)/COPpl(i,1);
                Qhpscpl(i,1)=HeatingFuhe(i,1)-Phppl(i,1);
            elseif Xhp(i,1)>X(j,4)
                Xhp02(i,1)=fix(X(j,4));
                Qdeloss(i,1)=HeatingFuhe(i,1)-Qhpsh(i,1)*Xhp02(i,1);                          %Qdeloss����ȱʧ������
                Xloss(i,1)=Qdeloss(i,1)/HeatingFuhe(i,1);
                Phppl(i,1)=Php(i,1)*Xhp02(i,1);
                Qhpscpl(i,1)=Xhp02(i,1)*Qhpsc(i,1);
            else
                Qdeloss(i,1)=HeatingFuhe(i,1);                          %Qdeloss����ȱʧ������
                Xloss(i,1)=1;
                Phppl(i,1)=0;
                Qhpscpl(i,1)=0;
            end                                                                      %Xhp��ʾ���յ��ȱ���������С���㣬���2 X(j,4)��ʾ���������߽�ֵ
            
            %�������ܵ�����
            Xge01(i,1)=Qhpscpl(i,1)/(6*0.95);
            Xge(i,1)=fix(Xge(i,1));
            
            if Xge01(i,1)-Xge(i,1) > 0                                                %�Ե��������ȡ��������+1
                Xge(i,1)=fix(Xge01(i,1))+1;
            else
                Xge(i,1)=Xge01(i,1);
            end
            %%--�ȱü������--%%
            
            %%--�縺�ɼ���--%%
            Pbss(i,1)=Fuhe(i,1)+Phppl(i,1)-Ppv(i,1)*X(j,1)-Pwt(i,1)*X(j,2);
            
            %���س�ŵ�״̬�ж�
            if Pbss(i,1) < 0                                                         %pbssΪ��������Ҫ��磬��ʱ��pbss�������������޵Ĺ���
                SOC(i+1,1)=(1-0.03)*SOC(i,1)-Pbss(i,1)*0.9/(X(j,3)*Eb);
                More(i,1)=abs(Pbss(i,1));                                            %��¼����ĵ���
                elecnum(i,1)=0;
            elseif Pbss(i,1) > 0                                                     %pbssΪ��������Ҫ�ŵ�
                SOC(i+1,1)=(1-0.03)*SOC(i,1)-Pbss(i,1)/(X(j,3)*0.9*Eb);
                More(i,1)=0;                                                         %��¼������˷ѵ���
                elecnum(i,1)=abs(Pbss(i,1));
            else
                SOC(i+1,1)=SOC(i,1);
                More(i,1)=0;                                                         %��¼������˷ѵ���
                elecnum(i,1)=0;
            end
            
            %��������ж�
            if SOC(i+1,1)>=SOCmax
                SOC(i+1,1)=SOCmax;
                %�ȼ��㼴�����ʱ��������ʵ��繦��
                if SOC(i,1)~=SOCmax                                                      %~=�ǲ����ڵ���˼
                    Pbss(i,1)=(SOC(i+1,1)-(1-0.03)*SOC(i,1))*Eb*X(j,3)/(-0.9);
                else
                    Pbss(i,1)=0;
                end
                %��100%�������Ϊ����������
                if More(i,1)>=abs(Pbss(i,1))
                    More(i,1)=More(i,1)-abs(Pbss(i,1));
                else
                    More(i,1)=0;
                end
                %��¼�������
            else
                More(i,1)=0;
            end
            
            %ȱ���ж�
            if SOC(i+1,1)<=SOCmin
                SOC(i+1,1)=SOCmin;
                %���㼴������ʱ��������ʵ�ŵ繦��
                if SOC(i,1)~=SOCmin
                    Pbss(i,1)=(SOC(i+1,1)-(1-0.03)*SOC(i,1))*Eb*(-0.9)*X(j,3);
                else
                    Pbss(i,1)=0;
                end
                %��ʵȱ����
                if elecnum(i,1)>=abs(Pbss(i,1))
                    elecnum(i,1)=elecnum(i,1)-abs(Pbss(i,1));
                else
                    elecnum(i,1)=0;
                end
            else
                elecnum(i,1)=0;
            end
            
            %�⹺��ķ�
            Celebuy(i,1)=elecnum(i,1)*DianJia(i,1);
            
            %�۵�ӯ��
            Celesell(i,1)=More(i,1)*Cse;
            %Emore(i,1)=More(i,1)/(Ppv(i,1)*X(j,1)+Pwt(i,1)*X(j,2)+0.001);
            %%--�縺�ɼ������--%%
            
        end
        %%%---12.00-16.00�������---%%%
        
        %%%---16.00-18.00,����ʱ�䣬��ֵ��ۣ�����ˮ��---%%%
        if d(i,1)>0.69 && d(i,1)<0.78
            
            %%--�ȱü���--%%
            Qtes(i,1)=((Ttes(i,1)-Tbui)*Cw*X(j,5)*mtes-X(j,5)*Ktes*Ates*(Ttes(i,1)-Temperature(i,1)))/3600;
            
            %��ֹ���������ȱø���
            if Qtes(i,1)<0
                Qtes(i,1)=0;
                
                %%--ˮ���¶ȼ���--%%
                Ttes(i+1,1)=Ktes*Ates*(Temperature(i,1)-Ttes(i,1))/(mtes*Cw)+Ttes(i,1);
                %%--ˮ���¶ȼ������--%%
            
                %%--�ȱü���--%%
                %��һ�μ����ȱ�����
                Xhp01(i,1)=HeatingFuhe(i,1)/Qhpsh(i,1);
                
                if Xhp01(i,1)<1
                    Xhp02(i,1)=0;
                else
                    Xhp02(i,1)=fix(Xhp01(i,1))+1;
                end
                
                Xhp(i,1)=Xhp02(i,1);
                
                %�ڶ��μ����ȱ����������ݱ߽�����,������ȱʧ���벿�ָ����µ��ȱù��ġ�COP
                if Xhp(i,1)<=X(j,4) && Xhp(i,1)>0
                    PL(i,1)=HeatingFuhe(i,1)/(Qhpsh(i,1)*Xhp02(i,1));                          %���1�����ڲ��ַ����µ��ȱ�
                    COPnom(i,1)=Qhpsh(i,1)/Php(i,1);
                    COPpl(i,1)=COPnom(i,1)*(-0.00006*PL(i,1)^6+0.0017*PL(i,1)^5-0.0181*PL(i,1)^4+0.096*PL(i,1)^3-0.2697*PL(i,1)^2+0.4276*PL(i,1)+0.7626);
                    Phppl(i,1)=HeatingFuhe(i,1)/COPpl(i,1);
                    Qhpscpl(i,1)=HeatingFuhe(i,1)-Phppl(i,1);
                elseif Xhp(i,1)>X(j,4)
                    Xhp02(i,1)=fix(X(j,4));
                    Qdeloss(i,1)=HeatingFuhe(i,1)-Qhpsh(i,1)*Xhp02(i,1);                          %Qdeloss����ȱʧ������
                    Xloss(i,1)=Qdeloss(i,1)/HeatingFuhe(i,1);
                    Phppl(i,1)=Php(i,1)*Xhp02(i,1);
                    Qhpscpl(i,1)=Xhp02(i,1)*Qhpsc(i,1);
                else
                    Qdeloss(i,1)=HeatingFuhe(i,1);                          %Qdeloss����ȱʧ������
                    Xloss(i,1)=1;
                    Phppl(i,1)=0;
                    Qhpscpl(i,1)=0;
                end                                                                      %Xhp��ʾ���յ��ȱ���������С���㣬���2 X(j,4)��ʾ���������߽�ֵ
            
                %�������ܵ�����
                Xge01(i,1)=Qhpscpl(i,1)/(6*0.95);
                Xge(i,1)=fix(Xge(i,1));
            
                if Xge01(i,1)-Xge(i,1) > 0                                                %�Ե��������ȡ��������+1
                    Xge(i,1)=fix(Xge01(i,1))+1;
                else
                    Xge(i,1)=Xge01(i,1);
                end
                %%--�ȱü������--%%
            end
            
            %���1/2/3
            if Qtes(i,1)>=HeatingFuhe(i,1)
                Ttes(i+1,1)=(X(j,5)*Ktes*Ates*(Temperature(i,1)-Ttes(i,1))-HeatingFuhe(i,1)*3600)/(X(j,5)*mtes*Cw)+Ttes(i,1);
                Xhp02(i,1)=0;
            elseif Qtes(i,1)<HeatingFuhe(i,1) && Qtes(i,1)>0
                Qhpbui(i,1)=HeatingFuhe(i,1)-Qtes(i,1);
                Ttes(i+1,1)=Tbui;
                
                %���2/3
                %��һ�μ����ȱ�����
                Xhp01(i,1)=Qhpbui(i,1)/Qhpsh(i,1);
                
                if Xhp01(i,1)<1
                    Xhp02(i,1)=0;
                else
                    Xhp02(i,1)=fix(Xhp01(i,1))+1;
                end
                
                Xhp(i,1)=Xhp02(i,1);
                
                %�ڶ��μ����ȱ����������ݱ߽�����,������ȱʧ���벿�ָ����µ��ȱù��ġ�COP
                if Xhp(i,1) <= X(j,4) && Xhp(i,1)>0
                    PL(i,1)=Qhpbui(i,1)/(Qhpsh(i,1)*Xhp02(i,1));                          %���1�����ڲ��ַ����µ��ȱ�
                    COPnom(i,1)=Qhpsh(i,1)/Php(i,1);
                    COPpl(i,1)=COPnom(i,1)*(-0.00006*PL(i,1)^6+0.0017*PL(i,1)^5-0.0181*PL(i,1)^4+0.096*PL(i,1)^3-0.2697*PL(i,1)^2+0.4276*PL(i,1)+0.7626);
                    Phppl(i,1)=Qhpbui(i,1)/COPpl(i,1);
                    Qhpscpl(i,1)=Qhpbui(i,1)-Phppl(i,1);
                elseif Xhp(i,1)> X(j,4)
                    Xhp02(i,1)=fix(X(j,4));
                    Qdeloss(i,1)=Qhpbui(i,1)-Qhpsh(i,1)*Xhp02(i,1);                          %Qdeloss����ȱʧ������
                    Xloss(i,1)=Qdeloss(i,1)/Qhpbui(i,1);
                    Phppl(i,1)=Php(i,1)*Xhp02(i,1);
                    Qhpscpl(i,1)=Xhp02(i,1)*Qhpsc(i,1);
                else
                    Qdeloss(i,1)=Qhpbui(i,1);                          %Qdeloss����ȱʧ������
                    Xloss(i,1)=1;
                    Phppl(i,1)=0;
                    Qhpscpl(i,1)=0;
                end                                                                      %Xhp��ʾ���յ��ȱ���������С���㣬���2 X(j,4)��ʾ���������߽�ֵ
                
                %������������
                Xge01(i,1)=Qhpscpl(i,1)/(6*0.95);
                Xge(i,1)=fix(Xge(i,1));
                
                if Xge01(i,1)-Xge(i,1) > 0                                               %�Ե��������ȡ��������+1
                    Xge(i,1)=fix(Xge01(i,1))+1;
                else
                    Xge(i,1)=Xge01(i,1);
                end
                
            end
            %%--�ȱü������--%%
            
            %%--�縺�ɼ���--%%
            Pbss(i,1)=Fuhe(i,1)+Phppl(i,1)-Ppv(i,1)*X(j,1)-Pwt(i,1)*X(j,2);
            
            %���س�ŵ�״̬�ж�
            if Pbss(i,1) < 0                                                            %pbssΪ��������Ҫ��磬��ʱ��pbss�������������޵Ĺ���
                SOC(i+1,1)=(1-0.03)*SOC(i,1)-Pbss(i,1)*0.9/(X(j,3)*Eb);
                More(i,1)=abs(Pbss(i,1));                                               %��¼����ĵ���
                elecnum(i,1)=0;
            elseif Pbss(i,1) > 0                                                        %pbssΪ��������Ҫ�ŵ�
                SOC(i+1,1)=(1-0.03)*SOC(i,1)-Pbss(i,1)/(X(j,3)*0.9*Eb);
                More(i,1)=0;                                                            %��¼������˷ѵ���
                elecnum(i,1)=abs(Pbss(i,1));
            else
                SOC(i+1,1)=SOC(i,1);
                More(i,1)=0;                                                            %��¼������˷ѵ���
                elecnum(i,1)=0;
            end
            
            %��������ж�
            if SOC(i+1,1)>=SOCmax
                SOC(i+1,1)=SOCmax;
                %�ȼ��㼴�����ʱ��������ʵ��繦��
                if SOC(i,1)~=SOCmax                                                         %~=�ǲ����ڵ���˼
                    Pbss(i,1)=(SOC(i+1,1)-(1-0.03)*SOC(i,1))*Eb*X(j,3)/(-0.9);
                else
                    Pbss(i,1)=0;
                end
                %��100%�������Ϊ����������
                if More(i,1)>=abs(Pbss(i,1))
                    More(i,1)=More(i,1)-abs(Pbss(i,1));
                else
                    More(i,1)=0;
                end
                %��¼�������
            else
                More(i,1)=0;
            end
            
            %ȱ���ж�
            if SOC(i+1,1)<=SOCmin
                SOC(i+1,1)=SOCmin;
                %���㼴������ʱ��������ʵ�ŵ繦��
                if SOC(i,1)~=SOCmin
                    Pbss(i,1)=(SOC(i+1,1)-(1-0.03)*SOC(i,1))*Eb*(-0.9)*X(j,3);
                else
                    Pbss(i,1)=0;
                end
                %��ʵȱ����
                if elecnum(i,1)>=abs(Pbss(i,1))
                    elecnum(i,1)=elecnum(i,1)-abs(Pbss(i,1));
                else
                    elecnum(i,1)=0;
                end
            else
                elecnum(i,1)=0;
            end
            
            %�⹺��ķ�
            Celebuy(i,1)=elecnum(i,1)*DianJia(i,1);
            
            %�۵�ӯ��
            Celesell(i,1)=More(i,1)*Cse;
            %Emore(i,1)=More(i,1)/(Ppv(i,1)*X(j,1)+Pwt(i,1)*X(j,2)+0.001);
            %%--�縺�ɼ������--%%
            
        end
        %%%---16.00-18.00�������---%%%
        
        
        %%%---18.00-23.00�����ϰ�ʱ�䣬�Ƿ�ȵ��---%%%
        if d(i,1)>0.78 && d(i,1)<0.97
            Ttes(i+1,1)=Ktes*Ates*(Temperature(i,1)-Ttes(i,1))/(mtes*Cw)+Ttes(i,1);
            
            %%%---�縺�ɼ���---%%%
            Pbss(i,1)=Fuhe(i,1)-Ppv(i,1)*X(j,1)-Pwt(i,1)*X(j,2);
            
            %---���س�ŵ�״̬�ж�---%
            if Pbss(i,1) < 0                                                        %pbssΪ��������Ҫ��磬��ʱ��pbss�������������޵Ĺ���
                SOC(i+1,1)=(1-0.03)*SOC(i,1)-Pbss(i,1)*0.9/(X(j,3)*Eb);
                More(i,1)=abs(Pbss(i,1));                                           %��¼����ĵ���
                elecnum(i,1)=0;
            elseif   Pbss(i,1) > 0                                                  %pbssΪ��������Ҫ�ŵ�
                SOC(i+1,1)=(1-0.03)*SOC(i,1)-Pbss(i,1)/(X(j,3)*0.9*Eb);
                More(i,1)=0;                                                        %��¼������˷ѵ���
                elecnum(i,1)=abs(Pbss(i,1));
            else
                SOC(i+1,1)=SOC(i,1);
                More(i,1)=0;                                                        %��¼������˷ѵ���
                elecnum(i,1)=0;
            end
            %---���س�ŵ�״̬�жϽ���---%
            
            %---��������ж�---%
            if  SOC(i+1,1)>=SOCmax
                SOC(i+1,1)=SOCmax;
                %�ȼ��㼴�����ʱ��������ʵ��繦��
                if  SOC(i,1)~=SOCmax                                                %~=�ǲ����ڵ���˼
                    Pbss(i,1)=(SOC(i+1,1)-(1-0.03)*SOC(i,1))*Eb*X(j,3)/(-0.9);
                else
                    Pbss(i,1)=0;
                end
                %��¼�������
                %��100%�������Ϊ����������
                if More(i,1)>=abs(Pbss(i,1))
                    More(i,1)=More(i,1)-abs(Pbss(i,1));                             %������õĺ�����
                else
                    More(i,1)=0;
                end
                %��¼�������
            else
                More(i,1)=0;
            end
            %---��������жϽ���---%
            
            %---ȱ���ж�---%
            if  SOC(i+1,1)<=SOCmin
                SOC(i+1,1)=SOCmin;
                %���㼴������ʱ��������ʵ�ŵ繦��
                if  SOC(i,1)~=SOCmin
                    Pbss(i,1)=(SOC(i+1,1)-(1-0.03)*SOC(i,1))*Eb*(-0.9)*X(j,3);
                else
                    Pbss(i,1)=0;
                end
                %�ڼ���ʱ����ʵȱ����
                if elecnum(i,1)>=abs(Pbss(i,1))
                    elecnum(i,1)=elecnum(i,1)-abs(Pbss(i,1));                       %������õĺ�����
                else
                    elecnum(i,1)=0;
                end
            else
                elecnum(i,1)=0;
            end
            
            %�⹺��ķ�
            Celebuy(i,1)=0;
            
            %�۵�ӯ��
            Celesell(i,1)=More(i,1)*Cse;
            %Emore(i,1)=More(i,1)/(Ppv(i,1)*X(j,1)+Pwt(i,1)*X(j,2)+0.001);
            
            %---ȱ���жϽ���---%
            %%%---��ģ�����---%%%
        end
        %%%---18.00-23.00�������---%%%
        
        %%%---�˼��ڼ������---%%%
    end
    
    %% ���ɼ���2521-3624h
    for i=2521:3624
        
        %%%---�縺�ɼ���---%%%
        Pbss(i,1)=Fuhe(i,1)-Ppv(i,1)*X(j,1)-Pwt(i,1)*X(j,2);
        
        %---���س�ŵ�״̬�ж�---%
        if Pbss(i,1) < 0                                                        %pbssΪ��������Ҫ��磬��ʱ��pbss�������������޵Ĺ���
            SOC(i+1,1)=(1-0.03)*SOC(i,1)-Pbss(i,1)*0.9/(X(j,3)*Eb);
            More(i,1)=abs(Pbss(i,1));                                           %��¼����ĵ���
            elecnum(i,1)=0;
        elseif   Pbss(i,1) > 0                                                  %pbssΪ��������Ҫ�ŵ�
            SOC(i+1,1)=(1-0.03)*SOC(i,1)-Pbss(i,1)/(X(j,3)*0.9*Eb);
            More(i,1)=0;                                                        %��¼������˷ѵ���
            elecnum(i,1)=abs(Pbss(i,1));
        else
            SOC(i+1,1)=SOC(i,1);
            More(i,1)=0;                                                        %��¼������˷ѵ���
            elecnum(i,1)=0;
        end
        %---���س�ŵ�״̬�жϽ���---%
        
        %---��������ж�---%
        if  SOC(i+1,1)>=SOCmax
            SOC(i+1,1)=SOCmax;
            %�ȼ��㼴�����ʱ��������ʵ��繦��
            if  SOC(i,1)~=SOCmax                                                %~=�ǲ����ڵ���˼
                Pbss(i,1)=(SOC(i+1,1)-(1-0.03)*SOC(i,1))*Eb*X(j,3)/(-0.9);
            else
                Pbss(i,1)=0;
            end
            %��¼�������
            %��100%�������Ϊ����������
            if More(i,1)>=abs(Pbss(i,1))
                More(i,1)=More(i,1)-abs(Pbss(i,1));                             %������õĺ�����
            else
                More(i,1)=0;
            end
            %��¼�������
        else
            More(i,1)=0;
        end
        %---��������жϽ���---%
        
        %---ȱ���ж�---%
        if  SOC(i+1,1)<=SOCmin
            SOC(i+1,1)=SOCmin;
            %���㼴������ʱ��������ʵ�ŵ繦��
            if  SOC(i,1)~=SOCmin
                Pbss(i,1)=(SOC(i+1,1)-(1-0.03)*SOC(i,1))*Eb*(-0.9)*X(j,3);
            else
                Pbss(i,1)=0;
            end
            %�ڼ���ʱ����ʵȱ����
            if elecnum(i,1)>=abs(Pbss(i,1))
                elecnum(i,1)=elecnum(i,1)-abs(Pbss(i,1));                       %������õĺ�����
            else
                elecnum(i,1)=0;
            end
        else
            elecnum(i,1)=0;
        end
        
        %�⹺��ķ�
        Celebuy(i,1)=elecnum(i,1)*DianJia(i,1);
        
        %�۵�ӯ��
        Celesell(i,1)=More(i,1)*Cse;
        %Emore(i,1)=More(i,1)/(Ppv(i,1)*X(j,1)+Pwt(i,1)*X(j,2)+0.001);
        
        %---ȱ���жϽ���---%
        %%%---��ģ�����---%%%
        
    end
    
    %% ���伾��3625h-5832h
    for i=3625:5832
        
        a(i,1)=i;
        b(i,1)=a(i,1)/24;
        c(i,1)=fix(b(i,1));
        d(i,1)=b(i,1)-c(i,1);
        
        %%%---23.00-07.00�����ʱ�䣬�ȱ���ˮ������---%%%
        if d(i,1)<0.31
            
            if Ttes(i,1)>Ttessdcc                                                                                       %ˮ���¶ȿ�����ͣ�ж�
                Qhptes(i,1)=abs((Ttessdc-Ttes(i,1))*Cw*X(j,5)*mtes+X(j,5)*Ktes*Ates*(Ttes(i,1)-Temperature(i,1)))/3600;      %�ȱ÷�����为��
                
                %%---�ȱò��ָ��ɼ�ȱʧ����---%%
                
                %��һ�μ����ȱ�����
                Xhp01(i,1)=Qhptes(i,1)/Qhpsc(i,1);
                
                if Xhp01(i,1)<1
                    Xhp02(i,1)=0;
                else
                    Xhp02(i,1)=fix(Xhp01(i,1))+1;
                end
                
                Xhp(i,1)=Xhp02(i,1);
                
                %�ڶ��μ����ȱ����������ݱ߽�����,������ȱʧ���벿�ָ����µ��ȱù��ġ�COP
                if Xhp(i,1)<=X(j,4) && Xhp(i,1)>0
                    PL(i,1)=Qhptes(i,1)/(Qhpsc(i,1)*Xhp02(i,1));                          %���1�����ڲ��ַ����µ��ȱ�
                    COPnom(i,1)=Qhpsc(i,1)/Php(i,1);
                    COPpl(i,1)=COPnom(i,1)*(-0.00006*PL(i,1)^6+0.0017*PL(i,1)^5-0.0181*PL(i,1)^4+0.096*PL(i,1)^3-0.2697*PL(i,1)^2+0.4276*PL(i,1)+0.7626);
                    Phppl(i,1)=Qhptes(i,1)/COPpl(i,1);                                    %���㲿�ָ����¸��ȱõĹ���
                    Qhpshpl(i,1)=Qhptes(i,1)+Phppl(i,1);
                    Ttes(i+1,1)=Ttessdc;
                elseif Xhp(i,1)>X(j,4)
                    Xhp02(i,1)=fix(X(j,4));
                    %Qdeloss(i,1)=Qhptes(i,1)-Qhpsh(i,1)*Xhp02(i,1);                          %Qdeloss����ȱʧ������
                    %Xloss(i,1)=Qdeloss(i,1)/Qhptes(i,1);
                    Phppl(i,1)=Php(i,1)*Xhp02(i,1);
                    Qhpshpl(i,1)=Xhp02(i,1)*Qhpsh(i,1);
                    Ttes(i+1,1)=(-Xhp02(i,1)*Qhpsc(i,1)*3600-X(j,5)*Ktes*Ates*(Ttes(i,1)-Temperature(i,1)))/(mtes*Cw*X(j,5))+Ttes(i,1);
                else
                    Phppl(i,1)=0;
                    Qhpshpl(i,1)=0;
                    Ttes(i+1,1)=Ktes*Ates*(Temperature(i,1)-Ttes(i,1))/(mtes*Cw)+Ttes(i,1);
                end                                                                      %Xhp��ʾ���յ��ȱ���������С���㣬���2 X(j,4)��ʾ���������߽�ֵ
                
                Phptes(i,1)=Phppl(i,1);
                
                %������������
                Xge01(i,1)=Qhpshpl(i,1)/(6*0.95);
                Xge(i,1)=fix(Xge01(i,1));
                if Xge01(i,1)-Xge(i,1) > 0                                                %�Ե��������ȡ��������+1
                    Xge(i,1)=fix(Xge01(i,1))+1;
                else
                    Xge(i,1)=Xge01(i,1);
                end
                
            else
                Ttes(i+1,1)=Ktes*Ates*(Temperature(i,1)-Ttes(i,1))/(mtes*Cw)+Ttes(i,1);
            end
            
            %%---�ȱò��ָ��ɼ�ȱʧ�������---%%
            
            %%%---�縺�ɼ���---%%%
            Pbss(i,1)=Fuhe(i,1)-Ppv(i,1)*X(j,1)-Pwt(i,1)*X(j,2);
            
            %---���س�ŵ�״̬�ж�---%
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
            %---���س�ŵ�״̬�жϽ���---%
            
            %---��������ж�---%
            if  SOC(i+1,1)>=SOCmax
                SOC(i+1,1)=SOCmax;
                %�ȼ��㼴�����ʱ��������ʵ��繦��
                if  SOC(i,1)~=SOCmax
                    Pbss(i,1)=(SOC(i+1,1)-(1-0.03)*SOC(i,1))*Eb*X(j,3)/(-0.9);
                else
                    Pbss(i,1)=0;
                end
                %��¼�������
                %��100%�������Ϊ����������
                if More(i,1)>=abs(Pbss(i,1))
                    More(i,1)=More(i,1)-abs(Pbss(i,1));
                else
                    More(i,1)=0;
                end
                %��¼�������
            else
                More(i,1)=0;
            end
            %---��������жϽ���---%
            
            %---ȱ���ж�---%
            if  SOC(i+1,1)<=SOCmin
                SOC(i+1,1)=SOCmin;
                %���㼴������ʱ��������ʵ�ŵ繦��
                if  SOC(i,1)~=SOCmin
                    Pbss(i,1)=(SOC(i+1,1)-(1-0.03)*SOC(i,1))*Eb*(-0.9)*X(j,3);
                else
                    Pbss(i,1)=0;
                end
                %�ڼ���ʱ����ʵȱ����
                if elecnum(i,1)>=abs(Pbss(i,1))
                    elecnum(i,1)=elecnum(i,1)-abs(Pbss(i,1));
                else
                    elecnum(i,1)=0;
                end
            else
                elecnum(i,1)=0;
            end
            
            %�⹺��ķ�
            Celebuy(i,1)=Phptes(i,1)*DianJia(i,1);
            
            %�۵�ӯ��
            Celesell(i,1)=More(i,1)*Cse;
            %Emore(i,1)=More(i,1)/(Ppv(i,1)*X(j,1)+Pwt(i,1)*X(j,2)+0.001);
            
        end
        
        %%%---23.00-07.00�������---%%%
        
        %%%---07.00-08.00���ǹ���ʱ�䣬�Ƿ�ȵ��---%%%
        if d(i,1)>0.31 && d(i,1)<0.36
            Ttes(i+1,1)=Ktes*Ates*(Temperature(i,1)-Ttes(i,1))/(mtes*Cw)+Ttes(i,1);
            
            %%%---�縺�ɼ���---%%%
            Pbss(i,1)=Fuhe(i,1)-Ppv(i,1)*X(j,1)-Pwt(i,1)*X(j,2);
            
            %---���س�ŵ�״̬�ж�---%
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
            %---���س�ŵ�״̬�жϽ���---%
            
            %---��������ж�---%
            if  SOC(i+1,1)>=SOCmax
                SOC(i+1,1)=SOCmax;
                %�ȼ��㼴�����ʱ��������ʵ��繦��
                if  SOC(i,1)~=SOCmax
                    Pbss(i,1)=(SOC(i+1,1)-(1-0.03)*SOC(i,1))*Eb*X(j,3)/(-0.9);
                else
                    Pbss(i,1)=0;
                end
                %��¼�������
                %��100%�������Ϊ����������
                if More(i,1)>=abs(Pbss(i,1))
                    More(i,1)=More(i,1)-abs(Pbss(i,1));
                else
                    More(i,1)=0;
                end
                %��¼�������
            else
                More(i,1)=0;
            end
            %---��������жϽ���---%
            
            %---ȱ���ж�---%
            if  SOC(i+1,1)<=SOCmin
                SOC(i+1,1)=SOCmin;
                %���㼴������ʱ��������ʵ�ŵ繦��
                if  SOC(i,1)~=SOCmin
                    Pbss(i,1)=(SOC(i+1,1)-(1-0.03)*SOC(i,1))*Eb*(-0.9)*X(j,3);
                else
                    Pbss(i,1)=0;
                end
                %�ڼ���ʱ����ʵȱ����
                if elecnum(i,1)>=abs(Pbss(i,1))
                    elecnum(i,1)=elecnum(i,1)-abs(Pbss(i,1));
                else
                    elecnum(i,1)=0;
                end
            else
                elecnum(i,1)=0;
            end
            
            %�⹺��ķ�
            Celebuy(i,1)=0;
            
            %�۵�ӯ��
            Celesell(i,1)=More(i,1)*Cse;
            %Emore(i,1)=More(i,1)/(Ppv(i,1)*X(j,1)+Pwt(i,1)*X(j,2)+0.001);
            
        end
        %%%---07.00-08.00�������---%%%
        
        %%%---8.00-9.00,����ʱ�䣬ƽ�����---%%%
        if d(i,1)>0.36 && d(i,1)<0.39
            
            %%--ˮ���¶ȼ���--%%
            Ttes(i+1,1)=Ktes*Ates*(Temperature(i,1)-Ttes(i,1))/(mtes*Cw)+Ttes(i,1);
            %%--ˮ���¶ȼ������--%%
            
            %%--�ȱü���--%%
            
            %��һ�μ����ȱ�����
            Xhp01(i,1)=CoolingFuhe(i,1)/Qhpsc(i,1);
            
            if Xhp01(i,1)<1
                Xhp02(i,1)=0;
            else
                Xhp02(i,1)=fix(Xhp01(i,1))+1;
            end
            
            Xhp(i,1)=Xhp02(i,1);
                
            %�ڶ��μ����ȱ����������ݱ߽�����,������ȱʧ���벿�ָ����µ��ȱù��ġ�COP
            if Xhp(i,1)<=X(j,4) && Xhp(i,1)>0
                PL(i,1)=CoolingFuhe(i,1)/(Qhpsc(i,1)*Xhp02(i,1));                          %���1�����ڲ��ַ����µ��ȱ�
                COPnom(i,1)=Qhpsc(i,1)/Php(i,1);
                COPpl(i,1)=COPnom(i,1)*(-0.00006*PL(i,1)^6+0.0017*PL(i,1)^5-0.0181*PL(i,1)^4+0.096*PL(i,1)^3-0.2697*PL(i,1)^2+0.4276*PL(i,1)+0.7626);
                Phppl(i,1)=CoolingFuhe(i,1)/COPpl(i,1);
                Qhpshpl(i,1)=CoolingFuhe(i,1)+Phppl(i,1);
            elseif Xhp(i,1)>X(j,4)
                Xhp02(i,1)=fix(X(j,4));
                Qdeloss(i,1)=CoolingFuhe(i,1)-Qhpsc(i,1)*Xhp02(i,1);                          %Qdeloss����ȱʧ������
                Xloss(i,1)=Qdeloss(i,1)/CoolingFuhe(i,1);
                Phppl(i,1)=Php(i,1)*Xhp02(i,1);
                Qhpshpl(i,1)=Xhp02(i,1)*Qhpsh(i,1);
            else
                Qdeloss(i,1)=CoolingFuhe(i,1);                          %Qdeloss����ȱʧ������
                Xloss(i,1)=1;
                Phppl(i,1)=0;
                Qhpshpl(i,1)=0;
            end                                                                      %Xhp��ʾ���յ��ȱ���������С���㣬���2 X(j,4)��ʾ���������߽�ֵ
            
            %������������
            Xge01(i,1)=Qhpshpl(i,1)/(6*0.95);
            Xge(i,1)=fix(Xge01(i,1));
            
            if Xge01(i,1)-Xge(i,1) > 0                                                %�Ե��������ȡ��������+1
                Xge(i,1)=fix(Xge01(i,1))+1;
            else
                Xge(i,1)=Xge01(i,1);
            end
            %%---�ȱü������---%%
            
            %%--�縺�ɼ���--%%
            Pbss(i,1)=Fuhe(i,1)+Phppl(i,1)-Ppv(i,1)*X(j,1)-Pwt(i,1)*X(j,2);
            
            %���س�ŵ�״̬�ж�
            if Pbss(i,1) < 0                                                         %pbssΪ��������Ҫ��磬��ʱ��pbss�������������޵Ĺ���
                SOC(i+1,1)=(1-0.03)*SOC(i,1)-Pbss(i,1)*0.9/(X(j,3)*Eb);
                More(i,1)=abs(Pbss(i,1));                                            %��¼����ĵ���
                elecnum(i,1)=0;
            elseif Pbss(i,1) > 0                                                     %pbssΪ��������Ҫ�ŵ�
                SOC(i+1,1)=(1-0.03)*SOC(i,1)-Pbss(i,1)/(X(j,3)*0.9*Eb);
                More(i,1)=0;                                                         %��¼������˷ѵ���
                elecnum(i,1)=abs(Pbss(i,1));
            else
                SOC(i+1,1)=SOC(i,1);
                More(i,1)=0;                                                         %��¼������˷ѵ���
                elecnum(i,1)=0;
            end
            
            %��������ж�
            if SOC(i+1,1)>=SOCmax
                SOC(i+1,1)=SOCmax;
                %�ȼ��㼴�����ʱ��������ʵ��繦��
                if SOC(i,1)~=SOCmax                                                      %~=�ǲ����ڵ���˼
                    Pbss(i,1)=(SOC(i+1,1)-(1-0.03)*SOC(i,1))*Eb*X(j,3)/(-0.9);
                else
                    Pbss(i,1)=0;
                end
                %��100%�������Ϊ����������
                if More(i,1)>=abs(Pbss(i,1))
                    More(i,1)=More(i,1)-abs(Pbss(i,1));
                else
                    More(i,1)=0;
                end
                %��¼�������
            else
                More(i,1)=0;
            end
            
            %ȱ���ж�
            if SOC(i+1,1)<=SOCmin
                SOC(i+1,1)=SOCmin;
                %���㼴������ʱ��������ʵ�ŵ繦��
                if SOC(i,1)~=SOCmin
                    Pbss(i,1)=(SOC(i+1,1)-(1-0.03)*SOC(i,1))*Eb*(-0.9)*X(j,3);
                else
                    Pbss(i,1)=0;
                end
                %��ʵȱ����
                if elecnum(i,1)>=abs(Pbss(i,1))
                    elecnum(i,1)=elecnum(i,1)-abs(Pbss(i,1));
                else
                    elecnum(i,1)=0;
                end
            else
                elecnum(i,1)=0;
            end
            
            %�⹺��ķ�
            Celebuy(i,1)=elecnum(i,1)*DianJia(i,1);
            
            %�۵�ӯ��
            Celesell(i,1)=More(i,1)*Cse;
            %Emore(i,1)=More(i,1)/(Ppv(i,1)*X(j,1)+Pwt(i,1)*X(j,2)+0.001);
            %%--�縺�ɼ������--%%
            
        end
        %%%---8.00-9.00�������---%%%
        
        %%%---9.00-12.00,����ʱ�䣬��ֵ��ۣ�����ˮ��---%%%
        if d(i,1)>0.39 && d(i,1)<0.51
            
            %%--�ȱü���--%%
            Qtes(i,1)=((Tbuc-Ttes(i,1))*Cw*X(j,5)*mtes+X(j,5)*Ktes*Ates*(Ttes(i,1)-Temperature(i,1)))/3600;
            
            if Qtes(i,1)<=0
                Qtes(i,1)=0;
                
                %%--ˮ���¶ȼ���--%%
                Ttes(i+1,1)=Ktes*Ates*(Temperature(i,1)-Ttes(i,1))/(mtes*Cw)+Ttes(i,1);
                %%--ˮ���¶ȼ������--%%
            
                %%--�ȱü���--%%
            
                %��һ�μ����ȱ�����
                Xhp01(i,1)=CoolingFuhe(i,1)/Qhpsc(i,1);
                
                if Xhp01(i,1)<1
                    Xhp02(i,1)=0;
                else
                    Xhp02(i,1)=fix(Xhp01(i,1))+1;
                end
                
                Xhp(i,1)=Xhp02(i,1);
                
                %�ڶ��μ����ȱ����������ݱ߽�����,������ȱʧ���벿�ָ����µ��ȱù��ġ�COP
                if Xhp(i,1)<=X(j,4) && Xhp(i,1)>0
                    PL(i,1)=CoolingFuhe(i,1)/(Qhpsc(i,1)*Xhp02(i,1));                          %���1�����ڲ��ַ����µ��ȱ�
                    COPnom(i,1)=Qhpsc(i,1)/Php(i,1);
                    COPpl(i,1)=COPnom(i,1)*(-0.00006*PL(i,1)^6+0.0017*PL(i,1)^5-0.0181*PL(i,1)^4+0.096*PL(i,1)^3-0.2697*PL(i,1)^2+0.4276*PL(i,1)+0.7626);
                    Phppl(i,1)=CoolingFuhe(i,1)/COPpl(i,1);
                    Qhpshpl(i,1)=CoolingFuhe(i,1)+Phppl(i,1);
                elseif Xhp(i,1)>X(j,4)
                    Xhp02(i,1)=fix(X(j,4));
                    Qdeloss(i,1)=CoolingFuhe(i,1)-Qhpsc(i,1)*Xhp02(i,1);                          %Qdeloss����ȱʧ������
                    Xloss(i,1)=Qdeloss(i,1)/CoolingFuhe(i,1);
                    Phppl(i,1)=Php(i,1)*Xhp02(i,1);
                    Qhpshpl(i,1)=Xhp02(i,1)*Qhpsh(i,1);
                else
                    Qdeloss(i,1)=CoolingFuhe(i,1);                          %Qdeloss����ȱʧ������
                    Xloss(i,1)=1;
                    Phppl(i,1)=0;
                    Qhpshpl(i,1)=0;
                end                                                                      %Xhp��ʾ���յ��ȱ���������С���㣬���2 X(j,4)��ʾ���������߽�ֵ
            
                %������������
                Xge01(i,1)=Qhpshpl(i,1)/(6*0.95);
                Xge(i,1)=fix(Xge01(i,1));
            
                if Xge01(i,1)-Xge(i,1) > 0                                                %�Ե��������ȡ��������+1
                    Xge(i,1)=fix(Xge01(i,1))+1;
                else
                    Xge(i,1)=Xge01(i,1);
                end
                %%---�ȱü������---%%

            end
            
            %���1/2/3
            if Qtes(i,1)>=CoolingFuhe(i,1)
                Ttes(i+1,1)=(X(j,5)*Ktes*Ates*(Temperature(i,1)-Ttes(i,1))+CoolingFuhe(i,1)*3600)/(X(j,5)*mtes*Cw)+Ttes(i,1);
                Xhp02(i,1)=0;
            elseif Qtes(i,1)<CoolingFuhe(i,1) && Qtes(i,1)>0
                Qhpbui(i,1)=CoolingFuhe(i,1)-Qtes(i,1);
                Ttes(i+1,1)=Tbuc;
                
                %���2/3
                %��һ�μ����ȱ�����
                Xhp01(i,1)=Qhpbui(i,1)/Qhpsc(i,1);
                
                if Xhp01(i,1)<1
                    Xhp02(i,1)=0;
                else
                    Xhp02(i,1)=fix(Xhp01(i,1))+1;
                end
                
                Xhp(i,1)=Xhp02(i,1);
                
                %�ڶ��μ����ȱ����������ݱ߽�����,������ȱʧ���벿�ָ����µ��ȱù��ġ�COP
                if Xhp(i,1)<=X(j,4) && Xhp(i,1)>0
                    PL(i,1)=Qhpbui(i,1)/(Qhpsc(i,1)*Xhp02(i,1));                          %���1�����ڲ��ַ����µ��ȱ�
                    COPnom(i,1)=Qhpsc(i,1)/Php(i,1);
                    COPpl(i,1)=COPnom(i,1)*(-0.00006*PL(i,1)^6+0.0017*PL(i,1)^5-0.0181*PL(i,1)^4+0.096*PL(i,1)^3-0.2697*PL(i,1)^2+0.4276*PL(i,1)+0.7626);
                    Phppl(i,1)=Qhpbui(i,1)/COPpl(i,1);
                    Qhpshpl(i,1)=Qhpbui(i,1)+Phppl(i,1);
                elseif Xhp(i,1)>X(j,4)
                    Xhp02(i,1)=fix(X(j,4));
                    Qdeloss(i,1)=Qhpbui(i,1)-Qhpsc(i,1)*Xhp02(i,1);                          %Qdeloss����ȱʧ������
                    Xloss(i,1)=Qdeloss(i,1)/Qhpbui(i,1);
                    Phppl(i,1)=Php(i,1)*Xhp02(i,1);
                    Qhpshpl(i,1)=Xhp02(i,1)*Qhpsh(i,1);
                else
                    Qdeloss(i,1)=Qhpbui(i,1);                          %Qdeloss����ȱʧ������
                    Xloss(i,1)=1;
                    Phppl(i,1)=0;
                    Qhpshpl(i,1)=0;
                end                                                                      %Xhp��ʾ���յ��ȱ���������С���㣬���2 X(j,4)��ʾ���������߽�ֵ
                
                %������������
                Xge01(i,1)=Qhpshpl(i,1)/(6*0.95);
                Xge(i,1)=fix(Xge(i,1));
                
                if Xge01(i,1)-Xge(i,1) > 0                                               %�Ե��������ȡ��������+1
                    Xge(i,1)=fix(Xge01(i,1))+1;
                else
                    Xge(i,1)=Xge01(i,1);
                end
                
            end
            %%--�ȱü������--%%
            
            %%--�縺�ɼ���--%%
            Pbss(i,1)=Fuhe(i,1)+Phppl(i,1)-Ppv(i,1)*X(j,1)-Pwt(i,1)*X(j,2);
            
            %���س�ŵ�״̬�ж�
            if Pbss(i,1) < 0                                                            %pbssΪ��������Ҫ��磬��ʱ��pbss�������������޵Ĺ���
                SOC(i+1,1)=(1-0.03)*SOC(i,1)-Pbss(i,1)*0.9/(X(j,3)*Eb);
                More(i,1)=abs(Pbss(i,1));                                               %��¼����ĵ���
                elecnum(i,1)=0;
            elseif Pbss(i,1) > 0                                                        %pbssΪ��������Ҫ�ŵ�
                SOC(i+1,1)=(1-0.03)*SOC(i,1)-Pbss(i,1)/(X(j,3)*0.9*Eb);
                More(i,1)=0;                                                            %��¼������˷ѵ���
                elecnum(i,1)=abs(Pbss(i,1));
            else
                SOC(i+1,1)=SOC(i,1);
                More(i,1)=0;                                                            %��¼������˷ѵ���
                elecnum(i,1)=0;
            end
            
            %��������ж�
            if SOC(i+1,1)>=SOCmax
                SOC(i+1,1)=SOCmax;
                %�ȼ��㼴�����ʱ��������ʵ��繦��
                if SOC(i,1)~=SOCmax                                                         %~=�ǲ����ڵ���˼
                    Pbss(i,1)=(SOC(i+1,1)-(1-0.03)*SOC(i,1))*Eb*X(j,3)/(-0.9);
                else
                    Pbss(i,1)=0;
                end
                %��100%�������Ϊ����������
                if More(i,1)>=abs(Pbss(i,1))
                    More(i,1)=More(i,1)-abs(Pbss(i,1));
                else
                    More(i,1)=0;
                end
                %��¼�������
            else
                More(i,1)=0;
            end
            
            %ȱ���ж�
            if SOC(i+1,1)<=SOCmin
                SOC(i+1,1)=SOCmin;
                %���㼴������ʱ��������ʵ�ŵ繦��
                if SOC(i,1)~=SOCmin
                    Pbss(i,1)=(SOC(i+1,1)-(1-0.03)*SOC(i,1))*Eb*(-0.9)*X(j,3);
                else
                    Pbss(i,1)=0;
                end
                %��ʵȱ����
                if elecnum(i,1)>=abs(Pbss(i,1))
                    elecnum(i,1)=elecnum(i,1)-abs(Pbss(i,1));
                else
                    elecnum(i,1)=0;
                end
            else
                elecnum(i,1)=0;
            end
            
            %�⹺��ķ�
            Celebuy(i,1)=elecnum(i,1)*DianJia(i,1);
            
            %�۵�ӯ��
            Celesell(i,1)=More(i,1)*Cse;
            %Emore(i,1)=More(i,1)/(Ppv(i,1)*X(j,1)+Pwt(i,1)*X(j,2)+0.001);
            %%--�縺�ɼ������--%%
            
        end
        %%%---9.00-12.00�������---%%%
        
        %%%---12.00-16.00,����ʱ�䣬ƽ�����---%%%
        if d(i,1)>0.51 && d(i,1)<0.69
            
            %%--ˮ���¶ȼ���--%%
            Ttes(i+1,1)=Ktes*Ates*(Temperature(i,1)-Ttes(i,1))/(mtes*Cw)+Ttes(i,1);
            %%--ˮ���¶ȼ������--%%
            
            %%--�ȱü���--%%
            
            %��һ�μ����ȱ�����
            Xhp01(i,1)=CoolingFuhe(i,1)/Qhpsc(i,1);
            
            if Xhp01(i,1)<1
                Xhp02(i,1)=0;
            else
                Xhp02(i,1)=fix(Xhp01(i,1))+1;
            end
            
            Xhp(i,1)=Xhp02(i,1);
                
            %�ڶ��μ����ȱ����������ݱ߽�����,������ȱʧ���벿�ָ����µ��ȱù��ġ�COP
            if Xhp(i,1)<=X(j,4) && Xhp(i,1)>0
                PL(i,1)=CoolingFuhe(i,1)/(Qhpsc(i,1)*Xhp02(i,1));                          %���1�����ڲ��ַ����µ��ȱ�
                COPnom(i,1)=Qhpsc(i,1)/Php(i,1);
                COPpl(i,1)=COPnom(i,1)*(-0.00006*PL(i,1)^6+0.0017*PL(i,1)^5-0.0181*PL(i,1)^4+0.096*PL(i,1)^3-0.2697*PL(i,1)^2+0.4276*PL(i,1)+0.7626);
                Phppl(i,1)=CoolingFuhe(i,1)/COPpl(i,1);
                Qhpshpl(i,1)=CoolingFuhe(i,1)+Phppl(i,1);
            elseif Xhp(i,1)>X(j,4)
                Xhp02(i,1)=fix(X(j,4));
                Qdeloss(i,1)=CoolingFuhe(i,1)-Qhpsc(i,1)*Xhp02(i,1);                          %Qdeloss����ȱʧ������
                Xloss(i,1)=Qdeloss(i,1)/CoolingFuhe(i,1);
                Phppl(i,1)=Php(i,1)*Xhp02(i,1);
                Qhpshpl(i,1)=Xhp02(i,1)*Qhpsh(i,1);
            else
                Qdeloss(i,1)=CoolingFuhe(i,1);                          %Qdeloss����ȱʧ������
                Xloss(i,1)=1;
                Phppl(i,1)=0;
                Qhpshpl(i,1)=0;
            end                                                                      %Xhp��ʾ���յ��ȱ���������С���㣬���2 X(j,4)��ʾ���������߽�ֵ
            
            %������������
            Xge01(i,1)=Qhpshpl(i,1)/(6*0.95);
            Xge(i,1)=fix(Xge01(i,1));
            
            if Xge01(i,1)-Xge(i,1) > 0                                                %�Ե��������ȡ��������+1
                Xge(i,1)=fix(Xge01(i,1))+1;
            else
                Xge(i,1)=Xge01(i,1);
            end
            %%---�ȱü������---%%
            
            %%--�縺�ɼ���--%%
            Pbss(i,1)=Fuhe(i,1)+Phppl(i,1)-Ppv(i,1)*X(j,1)-Pwt(i,1)*X(j,2);
            
            %���س�ŵ�״̬�ж�
            if Pbss(i,1) < 0                                                         %pbssΪ��������Ҫ��磬��ʱ��pbss�������������޵Ĺ���
                SOC(i+1,1)=(1-0.03)*SOC(i,1)-Pbss(i,1)*0.9/(X(j,3)*Eb);
                More(i,1)=abs(Pbss(i,1));                                            %��¼����ĵ���
                elecnum(i,1)=0;
            elseif Pbss(i,1) > 0                                                     %pbssΪ��������Ҫ�ŵ�
                SOC(i+1,1)=(1-0.03)*SOC(i,1)-Pbss(i,1)/(X(j,3)*0.9*Eb);
                More(i,1)=0;                                                         %��¼������˷ѵ���
                elecnum(i,1)=abs(Pbss(i,1));
            else
                SOC(i+1,1)=SOC(i,1);
                More(i,1)=0;                                                         %��¼������˷ѵ���
                elecnum(i,1)=0;
            end
            
            %��������ж�
            if SOC(i+1,1)>=SOCmax
                SOC(i+1,1)=SOCmax;
                %�ȼ��㼴�����ʱ��������ʵ��繦��
                if SOC(i,1)~=SOCmax                                                      %~=�ǲ����ڵ���˼
                    Pbss(i,1)=(SOC(i+1,1)-(1-0.03)*SOC(i,1))*Eb*X(j,3)/(-0.9);
                else
                    Pbss(i,1)=0;
                end
                %��100%�������Ϊ����������
                if More(i,1)>=abs(Pbss(i,1))
                    More(i,1)=More(i,1)-abs(Pbss(i,1));
                else
                    More(i,1)=0;
                end
                %��¼�������
            else
                More(i,1)=0;
            end
            
            %ȱ���ж�
            if SOC(i+1,1)<=SOCmin
                SOC(i+1,1)=SOCmin;
                %���㼴������ʱ��������ʵ�ŵ繦��
                if SOC(i,1)~=SOCmin
                    Pbss(i,1)=(SOC(i+1,1)-(1-0.03)*SOC(i,1))*Eb*(-0.9)*X(j,3);
                else
                    Pbss(i,1)=0;
                end
                %��ʵȱ����
                if elecnum(i,1)>=abs(Pbss(i,1))
                    elecnum(i,1)=elecnum(i,1)-abs(Pbss(i,1));
                else
                    elecnum(i,1)=0;
                end
            else
                elecnum(i,1)=0;
            end
            
            %�⹺��ķ�
            Celebuy(i,1)=elecnum(i,1)*DianJia(i,1);
            
            %�۵�ӯ��
            Celesell(i,1)=More(i,1)*Cse;
            %Emore(i,1)=More(i,1)/(Ppv(i,1)*X(j,1)+Pwt(i,1)*X(j,2)+0.001);
            %%--�縺�ɼ������--%%
            
        end
        %%%---12.00-16.00�������---%%%
        
        %%%---16.00-18.00,����ʱ�䣬��ֵ��ۣ�����ˮ��---%%%
        if d(i,1)>0.69 && d(i,1)<0.78
            
            %%--�ȱü���--%%
            Qtes(i,1)=((Tbuc-Ttes(i,1))*Cw*X(j,5)*mtes+X(j,5)*Ktes*Ates*(Ttes(i,1)-Temperature(i,1)))/3600;
            
            if Qtes(i,1)<=0
                Qtes(i,1)=0;
                
                %%--ˮ���¶ȼ���--%%
                Ttes(i+1,1)=Ktes*Ates*(Temperature(i,1)-Ttes(i,1))/(mtes*Cw)+Ttes(i,1);
                %%--ˮ���¶ȼ������--%%
            
                %%--�ȱü���--%%
            
                %��һ�μ����ȱ�����
                Xhp01(i,1)=CoolingFuhe(i,1)/Qhpsc(i,1);
                
                if Xhp01(i,1)<1
                    Xhp02(i,1)=0;
                else
                    Xhp02(i,1)=fix(Xhp01(i,1))+1;
                end
                
                Xhp(i,1)=Xhp02(i,1);
                
                %�ڶ��μ����ȱ����������ݱ߽�����,������ȱʧ���벿�ָ����µ��ȱù��ġ�COP
                if Xhp(i,1)<=X(j,4) && Xhp(i,1)>0
                    PL(i,1)=CoolingFuhe(i,1)/(Qhpsc(i,1)*Xhp02(i,1));                          %���1�����ڲ��ַ����µ��ȱ�
                    COPnom(i,1)=Qhpsc(i,1)/Php(i,1);
                    COPpl(i,1)=COPnom(i,1)*(-0.00006*PL(i,1)^6+0.0017*PL(i,1)^5-0.0181*PL(i,1)^4+0.096*PL(i,1)^3-0.2697*PL(i,1)^2+0.4276*PL(i,1)+0.7626);
                    Phppl(i,1)=CoolingFuhe(i,1)/COPpl(i,1);
                    Qhpshpl(i,1)=CoolingFuhe(i,1)+Phppl(i,1);
                elseif Xhp(i,1)>X(j,4)
                    Xhp02(i,1)=fix(X(j,4));
                    Qdeloss(i,1)=CoolingFuhe(i,1)-Qhpsc(i,1)*Xhp02(i,1);                          %Qdeloss����ȱʧ������
                    Xloss(i,1)=Qdeloss(i,1)/CoolingFuhe(i,1);
                    Phppl(i,1)=Php(i,1)*Xhp02(i,1);
                    Qhpshpl(i,1)=Xhp02(i,1)*Qhpsh(i,1);
                else
                    Qdeloss(i,1)=CoolingFuhe(i,1);                          %Qdeloss����ȱʧ������
                    Xloss(i,1)=1;
                    Phppl(i,1)=0;
                    Qhpshpl(i,1)=0;
                end                                                                      %Xhp��ʾ���յ��ȱ���������С���㣬���2 X(j,4)��ʾ���������߽�ֵ
            
                %������������
                Xge01(i,1)=Qhpshpl(i,1)/(6*0.95);
                Xge(i,1)=fix(Xge01(i,1));
            
                if Xge01(i,1)-Xge(i,1) > 0                                                %�Ե��������ȡ��������+1
                    Xge(i,1)=fix(Xge01(i,1))+1;
                else
                    Xge(i,1)=Xge01(i,1);
                end
                %%---�ȱü������---%%

            end
            
            %���1/2/3
            if Qtes(i,1)>=CoolingFuhe(i,1)
                Ttes(i+1,1)=(X(j,5)*Ktes*Ates*(Temperature(i,1)-Ttes(i,1))+CoolingFuhe(i,1)*3600)/(X(j,5)*mtes*Cw)+Ttes(i,1);
                Xhp02(i,1)=0;
            elseif Qtes(i,1)<CoolingFuhe(i,1) && Qtes(i,1)>0
                Qhpbui(i,1)=CoolingFuhe(i,1)-Qtes(i,1);
                Ttes(i+1,1)=Tbuc;
                
                %���2/3
                %��һ�μ����ȱ�����
                Xhp01(i,1)=Qhpbui(i,1)/Qhpsc(i,1);
                
                if Xhp01(i,1)<1
                    Xhp02(i,1)=0;
                else
                    Xhp02(i,1)=fix(Xhp01(i,1))+1;
                end
                
                Xhp(i,1)=Xhp02(i,1);
                
                %�ڶ��μ����ȱ����������ݱ߽�����,������ȱʧ���벿�ָ����µ��ȱù��ġ�COP
                if Xhp(i,1)<=X(j,4) && Xhp(i,1)>0
                    PL(i,1)=Qhpbui(i,1)/(Qhpsc(i,1)*Xhp02(i,1));                          %���1�����ڲ��ַ����µ��ȱ�
                    COPnom(i,1)=Qhpsc(i,1)/Php(i,1);
                    COPpl(i,1)=COPnom(i,1)*(-0.00006*PL(i,1)^6+0.0017*PL(i,1)^5-0.0181*PL(i,1)^4+0.096*PL(i,1)^3-0.2697*PL(i,1)^2+0.4276*PL(i,1)+0.7626);
                    Phppl(i,1)=Qhpbui(i,1)/COPpl(i,1);
                    Qhpshpl(i,1)=Qhpbui(i,1)+Phppl(i,1);
                elseif Xhp(i,1)>X(j,4)
                    Xhp02(i,1)=fix(X(j,4));
                    Qdeloss(i,1)=Qhpbui(i,1)-Qhpsc(i,1)*Xhp02(i,1);                          %Qdeloss����ȱʧ������
                    Xloss(i,1)=Qdeloss(i,1)/Qhpbui(i,1);
                    Phppl(i,1)=Php(i,1)*Xhp02(i,1);
                    Qhpshpl(i,1)=Xhp02(i,1)*Qhpsh(i,1);
                else
                    Qdeloss(i,1)=Qhpbui(i,1);                          %Qdeloss����ȱʧ������
                    Xloss(i,1)=1;
                    Phppl(i,1)=0;
                    Qhpshpl(i,1)=0;
                end                                                                      %Xhp��ʾ���յ��ȱ���������С���㣬���2 X(j,4)��ʾ���������߽�ֵ
                
                %������������
                Xge01(i,1)=Qhpshpl(i,1)/(6*0.95);
                Xge(i,1)=fix(Xge(i,1));
                
                if Xge01(i,1)-Xge(i,1) > 0                                               %�Ե��������ȡ��������+1
                    Xge(i,1)=fix(Xge01(i,1))+1;
                else
                    Xge(i,1)=Xge01(i,1);
                end
                
            end
            %%--�ȱü������--%%
            
            %%--�縺�ɼ���--%%
            Pbss(i,1)=Fuhe(i,1)+Phppl(i,1)-Ppv(i,1)*X(j,1)-Pwt(i,1)*X(j,2);
            
            %���س�ŵ�״̬�ж�
            if Pbss(i,1) < 0                                                            %pbssΪ��������Ҫ��磬��ʱ��pbss�������������޵Ĺ���
                SOC(i+1,1)=(1-0.03)*SOC(i,1)-Pbss(i,1)*0.9/(X(j,3)*Eb);
                More(i,1)=abs(Pbss(i,1));                                               %��¼����ĵ���
                elecnum(i,1)=0;
            elseif Pbss(i,1) > 0                                                        %pbssΪ��������Ҫ�ŵ�
                SOC(i+1,1)=(1-0.03)*SOC(i,1)-Pbss(i,1)/(X(j,3)*0.9*Eb);
                More(i,1)=0;                                                            %��¼������˷ѵ���
                elecnum(i,1)=abs(Pbss(i,1));
            else
                SOC(i+1,1)=SOC(i,1);
                More(i,1)=0;                                                            %��¼������˷ѵ���
                elecnum(i,1)=0;
            end
            
            %��������ж�
            if SOC(i+1,1)>=SOCmax
                SOC(i+1,1)=SOCmax;
                %�ȼ��㼴�����ʱ��������ʵ��繦��
                if SOC(i,1)~=SOCmax                                                         %~=�ǲ����ڵ���˼
                    Pbss(i,1)=(SOC(i+1,1)-(1-0.03)*SOC(i,1))*Eb*X(j,3)/(-0.9);
                else
                    Pbss(i,1)=0;
                end
                %��100%�������Ϊ����������
                if More(i,1)>=abs(Pbss(i,1))
                    More(i,1)=More(i,1)-abs(Pbss(i,1));
                else
                    More(i,1)=0;
                end
                %��¼�������
            else
                More(i,1)=0;
            end
            
            %ȱ���ж�
            if SOC(i+1,1)<=SOCmin
                SOC(i+1,1)=SOCmin;
                %���㼴������ʱ��������ʵ�ŵ繦��
                if SOC(i,1)~=SOCmin
                    Pbss(i,1)=(SOC(i+1,1)-(1-0.03)*SOC(i,1))*Eb*(-0.9)*X(j,3);
                else
                    Pbss(i,1)=0;
                end
                %��ʵȱ����
                if elecnum(i,1)>=abs(Pbss(i,1))
                    elecnum(i,1)=elecnum(i,1)-abs(Pbss(i,1));
                else
                    elecnum(i,1)=0;
                end
            else
                elecnum(i,1)=0;
            end
            
            %�⹺��ķ�
            Celebuy(i,1)=elecnum(i,1)*DianJia(i,1);
            
            %�۵�ӯ��
            Celesell(i,1)=More(i,1)*Cse;
            %Emore(i,1)=More(i,1)/(Ppv(i,1)*X(j,1)+Pwt(i,1)*X(j,2)+0.001);
            %%--�縺�ɼ������--%%
            
        end
        %%%---16.00-18.00�������---%%%
        
        
        %%%---18.00-23.00�����ϰ�ʱ�䣬�Ƿ�ֵ���---%%%
        if d(i,1)>0.78 && d(i,1)<0.97
            Ttes(i+1,1)=Ktes*Ates*(Temperature(i,1)-Ttes(i,1))/(mtes*Cw)+Ttes(i,1);
            
            %%%---�縺�ɼ���---%%%
            Pbss(i,1)=Fuhe(i,1)-Ppv(i,1)*X(j,1)-Pwt(i,1)*X(j,2);
            
            %---���س�ŵ�״̬�ж�---%
            if Pbss(i,1) < 0                                                        %pbssΪ��������Ҫ��磬��ʱ��pbss�������������޵Ĺ���
                SOC(i+1,1)=(1-0.03)*SOC(i,1)-Pbss(i,1)*0.9/(X(j,3)*Eb);
                More(i,1)=abs(Pbss(i,1));                                           %��¼����ĵ���
                elecnum(i,1)=0;
            elseif   Pbss(i,1) > 0                                                  %pbssΪ��������Ҫ�ŵ�
                SOC(i+1,1)=(1-0.03)*SOC(i,1)-Pbss(i,1)/(X(j,3)*0.9*Eb);
                More(i,1)=0;                                                        %��¼������˷ѵ���
                elecnum(i,1)=abs(Pbss(i,1));
            else
                SOC(i+1,1)=SOC(i,1);
                More(i,1)=0;                                                        %��¼������˷ѵ���
                elecnum(i,1)=0;
            end
            %---���س�ŵ�״̬�жϽ���---%
            
            %---��������ж�---%
            if  SOC(i+1,1)>=SOCmax
                SOC(i+1,1)=SOCmax;
                %�ȼ��㼴�����ʱ��������ʵ��繦��
                if  SOC(i,1)~=SOCmax                                                %~=�ǲ����ڵ���˼
                    Pbss(i,1)=(SOC(i+1,1)-(1-0.03)*SOC(i,1))*Eb*X(j,3)/(-0.9);
                else
                    Pbss(i,1)=0;
                end
                %��¼�������
                %��100%�������Ϊ����������
                if More(i,1)>=abs(Pbss(i,1))
                    More(i,1)=More(i,1)-abs(Pbss(i,1));                             %������õĺ�����
                else
                    More(i,1)=0;
                end
                %��¼�������
            else
                More(i,1)=0;
            end
            %---��������жϽ���---%
            
            %---ȱ���ж�---%
            if  SOC(i+1,1)<=SOCmin
                SOC(i+1,1)=SOCmin;
                %���㼴������ʱ��������ʵ�ŵ繦��
                if  SOC(i,1)~=SOCmin
                    Pbss(i,1)=(SOC(i+1,1)-(1-0.03)*SOC(i,1))*Eb*(-0.9)*X(j,3);
                else
                    Pbss(i,1)=0;
                end
                %�ڼ���ʱ����ʵȱ����
                if elecnum(i,1)>=abs(Pbss(i,1))
                    elecnum(i,1)=elecnum(i,1)-abs(Pbss(i,1));                       %������õĺ�����
                else
                    elecnum(i,1)=0;
                end
            else
                elecnum(i,1)=0;
            end
            
            %�⹺��ķ�
            Celebuy(i,1)=0;
            
            %�۵�ӯ��
            Celesell(i,1)=More(i,1)*Cse;
            %Emore(i,1)=More(i,1)/(Ppv(i,1)*X(j,1)+Pwt(i,1)*X(j,2)+0.001);
            
        end
        %%%---18.00-23.00�������---%%%
        
    end
    %%%---���伾�������---%%%
    %% ���ɼ�
    for i=5833:7631
        
        %%%---�縺�ɼ���---%%%
        Pbss(i,1)=Fuhe(i,1)-Ppv(i,1)*X(j,1)-Pwt(i,1)*X(j,2);
        
        %---���س�ŵ�״̬�ж�---%
        if Pbss(i,1) < 0                                                        %pbssΪ��������Ҫ��磬��ʱ��pbss�������������޵Ĺ���
            SOC(i+1,1)=(1-0.03)*SOC(i,1)-Pbss(i,1)*0.9/(X(j,3)*Eb);
            More(i,1)=abs(Pbss(i,1));                                           %��¼����ĵ���
            elecnum(i,1)=0;
        elseif   Pbss(i,1) > 0                                                  %pbssΪ��������Ҫ�ŵ�
            SOC(i+1,1)=(1-0.03)*SOC(i,1)-Pbss(i,1)/(X(j,3)*0.9*Eb);
            More(i,1)=0;                                                        %��¼������˷ѵ���
            elecnum(i,1)=abs(Pbss(i,1));
        else
            SOC(i+1,1)=SOC(i,1);
            More(i,1)=0;                                                        %��¼������˷ѵ���
            elecnum(i,1)=0;
        end
        %---���س�ŵ�״̬�жϽ���---%
        
        %---��������ж�---%
        if  SOC(i+1,1)>=SOCmax
            SOC(i+1,1)=SOCmax;
            %�ȼ��㼴�����ʱ��������ʵ��繦��
            if  SOC(i,1)~=SOCmax                                                %~=�ǲ����ڵ���˼
                Pbss(i,1)=(SOC(i+1,1)-(1-0.03)*SOC(i,1))*Eb*X(j,3)/(-0.9);
            else
                Pbss(i,1)=0;
            end
            %��¼�������
            %��100%�������Ϊ����������
            if More(i,1)>=abs(Pbss(i,1))
                More(i,1)=More(i,1)-abs(Pbss(i,1));                             %������õĺ�����
            else
                More(i,1)=0;
            end
            %��¼�������
        else
            More(i,1)=0;
        end
        %---��������жϽ���---%
        
        %---ȱ���ж�---%
        if  SOC(i+1,1)<=SOCmin
            SOC(i+1,1)=SOCmin;
            %���㼴������ʱ��������ʵ�ŵ繦��
            if  SOC(i,1)~=SOCmin
                Pbss(i,1)=(SOC(i+1,1)-(1-0.03)*SOC(i,1))*Eb*(-0.9)*X(j,3);
            else
                Pbss(i,1)=0;
            end
            %�ڼ���ʱ����ʵȱ����
            if elecnum(i,1)>=abs(Pbss(i,1))
                elecnum(i,1)=elecnum(i,1)-abs(Pbss(i,1));                       %������õĺ�����
            else
                elecnum(i,1)=0;
            end
        else
            elecnum(i,1)=0;
        end
        
        %�⹺��ķ�
        Celebuy(i,1)=elecnum(i,1)*DianJia(i,1);
        
        %�۵�ӯ��
        Celesell(i,1)=More(i,1)*Cse;
        %Emore(i,1)=More(i,1)/(Ppv(i,1)*X(j,1)+Pwt(i,1)*X(j,2)+0.001);
        
        %---ȱ���жϽ���---%
        %%%---��ģ�����---%%%
    end
    
    %% ��ů�� part02
    for i=7632:8760
        
        a(i,1)=i;
        b(i,1)=a(i,1)/24;
        c(i,1)=fix(b(i,1));
        d(i,1)=b(i,1)-c(i,1);
        
        %%%---23.00-07.00,���ʱ�䣬�ȱ���ˮ�乩��---%%%
        if d(i,1)<0.31
            
            if Ttes(i,1)<Ttessdhh                                                                                            %ˮ���¶ȿ�����ͣ�жϣ�ˮ���¶�<53
                Qhptes(i,1)=((Ttessdh-Ttes(i,1))*Cw*X(j,5)*mtes+X(j,5)*Ktes*Ates*(Ttes(i,1)-Temperature(i,1)))/3600;                %�ȱ÷��ʱ�ڹ����ȸ���
                
                %%---�ȱò��ָ��ɼ�ȱʧ����---%%
                
                %��һ�μ����ȱ�����
                Xhp01(i,1)=Qhptes(i,1)/Qhpsh(i,1);
                
                if Xhp01(i,1)<1
                    Xhp02(i,1)=0;
                else
                    Xhp02(i,1)=fix(Xhp01(i,1))+1;
                end
                
                Xhp(i,1)=Xhp02(i,1);
                
                %�ڶ��μ����ȱ����������ݱ߽�����,������ȱʧ���벿�ָ����µ��ȱù��ġ�COP
                if Xhp(i,1)<=X(j,4) && Xhp(i,1)>0
                    PL(i,1)=Qhptes(i,1)/(Qhpsh(i,1)*Xhp02(i,1));                          %���1�����ڲ��ַ����µ��ȱ�
                    COPnom(i,1)=Qhpsh(i,1)/Php(i,1);
                    COPpl(i,1)=COPnom(i,1)*(-0.00006*PL(i,1)^6+0.0017*PL(i,1)^5-0.0181*PL(i,1)^4+0.096*PL(i,1)^3-0.2697*PL(i,1)^2+0.4276*PL(i,1)+0.7626);
                    Phppl(i,1)=Qhptes(i,1)/COPpl(i,1);                                    %���㲿�ָ����¸��ȱõĹ���
                    Qhpscpl(i,1)=Qhptes(i,1)-Phppl(i,1);
                    Ttes(i+1,1)=Ttessdh;
                elseif Xhp(i,1)>X(j,4)
                    Xhp02(i,1)=fix(X(j,4));
                    %Qdeloss(i,1)=Qhptes(i,1)-Qhpsh(i,1)*Xhp02(i,1);                          %Qdeloss����ȱʧ������
                    %Xloss(i,1)=Qdeloss(i,1)/Qhptes(i,1);
                    Phppl(i,1)=Php(i,1)*Xhp02(i,1);
                    Qhpscpl(i,1)=Xhp02(i,1)*Qhpsc(i,1);
                    Ttes(i+1,1)=(Xhp02(i,1)*Qhpsh(i,1)*3600-X(j,5)*Ktes*Ates*(Ttes(i,1)-Temperature(i,1)))/(mtes*Cw*X(j,5))+Ttes(i,1);
                else
                    Phppl(i,1)=0;
                    Qhpscpl(i,1)=0;
                    Ttes(i+1,1)=Ktes*Ates*(Temperature(i,1)-Ttes(i,1))/(mtes*Cw)+Ttes(i,1);
                end                                                                      %Xhp��ʾ���յ��ȱ���������С���㣬���2 X(j,4)��ʾ���������߽�ֵ
                
                Phptes(i,1)=Phppl(i,1);
                
                %������������
                Xge01(i,1)=Qhpscpl(i,1)/(6*0.95);                                        %����ܴ���Ч������Ϊ95%
                Xge(i,1)=fix(Xge01(i,1));
                if Xge01(i,1)-Xge(i,1) > 0                                                 %�Ե��������ȡ��������+1
                    Xge(i,1)=fix(Xge01(i,1))+1;
                else
                    Xge(i,1)=Xge01(i,1);
                end
                
            else
                Ttes(i+1,1)=Ktes*Ates*(Temperature(i,1)-Ttes(i,1))/(mtes*Cw)+Ttes(i,1);
            end
            
            
            %%%---�縺�ɼ���---%%%
            Pbss(i,1)=Fuhe(i,1)-Ppv(i,1)*X(j,1)-Pwt(i,1)*X(j,2);
            
            %---���س�ŵ�״̬�ж�---%
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
            %---���س�ŵ�״̬�жϽ���---%
            
            %---��������ж�---%
            if  SOC(i+1,1)>=SOCmax
                SOC(i+1,1)=SOCmax;
                %�ȼ��㼴�����ʱ��������ʵ��繦��
                if  SOC(i,1)~=SOCmax
                    Pbss(i,1)=(SOC(i+1,1)-(1-0.03)*SOC(i,1))*Eb*X(j,3)/(-0.9);
                else
                    Pbss(i,1)=0;
                end
                %��¼�������
                %��100%�������Ϊ����������
                if More(i,1)>=abs(Pbss(i,1))
                    More(i,1)=More(i,1)-abs(Pbss(i,1));
                else
                    More(i,1)=0;
                end
                %��¼�������
            else
                More(i,1)=0;
            end
            %---��������жϽ���---%
            
            %---ȱ���ж�---%
            if  SOC(i+1,1)<=SOCmin
                SOC(i+1,1)=SOCmin;
                %���㼴������ʱ��������ʵ�ŵ繦��
                if  SOC(i,1)~=SOCmin
                    Pbss(i,1)=(SOC(i+1,1)-(1-0.03)*SOC(i,1))*Eb*(-0.9)*X(j,3);
                else
                    Pbss(i,1)=0;
                end
                %�ڼ���ʱ����ʵȱ����
                if elecnum(i,1)>=abs(Pbss(i,1))
                    elecnum(i,1)=elecnum(i,1)-abs(Pbss(i,1));
                else
                    elecnum(i,1)=0;
                end
            else
                elecnum(i,1)=0;
            end
            
            %�⹺��ķ�
            Celebuy(i,1)=Phptes(i,1)*DianJia(i,1);
            
            %�۵�ӯ��
            Celesell(i,1)=More(i,1)*Cse;
            %Emore(i,1)=More(i,1)/(Ppv(i,1)*X(j,1)+Pwt(i,1)*X(j,2)+0.001);
            
            %---ȱ���жϽ���---%
            %%%---��ģ�����---%%%
        end
        
        %%%---23.00-07.00�������---%%%
        
        %%%---07.00-08.00���㣬���ϰ�ʱ�䣬�Ƿ��ʱ��---%%%
        if d(i,1)>0.31 && d(i,1)<0.36
            
            Ttes(i+1,1)=Ktes*Ates*(Temperature(i,1)-Ttes(i,1))/(mtes*Cw)+Ttes(i,1);
            
            %%%---�縺�ɼ���---%%%
            Pbss(i,1)=Fuhe(i,1)-Ppv(i,1)*X(j,1)-Pwt(i,1)*X(j,2);
            
            %---���س�ŵ�״̬�ж�---%
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
            %---���س�ŵ�״̬�жϽ���---%
            
            %---��������ж�---%
            if  SOC(i+1,1)>=SOCmax
                SOC(i+1,1)=SOCmax;
                %�ȼ��㼴�����ʱ��������ʵ��繦��
                if  SOC(i,1)~=SOCmax
                    Pbss(i,1)=(SOC(i+1,1)-(1-0.03)*SOC(i,1))*Eb*X(j,3)/(-0.9);
                else
                    Pbss(i,1)=0;
                end
                %��¼�������
                %��100%�������Ϊ����������
                if More(i,1)>=abs(Pbss(i,1))
                    More(i,1)=More(i,1)-abs(Pbss(i,1));
                else
                    More(i,1)=0;
                end
                %��¼�������
            else
                More(i,1)=0;
            end
            %---��������жϽ���---%
            
            %---ȱ���ж�---%
            if  SOC(i+1,1)<=SOCmin
                SOC(i+1,1)=SOCmin;
                %���㼴������ʱ��������ʵ�ŵ繦��
                if  SOC(i,1)~=SOCmin
                    Pbss(i,1)=(SOC(i+1,1)-(1-0.03)*SOC(i,1))*Eb*(-0.9)*X(j,3);
                else
                    Pbss(i,1)=0;
                end
                %�ڼ���ʱ����ʵȱ����
                if elecnum(i,1)>=abs(Pbss(i,1))
                    elecnum(i,1)=elecnum(i,1)-abs(Pbss(i,1));
                else
                    elecnum(i,1)=0;
                end
            else
                elecnum(i,1)=0;
            end
            
            %�⹺��ķ�
            Celebuy(i,1)=0;
            
            %�۵�ӯ��
            Celesell(i,1)=More(i,1)*Cse;
            %Emore(i,1)=More(i,1)/(Ppv(i,1)*X(j,1)+Pwt(i,1)*X(j,2)+0.001);
            
        end
        %%%---07.00-08.00�������---%%%
        
        %%%---8.00-9.00,����ʱ�䣬ƽ�����---%%%
        if d(i,1)>0.36 && d(i,1)<0.39
            
            %%--ˮ���¶ȼ���--%%
            Ttes(i+1,1)=Ktes*Ates*(Temperature(i,1)-Ttes(i,1))/(mtes*Cw)+Ttes(i,1);
            %%--ˮ���¶ȼ������--%%
            
            %%--�ȱü���--%%
            %��һ�μ����ȱ�����
            Xhp01(i,1)=HeatingFuhe(i,1)/Qhpsh(i,1);
            
            if Xhp01(i,1)<1
                Xhp02(i,1)=0;
            else
                Xhp02(i,1)=fix(Xhp01(i,1))+1;
            end
            
            Xhp(i,1)=Xhp02(i,1);
                
            %�ڶ��μ����ȱ����������ݱ߽�����,������ȱʧ���벿�ָ����µ��ȱù��ġ�COP
            if Xhp(i,1)<=X(j,4) && Xhp(i,1)>0
                PL(i,1)=HeatingFuhe(i,1)/(Qhpsh(i,1)*Xhp02(i,1));                          %���1�����ڲ��ַ����µ��ȱ�
                COPnom(i,1)=Qhpsh(i,1)/Php(i,1);
                COPpl(i,1)=COPnom(i,1)*(-0.00006*PL(i,1)^6+0.0017*PL(i,1)^5-0.0181*PL(i,1)^4+0.096*PL(i,1)^3-0.2697*PL(i,1)^2+0.4276*PL(i,1)+0.7626);
                Phppl(i,1)=HeatingFuhe(i,1)/COPpl(i,1);
                Qhpscpl(i,1)=HeatingFuhe(i,1)-Phppl(i,1);
            elseif Xhp(i,1)>X(j,4)
                Xhp02(i,1)=fix(X(j,4));
                Qdeloss(i,1)=HeatingFuhe(i,1)-Qhpsh(i,1)*Xhp02(i,1);                          %Qdeloss����ȱʧ������
                Xloss(i,1)=Qdeloss(i,1)/HeatingFuhe(i,1);
                Phppl(i,1)=Php(i,1)*Xhp02(i,1);
                Qhpscpl(i,1)=Xhp02(i,1)*Qhpsc(i,1);
            else
                Qdeloss(i,1)=HeatingFuhe(i,1);                          %Qdeloss����ȱʧ������
                Xloss(i,1)=1;
                Phppl(i,1)=0;
                Qhpscpl(i,1)=0;
            end                                                                      %Xhp��ʾ���յ��ȱ���������С���㣬���2 X(j,4)��ʾ���������߽�ֵ
            
            %�������ܵ�����
            Xge01(i,1)=Qhpscpl(i,1)/(6*0.95);
            Xge(i,1)=fix(Xge(i,1));
            
            if Xge01(i,1)-Xge(i,1) > 0                                                %�Ե��������ȡ��������+1
                Xge(i,1)=fix(Xge01(i,1))+1;
            else
                Xge(i,1)=Xge01(i,1);
            end
            %%--�ȱü������--%%
            
            %%--�縺�ɼ���--%%
            Pbss(i,1)=Fuhe(i,1)+Phppl(i,1)-Ppv(i,1)*X(j,1)-Pwt(i,1)*X(j,2);
            
            %���س�ŵ�״̬�ж�
            if Pbss(i,1) < 0                                                         %pbssΪ��������Ҫ��磬��ʱ��pbss�������������޵Ĺ���
                SOC(i+1,1)=(1-0.03)*SOC(i,1)-Pbss(i,1)*0.9/(X(j,3)*Eb);
                More(i,1)=abs(Pbss(i,1));                                            %��¼����ĵ���
                elecnum(i,1)=0;
            elseif Pbss(i,1) > 0                                                     %pbssΪ��������Ҫ�ŵ�
                SOC(i+1,1)=(1-0.03)*SOC(i,1)-Pbss(i,1)/(X(j,3)*0.9*Eb);
                More(i,1)=0;                                                         %��¼������˷ѵ���
                elecnum(i,1)=abs(Pbss(i,1));
            else
                SOC(i+1,1)=SOC(i,1);
                More(i,1)=0;                                                         %��¼������˷ѵ���
                elecnum(i,1)=0;
            end
            
            %��������ж�
            if SOC(i+1,1)>=SOCmax
                SOC(i+1,1)=SOCmax;
                %�ȼ��㼴�����ʱ��������ʵ��繦��
                if SOC(i,1)~=SOCmax                                                      %~=�ǲ����ڵ���˼
                    Pbss(i,1)=(SOC(i+1,1)-(1-0.03)*SOC(i,1))*Eb*X(j,3)/(-0.9);
                else
                    Pbss(i,1)=0;
                end
                %��100%�������Ϊ����������
                if More(i,1)>=abs(Pbss(i,1))
                    More(i,1)=More(i,1)-abs(Pbss(i,1));
                else
                    More(i,1)=0;
                end
            else
                More(i,1)=0;
            end
            
            %ȱ���ж�
            if SOC(i+1,1)<=SOCmin
                SOC(i+1,1)=SOCmin;
                %���㼴������ʱ��������ʵ�ŵ繦��
                if SOC(i,1)~=SOCmin
                    Pbss(i,1)=(SOC(i+1,1)-(1-0.03)*SOC(i,1))*Eb*(-0.9)*X(j,3);
                else
                    Pbss(i,1)=0;
                end
                %��ʵȱ����
                if elecnum(i,1)>=abs(Pbss(i,1))
                    elecnum(i,1)=elecnum(i,1)-abs(Pbss(i,1));
                else
                    elecnum(i,1)=0;
                end
            else
                elecnum(i,1)=0;
            end
            
            %�⹺��ķ�
            Celebuy(i,1)=elecnum(i,1)*DianJia(i,1);
            
            %�۵�ӯ��
            Celesell(i,1)=More(i,1)*Cse;
            %Emore(i,1)=More(i,1)/(Ppv(i,1)*X(j,1)+Pwt(i,1)*X(j,2)+0.001);
            %%--�縺�ɼ������--%%
            
        end
        %%%---8.00-9.00�������---%%%
        
        %%%---9.00-12.00,����ʱ�䣬��ֵ��ۣ�����ˮ��---%%%
        if d(i,1)>0.39 && d(i,1)<0.51
            
            %%--�ȱü���--%%
            Qtes(i,1)=((Ttes(i,1)-Tbui)*Cw*X(j,5)*mtes-X(j,5)*Ktes*Ates*(Ttes(i,1)-Temperature(i,1)))/3600;
            
            %��ֹ���������ȱø���
            if Qtes(i,1)<= 0
                Qtes(i,1)=0;
                
                %%--ˮ���¶ȼ���--%%
                Ttes(i+1,1)=Ktes*Ates*(Temperature(i,1)-Ttes(i,1))/(mtes*Cw)+Ttes(i,1);
                %%--ˮ���¶ȼ������--%%
            
                %%--�ȱü���--%%
                %��һ�μ����ȱ�����
                Xhp01(i,1)=HeatingFuhe(i,1)/Qhpsh(i,1);
                
                if Xhp01(i,1)<1
                    Xhp02(i,1)=0;
                else
                    Xhp02(i,1)=fix(Xhp01(i,1))+1;
                end
                
                Xhp(i,1)=Xhp02(i,1);
                
                %�ڶ��μ����ȱ����������ݱ߽�����,������ȱʧ���벿�ָ����µ��ȱù��ġ�COP
                if Xhp(i,1)<=X(j,4) && Xhp(i,1)>0
                    PL(i,1)=HeatingFuhe(i,1)/(Qhpsh(i,1)*Xhp02(i,1));                          %���1�����ڲ��ַ����µ��ȱ�
                    COPnom(i,1)=Qhpsh(i,1)/Php(i,1);
                    COPpl(i,1)=COPnom(i,1)*(-0.00006*PL(i,1)^6+0.0017*PL(i,1)^5-0.0181*PL(i,1)^4+0.096*PL(i,1)^3-0.2697*PL(i,1)^2+0.4276*PL(i,1)+0.7626);
                    Phppl(i,1)=HeatingFuhe(i,1)/COPpl(i,1);
                    Qhpscpl(i,1)=HeatingFuhe(i,1)-Phppl(i,1);
                elseif Xhp(i,1)>X(j,4)
                    Xhp02(i,1)=fix(X(j,4));
                    Qdeloss(i,1)=HeatingFuhe(i,1)-Qhpsh(i,1)*Xhp02(i,1);                          %Qdeloss����ȱʧ������
                    Xloss(i,1)=Qdeloss(i,1)/HeatingFuhe(i,1);
                    Phppl(i,1)=Php(i,1)*Xhp02(i,1);
                    Qhpscpl(i,1)=Xhp02(i,1)*Qhpsc(i,1);
                else
                    Qdeloss(i,1)=HeatingFuhe(i,1);                          %Qdeloss����ȱʧ������
                    Xloss(i,1)=1;
                    Phppl(i,1)=0;
                    Qhpscpl(i,1)=0;
                end                                                                      %Xhp��ʾ���յ��ȱ���������С���㣬���2 X(j,4)��ʾ���������߽�ֵ
            
                %�������ܵ�����
                Xge01(i,1)=Qhpscpl(i,1)/(6*0.95);
                Xge(i,1)=fix(Xge(i,1));
            
                if Xge01(i,1)-Xge(i,1) > 0                                                %�Ե��������ȡ��������+1
                    Xge(i,1)=fix(Xge01(i,1))+1;
                else
                    Xge(i,1)=Xge01(i,1);
                end
                %%--�ȱü������--%%
            end
            
            %���1/2/3
            if Qtes(i,1) >= HeatingFuhe(i,1)
                Ttes(i+1,1)=(X(j,5)*Ktes*Ates*(Temperature(i,1)-Ttes(i,1))-HeatingFuhe(i,1)*3600)/(X(j,5)*mtes*Cw)+Ttes(i,1);
                Xhp02(i,1)=0;
            elseif Qtes(i,1)<HeatingFuhe(i,1) && Qtes(i,1)>0
                Qhpbui(i,1)=HeatingFuhe(i,1)-Qtes(i,1);
                Ttes(i+1,1)=Tbui;
                
                %���2/3
                %��һ�μ����ȱ�����
                Xhp01(i,1)=Qhpbui(i,1)/Qhpsh(i,1);
                
                if Xhp01(i,1)<1
                    Xhp02(i,1)=0;
                else
                    Xhp02(i,1)=fix(Xhp01(i,1))+1;
                end
                
                Xhp(i,1)=Xhp02(i,1);
                
                %�ڶ��μ����ȱ����������ݱ߽�����,������ȱʧ���벿�ָ����µ��ȱù��ġ�COP
                if Xhp(i,1) <= X(j,4) && Xhp(i,1)>0
                    PL(i,1)=Qhpbui(i,1)/(Qhpsh(i,1)*Xhp02(i,1));                          %���1�����ڲ��ַ����µ��ȱ�
                    COPnom(i,1)=Qhpsh(i,1)/Php(i,1);
                    COPpl(i,1)=COPnom(i,1)*(-0.00006*PL(i,1)^6+0.0017*PL(i,1)^5-0.0181*PL(i,1)^4+0.096*PL(i,1)^3-0.2697*PL(i,1)^2+0.4276*PL(i,1)+0.7626);
                    Phppl(i,1)=Qhpbui(i,1)/COPpl(i,1);
                    Qhpscpl(i,1)=Qhpbui(i,1)-Phppl(i,1);
                elseif Xhp(i,1)> X(j,4)
                    Xhp02(i,1)=fix(X(j,4));
                    Qdeloss(i,1)=Qhpbui(i,1)-Qhpsh(i,1)*Xhp02(i,1);                          %Qdeloss����ȱʧ������
                    Xloss(i,1)=Qdeloss(i,1)/Qhpbui(i,1);
                    Phppl(i,1)=Php(i,1)*Xhp02(i,1);
                    Qhpscpl(i,1)=Xhp02(i,1)*Qhpsc(i,1);
                else
                    Qdeloss(i,1)=Qhpbui(i,1);                          %Qdeloss����ȱʧ������
                    Xloss(i,1)=1;
                    Phppl(i,1)=0;
                    Qhpscpl(i,1)=0;
                end                                                                      %Xhp��ʾ���յ��ȱ���������С���㣬���2 X(j,4)��ʾ���������߽�ֵ
                
                %������������
                Xge01(i,1)=Qhpscpl(i,1)/(6*0.95);
                Xge(i,1)=fix(Xge(i,1));
                
                if Xge01(i,1)-Xge(i,1) > 0                                               %�Ե��������ȡ��������+1
                    Xge(i,1)=fix(Xge01(i,1))+1;
                else
                    Xge(i,1)=Xge01(i,1);
                end
                
            end
            %%--�ȱü������--%%
            
            %%--�縺�ɼ���--%%
            Pbss(i,1)=Fuhe(i,1)+Phppl(i,1)-Ppv(i,1)*X(j,1)-Pwt(i,1)*X(j,2);
            
            %���س�ŵ�״̬�ж�
            if Pbss(i,1) < 0                                                            %pbssΪ��������Ҫ��磬��ʱ��pbss�������������޵Ĺ���
                SOC(i+1,1)=(1-0.03)*SOC(i,1)-Pbss(i,1)*0.9/(X(j,3)*Eb);
                More(i,1)=abs(Pbss(i,1));                                               %��¼����ĵ���
                elecnum(i,1)=0;
            elseif Pbss(i,1) > 0                                                        %pbssΪ��������Ҫ�ŵ�
                SOC(i+1,1)=(1-0.03)*SOC(i,1)-Pbss(i,1)/(X(j,3)*0.9*Eb);
                More(i,1)=0;                                                            %��¼������˷ѵ���
                elecnum(i,1)=abs(Pbss(i,1));
            else
                SOC(i+1,1)=SOC(i,1);
                More(i,1)=0;                                                            %��¼������˷ѵ���
                elecnum(i,1)=0;
            end
            
            %��������ж�
            if SOC(i+1,1)>=SOCmax
                SOC(i+1,1)=SOCmax;
                %�ȼ��㼴�����ʱ��������ʵ��繦��
                if SOC(i,1)~=SOCmax                                                         %~=�ǲ����ڵ���˼
                    Pbss(i,1)=(SOC(i+1,1)-(1-0.03)*SOC(i,1))*Eb*X(j,3)/(-0.9);
                else
                    Pbss(i,1)=0;
                end
                %��100%�������Ϊ����������
                if More(i,1)>=abs(Pbss(i,1))
                    More(i,1)=More(i,1)-abs(Pbss(i,1));
                else
                    More(i,1)=0;
                end
                %��¼�������
            else
                More(i,1)=0;
            end
            
            %ȱ���ж�
            if SOC(i+1,1)<=SOCmin
                SOC(i+1,1)=SOCmin;
                %���㼴������ʱ��������ʵ�ŵ繦��
                if SOC(i,1)~=SOCmin
                    Pbss(i,1)=(SOC(i+1,1)-(1-0.03)*SOC(i,1))*Eb*(-0.9)*X(j,3);
                else
                    Pbss(i,1)=0;
                end
                %��ʵȱ����
                if elecnum(i,1)>=abs(Pbss(i,1))
                    elecnum(i,1)=elecnum(i,1)-abs(Pbss(i,1));
                else
                    elecnum(i,1)=0;
                end
            else
                elecnum(i,1)=0;
            end
            
            %�⹺��ķ�
            Celebuy(i,1)=elecnum(i,1)*DianJia(i,1);
            
            %�۵�ӯ��
            Celesell(i,1)=More(i,1)*Cse;
            %Emore(i,1)=More(i,1)/(Ppv(i,1)*X(j,1)+Pwt(i,1)*X(j,2)+0.001);
            %%--�縺�ɼ������--%%
            
        end
        %%%---9.00-12.00�������---%%%
        
        %%%---12.00-16.00������ʱ�䣬ƽ�����---%%%
        if d(i,1)>0.51 && d(i,1)<0.69
            
            %%--ˮ���¶ȼ���--%%
            Ttes(i+1,1)=Ktes*Ates*(Temperature(i,1)-Ttes(i,1))/(mtes*Cw)+Ttes(i,1);
            %%--ˮ���¶ȼ������--%%
            
            %%--�ȱü���--%%
            %��һ�μ����ȱ�����
            Xhp01(i,1)=HeatingFuhe(i,1)/Qhpsh(i,1);
            
            if Xhp01(i,1)<1
                Xhp02(i,1)=0;
            else
                Xhp02(i,1)=fix(Xhp01(i,1))+1;
            end
            
            Xhp(i,1)=Xhp02(i,1);
                
            %�ڶ��μ����ȱ����������ݱ߽�����,������ȱʧ���벿�ָ����µ��ȱù��ġ�COP
            if Xhp(i,1)<=X(j,4) && Xhp(i,1)>0
                PL(i,1)=HeatingFuhe(i,1)/(Qhpsh(i,1)*Xhp02(i,1));                          %���1�����ڲ��ַ����µ��ȱ�
                COPnom(i,1)=Qhpsh(i,1)/Php(i,1);
                COPpl(i,1)=COPnom(i,1)*(-0.00006*PL(i,1)^6+0.0017*PL(i,1)^5-0.0181*PL(i,1)^4+0.096*PL(i,1)^3-0.2697*PL(i,1)^2+0.4276*PL(i,1)+0.7626);
                Phppl(i,1)=HeatingFuhe(i,1)/COPpl(i,1);
                Qhpscpl(i,1)=HeatingFuhe(i,1)-Phppl(i,1);
            elseif Xhp(i,1)>X(j,4)
                Xhp02(i,1)=fix(X(j,4));
                Qdeloss(i,1)=HeatingFuhe(i,1)-Qhpsh(i,1)*Xhp02(i,1);                          %Qdeloss����ȱʧ������
                Xloss(i,1)=Qdeloss(i,1)/HeatingFuhe(i,1);
                Phppl(i,1)=Php(i,1)*Xhp02(i,1);
                Qhpscpl(i,1)=Xhp02(i,1)*Qhpsc(i,1);
            else
                Qdeloss(i,1)=HeatingFuhe(i,1);                          %Qdeloss����ȱʧ������
                Xloss(i,1)=1;
                Phppl(i,1)=0;
                Qhpscpl(i,1)=0;
            end                                                                      %Xhp��ʾ���յ��ȱ���������С���㣬���2 X(j,4)��ʾ���������߽�ֵ
            
            %�������ܵ�����
            Xge01(i,1)=Qhpscpl(i,1)/(6*0.95);
            Xge(i,1)=fix(Xge(i,1));
            
            if Xge01(i,1)-Xge(i,1) > 0                                                %�Ե��������ȡ��������+1
                Xge(i,1)=fix(Xge01(i,1))+1;
            else
                Xge(i,1)=Xge01(i,1);
            end
            %%--�ȱü������--%%
            
            %%--�縺�ɼ���--%%
            Pbss(i,1)=Fuhe(i,1)+Phppl(i,1)-Ppv(i,1)*X(j,1)-Pwt(i,1)*X(j,2);
            
            %���س�ŵ�״̬�ж�
            if Pbss(i,1) < 0                                                         %pbssΪ��������Ҫ��磬��ʱ��pbss�������������޵Ĺ���
                SOC(i+1,1)=(1-0.03)*SOC(i,1)-Pbss(i,1)*0.9/(X(j,3)*Eb);
                More(i,1)=abs(Pbss(i,1));                                            %��¼����ĵ���
                elecnum(i,1)=0;
            elseif Pbss(i,1) > 0                                                     %pbssΪ��������Ҫ�ŵ�
                SOC(i+1,1)=(1-0.03)*SOC(i,1)-Pbss(i,1)/(X(j,3)*0.9*Eb);
                More(i,1)=0;                                                         %��¼������˷ѵ���
                elecnum(i,1)=abs(Pbss(i,1));
            else
                SOC(i+1,1)=SOC(i,1);
                More(i,1)=0;                                                         %��¼������˷ѵ���
                elecnum(i,1)=0;
            end
            
            %��������ж�
            if SOC(i+1,1)>=SOCmax
                SOC(i+1,1)=SOCmax;
                %�ȼ��㼴�����ʱ��������ʵ��繦��
                if SOC(i,1)~=SOCmax                                                      %~=�ǲ����ڵ���˼
                    Pbss(i,1)=(SOC(i+1,1)-(1-0.03)*SOC(i,1))*Eb*X(j,3)/(-0.9);
                else
                    Pbss(i,1)=0;
                end
                %��100%�������Ϊ����������
                if More(i,1)>=abs(Pbss(i,1))
                    More(i,1)=More(i,1)-abs(Pbss(i,1));
                else
                    More(i,1)=0;
                end
                %��¼�������
            else
                More(i,1)=0;
            end
            
            %ȱ���ж�
            if SOC(i+1,1)<=SOCmin
                SOC(i+1,1)=SOCmin;
                %���㼴������ʱ��������ʵ�ŵ繦��
                if SOC(i,1)~=SOCmin
                    Pbss(i,1)=(SOC(i+1,1)-(1-0.03)*SOC(i,1))*Eb*(-0.9)*X(j,3);
                else
                    Pbss(i,1)=0;
                end
                %��ʵȱ����
                if elecnum(i,1)>=abs(Pbss(i,1))
                    elecnum(i,1)=elecnum(i,1)-abs(Pbss(i,1));
                else
                    elecnum(i,1)=0;
                end
            else
                elecnum(i,1)=0;
            end
            
            %�⹺��ķ�
            Celebuy(i,1)=elecnum(i,1)*DianJia(i,1);
            
            %�۵�ӯ��
            Celesell(i,1)=More(i,1)*Cse;
            %Emore(i,1)=More(i,1)/(Ppv(i,1)*X(j,1)+Pwt(i,1)*X(j,2)+0.001);
            %%--�縺�ɼ������--%%
            
        end
        %%%---12.00-16.00�������---%%%
        
        %%%---16.00-18.00,����ʱ�䣬��ֵ��ۣ�����ˮ��---%%%
        if d(i,1)>0.69 && d(i,1)<0.78
            
            %%--�ȱü���--%%
            Qtes(i,1)=((Ttes(i,1)-Tbui)*Cw*X(j,5)*mtes-X(j,5)*Ktes*Ates*(Ttes(i,1)-Temperature(i,1)))/3600;
            
            %��ֹ���������ȱø���
            if Qtes(i,1)<0
                Qtes(i,1)=0;
                
                %%--ˮ���¶ȼ���--%%
                Ttes(i+1,1)=Ktes*Ates*(Temperature(i,1)-Ttes(i,1))/(mtes*Cw)+Ttes(i,1);
                %%--ˮ���¶ȼ������--%%
            
                %%--�ȱü���--%%
                %��һ�μ����ȱ�����
                Xhp01(i,1)=HeatingFuhe(i,1)/Qhpsh(i,1);
                
                if Xhp01(i,1)<1
                    Xhp02(i,1)=0;
                else
                    Xhp02(i,1)=fix(Xhp01(i,1))+1;
                end
                
                Xhp(i,1)=Xhp02(i,1);
                
                %�ڶ��μ����ȱ����������ݱ߽�����,������ȱʧ���벿�ָ����µ��ȱù��ġ�COP
                if Xhp(i,1)<=X(j,4) && Xhp(i,1)>0
                    PL(i,1)=HeatingFuhe(i,1)/(Qhpsh(i,1)*Xhp02(i,1));                          %���1�����ڲ��ַ����µ��ȱ�
                    COPnom(i,1)=Qhpsh(i,1)/Php(i,1);
                    COPpl(i,1)=COPnom(i,1)*(-0.00006*PL(i,1)^6+0.0017*PL(i,1)^5-0.0181*PL(i,1)^4+0.096*PL(i,1)^3-0.2697*PL(i,1)^2+0.4276*PL(i,1)+0.7626);
                    Phppl(i,1)=HeatingFuhe(i,1)/COPpl(i,1);
                    Qhpscpl(i,1)=HeatingFuhe(i,1)-Phppl(i,1);
                elseif Xhp(i,1)>X(j,4)
                    Xhp02(i,1)=fix(X(j,4));
                    Qdeloss(i,1)=HeatingFuhe(i,1)-Qhpsh(i,1)*Xhp02(i,1);                          %Qdeloss����ȱʧ������
                    Xloss(i,1)=Qdeloss(i,1)/HeatingFuhe(i,1);
                    Phppl(i,1)=Php(i,1)*Xhp02(i,1);
                    Qhpscpl(i,1)=Xhp02(i,1)*Qhpsc(i,1);
                else
                    Qdeloss(i,1)=HeatingFuhe(i,1);                          %Qdeloss����ȱʧ������
                    Xloss(i,1)=1;
                    Phppl(i,1)=0;
                    Qhpscpl(i,1)=0;
                end                                                                      %Xhp��ʾ���յ��ȱ���������С���㣬���2 X(j,4)��ʾ���������߽�ֵ
            
                %�������ܵ�����
                Xge01(i,1)=Qhpscpl(i,1)/(6*0.95);
                Xge(i,1)=fix(Xge(i,1));
            
                if Xge01(i,1)-Xge(i,1) > 0                                                %�Ե��������ȡ��������+1
                    Xge(i,1)=fix(Xge01(i,1))+1;
                else
                    Xge(i,1)=Xge01(i,1);
                end
                %%--�ȱü������--%%
            end
            
            %���1/2/3
            if Qtes(i,1)>=HeatingFuhe(i,1)
                Ttes(i+1,1)=(X(j,5)*Ktes*Ates*(Temperature(i,1)-Ttes(i,1))-HeatingFuhe(i,1)*3600)/(X(j,5)*mtes*Cw)+Ttes(i,1);
                Xhp02(i,1)=0;
            elseif Qtes(i,1)<HeatingFuhe(i,1) && Qtes(i,1)>0
                Qhpbui(i,1)=HeatingFuhe(i,1)-Qtes(i,1);
                Ttes(i+1,1)=Tbui;
                
                %���2/3
                %��һ�μ����ȱ�����
                Xhp01(i,1)=Qhpbui(i,1)/Qhpsh(i,1);
                
                if Xhp01(i,1)<1
                    Xhp02(i,1)=0;
                else
                    Xhp02(i,1)=fix(Xhp01(i,1))+1;
                end
                
                Xhp(i,1)=Xhp02(i,1);
                
                %�ڶ��μ����ȱ����������ݱ߽�����,������ȱʧ���벿�ָ����µ��ȱù��ġ�COP
                if Xhp(i,1) <= X(j,4) && Xhp(i,1)>0
                    PL(i,1)=Qhpbui(i,1)/(Qhpsh(i,1)*Xhp02(i,1));                          %���1�����ڲ��ַ����µ��ȱ�
                    COPnom(i,1)=Qhpsh(i,1)/Php(i,1);
                    COPpl(i,1)=COPnom(i,1)*(-0.00006*PL(i,1)^6+0.0017*PL(i,1)^5-0.0181*PL(i,1)^4+0.096*PL(i,1)^3-0.2697*PL(i,1)^2+0.4276*PL(i,1)+0.7626);
                    Phppl(i,1)=Qhpbui(i,1)/COPpl(i,1);
                    Qhpscpl(i,1)=Qhpbui(i,1)-Phppl(i,1);
                elseif Xhp(i,1)> X(j,4)
                    Xhp02(i,1)=fix(X(j,4));
                    Qdeloss(i,1)=Qhpbui(i,1)-Qhpsh(i,1)*Xhp02(i,1);                          %Qdeloss����ȱʧ������
                    Xloss(i,1)=Qdeloss(i,1)/Qhpbui(i,1);
                    Phppl(i,1)=Php(i,1)*Xhp02(i,1);
                    Qhpscpl(i,1)=Xhp02(i,1)*Qhpsc(i,1);
                else
                    Qdeloss(i,1)=Qhpbui(i,1);                          %Qdeloss����ȱʧ������
                    Xloss(i,1)=1;
                    Phppl(i,1)=0;
                    Qhpscpl(i,1)=0;
                end                                                                      %Xhp��ʾ���յ��ȱ���������С���㣬���2 X(j,4)��ʾ���������߽�ֵ
                
                %������������
                Xge01(i,1)=Qhpscpl(i,1)/(6*0.95);
                Xge(i,1)=fix(Xge(i,1));
                
                if Xge01(i,1)-Xge(i,1) > 0                                               %�Ե��������ȡ��������+1
                    Xge(i,1)=fix(Xge01(i,1))+1;
                else
                    Xge(i,1)=Xge01(i,1);
                end
                
            end
            %%--�ȱü������--%%
            
            %%--�縺�ɼ���--%%
            Pbss(i,1)=Fuhe(i,1)+Phppl(i,1)-Ppv(i,1)*X(j,1)-Pwt(i,1)*X(j,2);
            
            %���س�ŵ�״̬�ж�
            if Pbss(i,1) < 0                                                            %pbssΪ��������Ҫ��磬��ʱ��pbss�������������޵Ĺ���
                SOC(i+1,1)=(1-0.03)*SOC(i,1)-Pbss(i,1)*0.9/(X(j,3)*Eb);
                More(i,1)=abs(Pbss(i,1));                                               %��¼����ĵ���
                elecnum(i,1)=0;
            elseif Pbss(i,1) > 0                                                        %pbssΪ��������Ҫ�ŵ�
                SOC(i+1,1)=(1-0.03)*SOC(i,1)-Pbss(i,1)/(X(j,3)*0.9*Eb);
                More(i,1)=0;                                                            %��¼������˷ѵ���
                elecnum(i,1)=abs(Pbss(i,1));
            else
                SOC(i+1,1)=SOC(i,1);
                More(i,1)=0;                                                            %��¼������˷ѵ���
                elecnum(i,1)=0;
            end
            
            %��������ж�
            if SOC(i+1,1)>=SOCmax
                SOC(i+1,1)=SOCmax;
                %�ȼ��㼴�����ʱ��������ʵ��繦��
                if SOC(i,1)~=SOCmax                                                         %~=�ǲ����ڵ���˼
                    Pbss(i,1)=(SOC(i+1,1)-(1-0.03)*SOC(i,1))*Eb*X(j,3)/(-0.9);
                else
                    Pbss(i,1)=0;
                end
                %��100%�������Ϊ����������
                if More(i,1)>=abs(Pbss(i,1))
                    More(i,1)=More(i,1)-abs(Pbss(i,1));
                else
                    More(i,1)=0;
                end
                %��¼�������
            else
                More(i,1)=0;
            end
            
            %ȱ���ж�
            if SOC(i+1,1)<=SOCmin
                SOC(i+1,1)=SOCmin;
                %���㼴������ʱ��������ʵ�ŵ繦��
                if SOC(i,1)~=SOCmin
                    Pbss(i,1)=(SOC(i+1,1)-(1-0.03)*SOC(i,1))*Eb*(-0.9)*X(j,3);
                else
                    Pbss(i,1)=0;
                end
                %��ʵȱ����
                if elecnum(i,1)>=abs(Pbss(i,1))
                    elecnum(i,1)=elecnum(i,1)-abs(Pbss(i,1));
                else
                    elecnum(i,1)=0;
                end
            else
                elecnum(i,1)=0;
            end
            
            %�⹺��ķ�
            Celebuy(i,1)=elecnum(i,1)*DianJia(i,1);
            
            %�۵�ӯ��
            Celesell(i,1)=More(i,1)*Cse;
            %Emore(i,1)=More(i,1)/(Ppv(i,1)*X(j,1)+Pwt(i,1)*X(j,2)+0.001);
            %%--�縺�ɼ������--%%
            
        end
        %%%---16.00-18.00�������---%%%
        
        
        %%%---18.00-23.00�����ϰ�ʱ�䣬�Ƿ�ȵ��---%%%
        if d(i,1)>0.78 && d(i,1)<0.97
            Ttes(i+1,1)=Ktes*Ates*(Temperature(i,1)-Ttes(i,1))/(mtes*Cw)+Ttes(i,1);
            
            %%%---�縺�ɼ���---%%%
            Pbss(i,1)=Fuhe(i,1)-Ppv(i,1)*X(j,1)-Pwt(i,1)*X(j,2);
            
            %---���س�ŵ�״̬�ж�---%
            if Pbss(i,1) < 0                                                        %pbssΪ��������Ҫ��磬��ʱ��pbss�������������޵Ĺ���
                SOC(i+1,1)=(1-0.03)*SOC(i,1)-Pbss(i,1)*0.9/(X(j,3)*Eb);
                More(i,1)=abs(Pbss(i,1));                                           %��¼����ĵ���
                elecnum(i,1)=0;
            elseif   Pbss(i,1) > 0                                                  %pbssΪ��������Ҫ�ŵ�
                SOC(i+1,1)=(1-0.03)*SOC(i,1)-Pbss(i,1)/(X(j,3)*0.9*Eb);
                More(i,1)=0;                                                        %��¼������˷ѵ���
                elecnum(i,1)=abs(Pbss(i,1));
            else
                SOC(i+1,1)=SOC(i,1);
                More(i,1)=0;                                                        %��¼������˷ѵ���
                elecnum(i,1)=0;
            end
            %---���س�ŵ�״̬�жϽ���---%
            
            %---��������ж�---%
            if  SOC(i+1,1)>=SOCmax
                SOC(i+1,1)=SOCmax;
                %�ȼ��㼴�����ʱ��������ʵ��繦��
                if  SOC(i,1)~=SOCmax                                                %~=�ǲ����ڵ���˼
                    Pbss(i,1)=(SOC(i+1,1)-(1-0.03)*SOC(i,1))*Eb*X(j,3)/(-0.9);
                else
                    Pbss(i,1)=0;
                end
                %��¼�������
                %��100%�������Ϊ����������
                if More(i,1)>=abs(Pbss(i,1))
                    More(i,1)=More(i,1)-abs(Pbss(i,1));                             %������õĺ�����
                else
                    More(i,1)=0;
                end
                %��¼�������
            else
                More(i,1)=0;
            end
            %---��������жϽ���---%
            
            %---ȱ���ж�---%
            if  SOC(i+1,1)<=SOCmin
                SOC(i+1,1)=SOCmin;
                %���㼴������ʱ��������ʵ�ŵ繦��
                if  SOC(i,1)~=SOCmin
                    Pbss(i,1)=(SOC(i+1,1)-(1-0.03)*SOC(i,1))*Eb*(-0.9)*X(j,3);
                else
                    Pbss(i,1)=0;
                end
                %�ڼ���ʱ����ʵȱ����
                if elecnum(i,1)>=abs(Pbss(i,1))
                    elecnum(i,1)=elecnum(i,1)-abs(Pbss(i,1));                       %������õĺ�����
                else
                    elecnum(i,1)=0;
                end
            else
                elecnum(i,1)=0;
            end
            
            %�⹺��ķ�
            Celebuy(i,1)=0;
            
            %�۵�ӯ��
            Celesell(i,1)=More(i,1)*Cse;
            %Emore(i,1)=More(i,1)/(Ppv(i,1)*X(j,1)+Pwt(i,1)*X(j,2)+0.001);
            
            %---ȱ���жϽ���---%
            %%%---��ģ�����---%%%
        end
        %%%---18.00-23.00�������---%%%
        
        %%%---�˼��ڼ������---%%%
    end
    
    
    %Ŀ�꺯������
    %Chengfa(j,1)=sum(LPS)*0.5*20;                                                          %1000:�ͷ�ϵ��k��   20����ѵ��ۣ�
    Xgem=max(Xge);
    LCC1(j,1)=Cpv*X(j,1)+Cwt*X(j,2)+Cbss*X(j,3)+Chp*X(j,4)+Cge*Xgem+Ctes*X(j,5);            %ϵͳ����
    Cs(j,1)=(sum(Celebuy)-sum(Celesell))*60;                                                %�������۵�
    LCC(j,1)=LCC1(j,1)+Cs(j,1);                                                             %Ŀ�꺯��1��ȫ��������ϵͳ���óɱ�
    CO2ot(j,1)=COpv*X(j,1)+COwt*X(j,2)+CObss*X(j,3)+COhp*X(j,4)+COge*Xgem+COtes*X(j,5);     %��ʼ̼�ŷ�
    CO2op(j,1)=(sum(elecnum)+sum(Phptes))*0.89*60;                                          %����̼�ŷ�
    CO2(j,1)=CO2ot(j,1)+CO2op(j,1);                                                         %Ŀ�꺯��2��ȫ��������̼�ŷ�
    EER(j,1)=(sum(Pwt)*X(j,2)+sum(Ppv)*X(j,1)-sum(More)+sum(Qhpscpl))/(sum(Pwt)*X(j,2)+sum(Ppv)*X(j,1)-sum(More)+sum(elecnum)+sum(Phptes)+sum(HeatingFuhe)+sum(CoolingFuhe));                    %Ŀ�꺯��3������
    %Qloss(j,1)=sum(Xloss)+sum(More)/(sum(Pwt)*X(j,2)+sum(Ppv)*X(j,1));                                                       %Լ��������Xloss(j,1)=0
    %Qloss(j,1)=sum(Xloss)+sum(More);
    Qloss(j,1)=sum(Qdeloss)+sum(More);
    %Qloss(j,1)=sum(More);
    
    %LPSP(j,1)=sum(LPS)/(sum(Fuhe)+sum(Php)*X(j,4));                                         %ȱ����
    %SPSP(j,1)=sum(More)/(sum(Fuhe)+sum(Php)*X(j,4));                                        %�˷ѵ���
    %SHSP(j,1)=(sum(MoreHeating)+sum(MoreCooling))/(sum(CoolingFuhe)+sum(HeatingFuhe));      %�˷���/����
    %SHSP(j,1)=sum(MoreHeating)/sum(HeatingFuhe)+sum(MoreCooling)/sum(CoolingFuhe);
    %LHSP(j,1)=(sum(LossHeating)+sum(LossCooling))/(sum(CoolingFuhe)+sum(HeatingFuhe));      %ȱ��/����
    %LHSP(j,1)=sum(LossHeating)/sum(HeatingFuhe)+sum(LossCooling)/sum(CoolingFuhe);
    %LGWMR(j,1)=sum(LossGeHeating)/(sum(Qhpscc)*X(j,4));                                     %����ȱʧ��
    %MGWMR(j,1)=sum(MoreGeHeating)/(sum(Qhpscc)*X(j,4));                                     %���ȹ�����
    %PopObj(j,1) = LCC(j,1)/100000000*2/3;
    PopObj(j,1) = LCC(j,1);
    PopObj(j,2) = CO2(j,1);
    PopObj(j,3) = 1-EER(j,1);
    PopObj(j,4) = Qloss(j,1);
end
    %%--�������ӽ⼯�������--%%
    

    
end


%F2=0.45*(LPSP+SPSP)+0.2*LHSP+0.2*SHSP+0.15*LSGSP;                                 %Ŀ�꺯��2��ϵͳ�ɿ���Ŀ��
%F2=0.35*LHSP+0.35*SHSP+0.3*LSGSP;
%F2=0.4*LPSP+0.3*LHSP+0.3*LGWMR;         


end
