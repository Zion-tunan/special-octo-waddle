%% ����Te��Tk����ֵ����������������
%% �����д��л־Զ

%% ��ͼ����
Tgeout=zeros(8760,1);
for i=1:8760
    Tgeout(i,1)=16.5+rand(1,1)-0.5;
    Te(i,1)=-7.818852817914556e-7*(i-4380)^2+5;
    Tk(i,1)=-4.170054836221096e-7*(i-4380)^2+46.5;
    x(i,1)=i;
    
        %��ů��part1
    for i = 1:2520
        Tk(i,1)=55;
        Te(i,1)=Tgeout(i,1)-10;
    end
    
    %���伾
    for i = 3625:5832
        Tk(i,1)=Tgeout(i,1)+10;
        Te(i,1)=0;
    end
    
    %��ů��part2
    for i = 7632:8760
        Tk(i,1)=55;
        Te(i,1)=Tgeout(i,1)-10;
    end
end


disp('��ʼ��ͼ...')
plot(x,Te(:,1),x,Tk(:,1));
title('�����¶��������¶�ȫ��ʵ��ֵ');
% axis([0 1.5 0 14]);
xlabel('\fontsize{12}Сʱ/h');
ylabel('\fontsize{12}�¶�/��');