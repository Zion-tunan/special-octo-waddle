%% ��Ŀ���Ż�������
%% �����д��л־Զ

%% �������������

clear
close all
clc

global Ppv Pwt Qhpsh Qhpsc Php Tk Tgeout
[Ppv,Pwt,Qhpsh,Qhpsc,Php,Tk,Tgeout] = PVWTHP(1);                            %������̫�������ȱó���
CostFunction=@(x) ZDD(x);                                                   %�ɱ�����

turn_figure = 'on';

param.Test_function = 'UF10';                                               %���Ժ���
Variable;                                                                   %�������
N = 1;                                                                      %���Դ���
%performance = [];                                                           %���ܱ��֣�����Ŀ��ֵ���򻯣�

for t_num = 1:N
%% �򡢸����ʼ����Ⱥ�������ʼ��Ӧ�Ⱥ������״θ��·�֧�䵵��

    pop = Code_fun(pop,param);                                              %��ʼ����Ⱥ+������Ӧ�Ⱥ���
    
    pareto = Update_gather(pop,param,pareto);                               %��ʼ���·�֧��ֵ���
    
    [pop.Gbest_Index,pop.Gbest_Value,pareto] = Update_Gbest(pareto,param);  %��ʼ���������ѡȡȫ�����Ž�
    
%% �󡢵���Ѱ��

    for t = 1:param.Maxgen
        disp(['��������',num2str(t),' / ',num2str(param.Maxgen)])
        
        %%% �ٶȸ���
        param.W = param.Wmin+(param.Wmax-param.Wmin)...
            *(param.Maxgen-t)/param.Maxgen;                                 %����Ȩ��
        
        pop = Velocity(pop,param);                                          %�ٶȸ���
        
        pop = Update_Pbest(pop,param.Size);                                 %���¸��弫ֵ
        
        pareto = Update_gather(pop,param,pareto);                           %���·�֧��⵵��
        
        %b = ones(100,1);
        %a = b/pareto.Value(:,3);
        
        %%% ���������ȫ�����Ž��ѡȡ
        [pop.Gbest_Index,pop.Gbest_Value,pareto] = Update_Gbest(pareto,param);
        switch turn_figure
            case 'on'
                if param.Obj == 2
                    plot(pareto.Value(:,1),pareto.Value(:,2),'ro')
                    xlabel('{\itF}_1({\itx})','FontSize',13,'FontWeight','bold','FontName','Times New Roman')
                    ylabel('{\itF}_2({\itx})','FontSize',13,'FontWeight','bold','FontName','Times New Roman')
                    hold on
                    plot(param.paretoFornt(:,1),param.paretoFornt(:,2),'b-')
                    hold off
                    legend('��õ�ǰ��','��ʵ��ǰ��')
                elseif param.Obj == 4
                    
                    figure(1);
                    plot3(pareto.Value(:,1),pareto.Value(:,2),pareto.Value(:,3),'co')
                    xlabel('{\itF}_1({\itx})','FontSize',13,'FontWeight','bold','FontName','Times New Roman')
                    ylabel('{\itF}_2({\itx})','FontSize',13,'FontWeight','bold','FontName','Times New Roman')
                    zlabel('{\itF}_3({\itx})','FontSize',13,'FontWeight','bold','FontName','Times New Roman')
                    
                    figure(2);
                    plot(pareto.Value(:,1),pareto.Value(:,2),'ro')
                    xlabel('{\itF}_1({\itx})','FontSize',13,'FontWeight','bold','FontName','Times New Roman')
                    ylabel('{\itF}_2({\itx})','FontSize',13,'FontWeight','bold','FontName','Times New Roman')
                    
                    figure(3);
                    plot(pareto.Value(:,1),pareto.Value(:,3),'go')
                    xlabel('{\itF}_1({\itx})','FontSize',13,'FontWeight','bold','FontName','Times New Roman')
                    ylabel('{\itF}_3({\itx})','FontSize',13,'FontWeight','bold','FontName','Times New Roman')
                    
                    figure(4);
                    plot(pareto.Value(:,2),pareto.Value(:,3),'ko')
                    xlabel('{\itF}_2({\itx})','FontSize',13,'FontWeight','bold','FontName','Times New Roman')
                    ylabel('{\itF}_3({\itx})','FontSize',13,'FontWeight','bold','FontName','Times New Roman')
                    %����1��
                    %x=pareto.Value(:,1);
                    %y=pareto.Value(:,2);
                    %z=pareto.Value(:,3);
                    %[X,Y]=meshgrid(x,y);
                    %Z = zeros(size(X));
                    %for i=1:length(X)
                    %   for j = 1:length(Y)
                    %       Z(i,j) = z(x==x(i) && y==y(j));
                    %   end
                    %end
                    %surf(X,Y,Z);
                    
                    %����2
                    %x=pareto.Value(:,1);
                    %y=pareto.Value(:,2);
                    %z=pareto.Value(:,3);
                    %[X,Y,Z]=griddata(x,y,z,linspace(min(x),max(x))',linspace(min(y),max(y)),'v4');
                    %mesh(X,Y,Z);
                    %contourf(X,Y,Z);
                    %pcolor(X,Y,Z);
                    
                    %scatter3(pareto.Value(:,1),pareto.Value(:,2),pareto.Value(:,3),'bo','filled')

                    %hold on
                    %scatter3(param.paretoFornt(:,1),param.paretoFornt(:,2),param.paretoFornt(:,3),'bo','filled')
                    %alpha(0.15)
                    %hold off
                end
                title(num2str(param.Test_function))
                grid on
                pause(0.1)
        end
    end
    
end


