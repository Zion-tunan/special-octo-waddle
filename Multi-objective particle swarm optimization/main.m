%% 三目标优化主程序
%% 程序编写：谢志远

%% Ⅰ、输入基本参数

clear
close all
clc

global Ppv Pwt Qhpsh Qhpsc Php Tk Tgeout
[Ppv,Pwt,Qhpsh,Qhpsc,Php,Tk,Tgeout] = PVWTHP(1);                            %风力、太阳能与热泵出力
CostFunction=@(x) ZDD(x);                                                   %成本函数

turn_figure = 'on';

param.Test_function = 'UF10';                                               %测试函数
Variable;                                                                   %定义变量
N = 1;                                                                      %测试次数
%performance = [];                                                           %性能表现（包括目标值正则化）

for t_num = 1:N
%% Ⅱ、赋予初始化种群、计算初始适应度函数并首次更新非支配档案

    pop = Code_fun(pop,param);                                              %初始化种群+计算适应度函数
    
    pareto = Update_gather(pop,param,pareto);                               %初始更新非支配街档案
    
    [pop.Gbest_Index,pop.Gbest_Value,pareto] = Update_Gbest(pareto,param);  %初始档案管理和选取全局最优解
    
%% Ⅲ、迭代寻优

    for t = 1:param.Maxgen
        disp(['迭代次数',num2str(t),' / ',num2str(param.Maxgen)])
        
        %%% 速度更新
        param.W = param.Wmin+(param.Wmax-param.Wmin)...
            *(param.Maxgen-t)/param.Maxgen;                                 %惯性权重
        
        pop = Velocity(pop,param);                                          %速度更新
        
        pop = Update_Pbest(pop,param.Size);                                 %更新个体极值
        
        pareto = Update_gather(pop,param,pareto);                           %更新非支配解档案
        
        %b = ones(100,1);
        %a = b/pareto.Value(:,3);
        
        %%% 档案管理和全局最优解的选取
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
                    legend('获得的前沿','真实的前沿')
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
                    %方法1：
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
                    
                    %方法2
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


