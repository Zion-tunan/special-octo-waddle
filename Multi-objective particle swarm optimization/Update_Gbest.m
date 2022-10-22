%% 计算粒子拥挤距离并更新全局最优值
%% 程序编写：谢志远

%% 标记与剔除
function [Gbest_Index,Gbest_Value,pareto] = Update_Gbest(pareto,param)

[temp_num,~] = size(pareto.Value);                    %非支配街档案大小，返回行数
temp_erroe = temp_num-pareto.Capacity;                %超出预设内存的程度
pareto.Grid = zeros(temp_num,1);                      %定义拥挤距离空矩阵

%%选定一个目标，并以此进行排序

[~,sequence] = sort(pareto.Value(:,1));               %默认升序排序,sequence是index
pareto.Value = pareto.Value(sequence,:);
pareto.Index = pareto.Index(sequence,:);

%---------------------------计算拥挤距离---------------------------------

pareto.Grid(1) = inf;
pareto.Grid(end) = inf;

if temp_num > 2
    % 计算拥挤距离
    it_num = temp_num-1;
    for t = 2:it_num
        % 计算解 t 的拥挤距离
        pareto.Grid(t) = sum(abs(pareto.Value(t-1,:)-pareto.Value(t+1,:)));
    end
end

temp_erroe = temp_num-pareto.Capacity;                                       %超出预设内存的程度

if temp_erroe > 0                                                            %超出预设容量
    for u = 1:temp_erroe
        [~,temp_Index] = min(pareto.Grid);                                   %剔除距离程度大的解,temp_Index反映最小值在各列第几行。1x1
        % 在原档案中剔除拥挤程度小的解
        pareto.Value(temp_Index,:) = [];                                     %最小值所在行直接剔除
        pareto.Index(temp_Index,:) = [];
        pareto.Grid(temp_Index) = [];
        
        if temp_Index ~= temp_num-1                                          %计算补位上来的拥挤距离
            pareto.Grid(temp_Index) = ...
                sum(abs(pareto.Value(temp_Index-1,:)-pareto.Value(temp_Index+1,:)));
        end
        
        if temp_Index ~= 2                                                   %计算上一个位置的拥挤距离
            pareto.Grid(temp_Index-1) = ...
                sum(abs(pareto.Value(temp_Index-2,:)-pareto.Value(temp_Index,:)));
        end
        
        temp_num = temp_num-1;
        
    end
end


%%  --------------选择Gbest--------------------
if rand < 0.8                                                                   %最大拥挤距离
    
    if temp_num > 2
        %%计算拥挤距离
        iteration = temp_num-1;
        pareto.Grid(1) = inf;
        pareto.Grid(end) = inf;
        for t=2:iteration
            %%%计算解t的拥挤距离
            pareto.Grid(t) = sum(abs(pareto.Value(t-1,:)-pareto.Value(t+1,:)));
        end
        
        temp_file =2:iteration;
        [~,temp_index] = max(pareto.Grid(temp_file));                           %%用距离最大的解
        
        %%%把最大的距离的解作为全局最优解
        Gbest_Index = repmat(pareto.Index(temp_file(temp_index),:),param.Size,1);
        Gbest_Value = repmat(pareto.Value(temp_file(temp_index),:),param.Size,1);
    else
        Gbest_Index = repmat(pareto.Index(randi(temp_num),:),param.Size,1);
        Gbest_Value = repmat(pareto.Value(randi(temp_num),:),param.Size,1);
    end
    
else
    
    if rand < 0.5                                                                %随机
        index1 = randi(temp_num,param.Size,1);
        Gbest_Index = pareto.Index(index1,:);
        Gbest_Value = pareto.Value(index1,:);
    else                                                                         %边界解
        temp_xp = randi([0 1],param.Size,1);
        index1 = temp_xp*(temp_num-1)+1;
        Gbest_Index = pareto.Index(index1,:);
        Gbest_Value = pareto.Value(index1,:);
    end
    
end



