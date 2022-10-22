%% ��������ӵ�����벢����ȫ������ֵ
%% �����д��л־Զ

%% ������޳�
function [Gbest_Index,Gbest_Value,pareto] = Update_Gbest(pareto,param)

[temp_num,~] = size(pareto.Value);                    %��֧��ֵ�����С����������
temp_erroe = temp_num-pareto.Capacity;                %����Ԥ���ڴ�ĳ̶�
pareto.Grid = zeros(temp_num,1);                      %����ӵ������վ���

%%ѡ��һ��Ŀ�꣬���Դ˽�������

[~,sequence] = sort(pareto.Value(:,1));               %Ĭ����������,sequence��index
pareto.Value = pareto.Value(sequence,:);
pareto.Index = pareto.Index(sequence,:);

%---------------------------����ӵ������---------------------------------

pareto.Grid(1) = inf;
pareto.Grid(end) = inf;

if temp_num > 2
    % ����ӵ������
    it_num = temp_num-1;
    for t = 2:it_num
        % ����� t ��ӵ������
        pareto.Grid(t) = sum(abs(pareto.Value(t-1,:)-pareto.Value(t+1,:)));
    end
end

temp_erroe = temp_num-pareto.Capacity;                                       %����Ԥ���ڴ�ĳ̶�

if temp_erroe > 0                                                            %����Ԥ������
    for u = 1:temp_erroe
        [~,temp_Index] = min(pareto.Grid);                                   %�޳�����̶ȴ�Ľ�,temp_Index��ӳ��Сֵ�ڸ��еڼ��С�1x1
        % ��ԭ�������޳�ӵ���̶�С�Ľ�
        pareto.Value(temp_Index,:) = [];                                     %��Сֵ������ֱ���޳�
        pareto.Index(temp_Index,:) = [];
        pareto.Grid(temp_Index) = [];
        
        if temp_Index ~= temp_num-1                                          %���㲹λ������ӵ������
            pareto.Grid(temp_Index) = ...
                sum(abs(pareto.Value(temp_Index-1,:)-pareto.Value(temp_Index+1,:)));
        end
        
        if temp_Index ~= 2                                                   %������һ��λ�õ�ӵ������
            pareto.Grid(temp_Index-1) = ...
                sum(abs(pareto.Value(temp_Index-2,:)-pareto.Value(temp_Index,:)));
        end
        
        temp_num = temp_num-1;
        
    end
end


%%  --------------ѡ��Gbest--------------------
if rand < 0.8                                                                   %���ӵ������
    
    if temp_num > 2
        %%����ӵ������
        iteration = temp_num-1;
        pareto.Grid(1) = inf;
        pareto.Grid(end) = inf;
        for t=2:iteration
            %%%�����t��ӵ������
            pareto.Grid(t) = sum(abs(pareto.Value(t-1,:)-pareto.Value(t+1,:)));
        end
        
        temp_file =2:iteration;
        [~,temp_index] = max(pareto.Grid(temp_file));                           %%�þ������Ľ�
        
        %%%�����ľ���Ľ���Ϊȫ�����Ž�
        Gbest_Index = repmat(pareto.Index(temp_file(temp_index),:),param.Size,1);
        Gbest_Value = repmat(pareto.Value(temp_file(temp_index),:),param.Size,1);
    else
        Gbest_Index = repmat(pareto.Index(randi(temp_num),:),param.Size,1);
        Gbest_Value = repmat(pareto.Value(randi(temp_num),:),param.Size,1);
    end
    
else
    
    if rand < 0.5                                                                %���
        index1 = randi(temp_num,param.Size,1);
        Gbest_Index = pareto.Index(index1,:);
        Gbest_Value = pareto.Value(index1,:);
    else                                                                         %�߽��
        temp_xp = randi([0 1],param.Size,1);
        index1 = temp_xp*(temp_num-1)+1;
        Gbest_Index = pareto.Index(index1,:);
        Gbest_Value = pareto.Value(index1,:);
    end
    
end



