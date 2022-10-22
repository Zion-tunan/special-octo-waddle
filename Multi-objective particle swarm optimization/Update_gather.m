%% 更新非支配档案
%% 程序编写：谢志远

%% 全部更新为非支配解

function pareto = Update_gather(pop,param,pareto)
%% 先筛选当前种群的非支配解
gather_Logic = ones(param.Size,1);                           %支配关系判断（0为被（受）支配

for u = 2:param.Size
    temp = repmat(pop.Value(u,:),u-1,1);                     %pop.Value(u,:)的第U行的全三列，三个值。变成u-1行x3列
    apflag = Dominate(temp,pop.Value(1:u-1,:));              %将pop.Value(u,:)的第U行的全三列的值与1:u-1行的全三列值进行比较插旗，返回1代表第u行数据支配其他行数据，返回-1代表被其他行数据支配。生成u-1 x 1矩阵
    gather_Logic(u) = -1*(sum(apflag==-1)>0)+1;              %即当1:u-1行的数据有被第u行数据支配时，gather_Logic(u)=0.
    temp1 = 1:u-1;                                           %生成[1 2 3 4 5 ... u-1]矩阵，即1 x u-1
    temp2 = apflag==1;                                       %将被第u行数据支配的位置标上1，其他位置0，返回一个u-1x1的只有0和1的矩阵
    gather_Logic(temp1(temp2)) = 0;                          %  ?什么意思  0是支配解   应该是第u行被支配时设为0
    %if pop.Value(u,4)>200
    %    gather_Logic(u)=0;
    %end
end
%最后非支配解gather_Logic=1

gather_Logic = logical(gather_Logic);                        %把非0值变为1，0还是0，即非支配解变为1，0是支配解.  这个logical是必须加的，否则下边式子不成立

%temp_Value = pop.Value(gather_Logic,:,m);                   %当前种群中非支配解
%temp_Index = pop.Index(gather_Logic,:,m);                   %当前种群中非支配解对应的决策变量

temp_Value = [pop.Value(gather_Logic,:);pareto.Value];       %将pop.Value中的非支配解去除
temp_Index = [pop.Index(gather_Logic,:);pareto.Index];

%% 当前种群中非支配解 和 原有档案的集合

[temp_Num,~] = size(temp_Value);                             %档案中解的数量

gather_Logic = ones(temp_Num,1);                             %生成新数量temp_Num x 1的矩阵

for u = 2:temp_Num
    temp = repmat(temp_Value(u,:),u-1,1);
    apflag = Dominate(temp,temp_Value(1:u-1,:));
    gather_Logic(u) = -1*(sum(apflag==-1)>0)+1;
    temp1 = 1:u-1;
    temp2 = apflag==1;
    gather_Logic(temp1(temp2)) = 0;
    %if pop.Value(u,4)>200
    %    gather_Logic(u)=0;
    %end
end

gather_Logic = logical(gather_Logic);

pareto.Value = temp_Value(gather_Logic,:);                    %当前种群中非支配解
pareto.Index = temp_Index(gather_Logic,:);                    %当前种群中非支配解对应的决策变量

