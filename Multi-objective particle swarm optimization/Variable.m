%% 赋予计算变量初始值
%% 程序编写：谢志远

%% 定义算法参数及分配内容

param = Test_fun(param);                        %种群个体维度和取值范围，并将Particle=Param.
param = ParetoFront_fun(param);                 %种群真实前沿类型，并将Particle=Param.
param.Size =200;                                %种群规模,100行的来源

param.Imax = repmat(param.Imax,param.Size,1);   %种群规模，如果A是一个3x4x5的矩阵，有B = repmat(A,2,3)则最后的矩阵是6x12x5，此时得到的param.size x 原来列数的矩阵
param.Imin = repmat(param.Imin,param.Size,1);   %种群规模，同上

param.Maxgen =40000;                             %进化代数
param.C = [2,2];                                %个体最优/全局最优的加速学习因子
param.Vmax = 0.5*(param.Imax-param.Imin);       %个体速度上限，100x30--100x5
param.Vmin = -param.Vmax;                       %个体速度下限，100x30--100x5
param.Wmax = 0.9;                               %最大权重因子
param.Wmin = 0.5;                               %最小权重因子
param.W = 0.75;                                 %权重因子

%%定义空矩阵
pop.Index = zeros(param.Size,param.Dim);        %个体位置空矩阵，param.Size x param.Dim,400x5
pop.V = pop.Index;                              %个体速度空矩阵，400x5
pop.Pbest_Index = pop.Index;                    %个体最优值空矩阵位置信息，400x5
pop.Gbest_Index = zeros(1,param.Dim);           %全局最优值位置空矩阵位置信息，1x5
pop.Value = zeros(param.Size,param.Obj);        %个体解集空矩阵，400x3
pop.Pbest_Value = zeros(param.Size,param.Obj);  %个体最优值空矩阵解集信息，400x3
pop.Gbest_Value = zeros(param.Size ,param.Obj); %全局最优值空矩阵解集信息，400x3

%%%多目标参数
pareto.Capacity = 100;                          %档案最大容量
pareto.Index = [];                              %档案的位置信息集合
pareto.Value = [];                              %档案的适应度函数值集合



