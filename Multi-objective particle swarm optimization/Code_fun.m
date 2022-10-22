%% 生成种群初始位置、速度与计算初始适应度值
%% 程序编写：谢志远

%% rand与Fitness_fun

function pop = Code_fun(pop,particle)

%%初始化个体位置和位置信息
pop.Index = particle.Imin+(particle.Imax-particle.Imin).*...
    rand(particle.Size,particle.Dim);                                %个体位置,500*4
pop.V = particle.Vmin+(particle.Vmax-particle.Vmin)...
    .*rand(particle.Size,particle.Dim);                              %个体速度,500*4


%%%计算适应度函数值
    pop.Value= Fitness_fun(pop.Index,particle.Dim,particle.Test_function );            %pop.Index给出某个具体的个体位置，particle.Dim给出变量数，particle.Test_function转到UF10。对应PopObj = Fitness_fun(X,X_num,Test_function)

pop.Pbest_Index = pop.Index;                                                           %个体最优位置,400x5
pop.Pbest_Value = pop.Value;                                                           %个体最优适应度函数值,400x3