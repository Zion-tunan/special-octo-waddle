%% ������Ⱥ��ʼλ�á��ٶ�������ʼ��Ӧ��ֵ
%% �����д��л־Զ

%% rand��Fitness_fun

function pop = Code_fun(pop,particle)

%%��ʼ������λ�ú�λ����Ϣ
pop.Index = particle.Imin+(particle.Imax-particle.Imin).*...
    rand(particle.Size,particle.Dim);                                %����λ��,500*4
pop.V = particle.Vmin+(particle.Vmax-particle.Vmin)...
    .*rand(particle.Size,particle.Dim);                              %�����ٶ�,500*4


%%%������Ӧ�Ⱥ���ֵ
    pop.Value= Fitness_fun(pop.Index,particle.Dim,particle.Test_function );            %pop.Index����ĳ������ĸ���λ�ã�particle.Dim������������particle.Test_functionת��UF10����ӦPopObj = Fitness_fun(X,X_num,Test_function)

pop.Pbest_Index = pop.Index;                                                           %��������λ��,400x5
pop.Pbest_Value = pop.Value;                                                           %����������Ӧ�Ⱥ���ֵ,400x3