%% ������������ʼֵ
%% �����д��л־Զ

%% �����㷨��������������

param = Test_fun(param);                        %��Ⱥ����ά�Ⱥ�ȡֵ��Χ������Particle=Param.
param = ParetoFront_fun(param);                 %��Ⱥ��ʵǰ�����ͣ�����Particle=Param.
param.Size =200;                                %��Ⱥ��ģ,100�е���Դ

param.Imax = repmat(param.Imax,param.Size,1);   %��Ⱥ��ģ�����A��һ��3x4x5�ľ�����B = repmat(A,2,3)�����ľ�����6x12x5����ʱ�õ���param.size x ԭ�������ľ���
param.Imin = repmat(param.Imin,param.Size,1);   %��Ⱥ��ģ��ͬ��

param.Maxgen =40000;                             %��������
param.C = [2,2];                                %��������/ȫ�����ŵļ���ѧϰ����
param.Vmax = 0.5*(param.Imax-param.Imin);       %�����ٶ����ޣ�100x30--100x5
param.Vmin = -param.Vmax;                       %�����ٶ����ޣ�100x30--100x5
param.Wmax = 0.9;                               %���Ȩ������
param.Wmin = 0.5;                               %��СȨ������
param.W = 0.75;                                 %Ȩ������

%%����վ���
pop.Index = zeros(param.Size,param.Dim);        %����λ�ÿվ���param.Size x param.Dim,400x5
pop.V = pop.Index;                              %�����ٶȿվ���400x5
pop.Pbest_Index = pop.Index;                    %��������ֵ�վ���λ����Ϣ��400x5
pop.Gbest_Index = zeros(1,param.Dim);           %ȫ������ֵλ�ÿվ���λ����Ϣ��1x5
pop.Value = zeros(param.Size,param.Obj);        %����⼯�վ���400x3
pop.Pbest_Value = zeros(param.Size,param.Obj);  %��������ֵ�վ���⼯��Ϣ��400x3
pop.Gbest_Value = zeros(param.Size ,param.Obj); %ȫ������ֵ�վ���⼯��Ϣ��400x3

%%%��Ŀ�����
pareto.Capacity = 100;                          %�����������
pareto.Index = [];                              %������λ����Ϣ����
pareto.Value = [];                              %��������Ӧ�Ⱥ���ֵ����



