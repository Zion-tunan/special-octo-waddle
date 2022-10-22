%% �ٶȸ���
%% �����д��л־Զ

%% ����λ�����ٶȣ��������Ӳ�������Ӧ��ֵ
function pop = Velocity(pop,particle)

%%�ٶȸ���
pop.V = particle.W*pop.V+particle.C(1)*rand(particle.Size,particle.Dim).*...
    (pop.Pbest_Index-pop.Index)+particle.C(2)*rand(particle.Size,particle.Dim)...
    .*(pop.Gbest_Index-pop.Index);
%% �����ٶ�
poc_Logical = pop.V>particle.Vmax;
pop.V(poc_Logical)=particle.Vmax(poc_Logical);
poc_Logical = pop.V<particle.Vmin;
pop.V(poc_Logical) = particle.Vmin(poc_Logical);
%%%����λ��
pop.Index = pop.V+pop.Index;
%%%����λ��
poc_Logicalupper = pop.Index>particle.Imax;
poc_Logicallower = pop.Index<particle.Imin;
if rand < 0.8
    pop.Index(poc_Logicallower) = particle.Imin(poc_Logicallower);
    pop.Index(poc_Logicalupper)= particle.Imax(poc_Logicalupper);
else
    X = particle.Imin+(particle.Imax-particle.Imin).*rand(particle.Size,particle.Dim);
    pop.Index(poc_Logicallower) = X(poc_Logicallower) ;
    pop.Index(poc_Logicalupper) = X(poc_Logicalupper) ;
end
%%%������Ӧ�Ⱥ���ֵ
    pop.Value = Fitness_fun(pop.Index,particle.Dim,particle.Test_function );

