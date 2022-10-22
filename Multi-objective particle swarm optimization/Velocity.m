%% 速度更新
%% 程序编写：谢志远

%% 限制位置与速度，更新粒子并计算适应度值
function pop = Velocity(pop,particle)

%%速度更新
pop.V = particle.W*pop.V+particle.C(1)*rand(particle.Size,particle.Dim).*...
    (pop.Pbest_Index-pop.Index)+particle.C(2)*rand(particle.Size,particle.Dim)...
    .*(pop.Gbest_Index-pop.Index);
%% 限制速度
poc_Logical = pop.V>particle.Vmax;
pop.V(poc_Logical)=particle.Vmax(poc_Logical);
poc_Logical = pop.V<particle.Vmin;
pop.V(poc_Logical) = particle.Vmin(poc_Logical);
%%%更新位置
pop.Index = pop.V+pop.Index;
%%%限制位置
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
%%%计算适应度函数值
    pop.Value = Fitness_fun(pop.Index,particle.Dim,particle.Test_function );

