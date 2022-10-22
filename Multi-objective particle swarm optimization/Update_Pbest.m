%% 产生个体最优值
%% 程序编写：谢志远

%% 通过支配关系替换
function pop = Update_Pbest(pop,param_Size)

doflag = Dominate(pop.Value,pop.Pbest_Value);                   % 判断支配关系
temp_Index = doflag==1;                                         % (A,B)--A(B)支配B(A)  doflag=1(-1)
pop.Pbest_Value(temp_Index,:) = pop.Value(temp_Index,:);        % 支配解替换个体最优
pop.Pbest_Index(temp_Index,:) = pop.Index(temp_Index,:);

temp_Index = doflag == 0;
usep = rand(param_Size,1)<0.5;
temp_Index = temp_Index&usep;

pop.Pbest_Value(temp_Index,:) = pop.Value(temp_Index,:);        % 支配解替换个体最优
pop.Pbest_Index(temp_Index,:) = pop.Index(temp_Index,:);
