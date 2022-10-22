%% 构建目标函数个数、变量维数与上下限
%% 程序编写：谢志远

%% 不同的情况列举

function particle = Test_fun(particle)

switch particle.Test_function
    case 'ZDT1'
        particle.Dim = 30;                                                %维度
        particle.Obj = 2;                                                 %目标个数（如果是约束问题则+1）
        particle.Imin = zeros(1,particle.Dim);                            %个体位置下限
        particle.Imax = ones(1,particle.Dim);                             %个体位置上限
    case 'ZDT2'
        particle.Dim = 30;                                                %维度
        particle.Obj = 2;                                                 %目标个数（如果是约束问题则+1）
        particle.Imin = zeros(1,particle.Dim);                            %个体位置下限
        particle.Imax = ones(1,particle.Dim);                             %个体位置上限
    case 'ZDT3'
        particle.Dim = 30;                                                %维度
        particle.Obj = 2;                                                 %目标个数（如果是约束问题则+1）
        particle.Imin = zeros(1,particle.Dim);                            %个体位置下限
        particle.Imax = ones(1,particle.Dim);                             %个体位置上限
    case 'ZDT4'
        particle.Dim = 10;                                                %维度
        particle.Obj = 2;                                                 %目标个数（如果是约束问题则+1）
        particle.Imin(1,1) = 0;                                           %个体位置下限
        particle.Imax(1,1) = 1;                                           %个体位置上限
        particle.Imin(1,2:particle.Dim) = -5*ones(1,particle.Dim-1);      %个体位置下限
        particle.Imax(1,2:particle.Dim) = 5*ones(1,particle.Dim-1);       %个体位置上限

    case 'ZDT6'
        particle.Dim = 10;                                                %维度
        particle.Obj = 2;                                                 %目标个数（如果是约束问题则+1）
        particle.Imin = zeros(1,particle.Dim);                            %个体位置下限
        particle.Imax = ones(1,particle.Dim);                             %个体位置上限
    case 'DTLZ1'
        particle.Obj = 3;    % 目标个数
        particle.Dim = 12;   % 变量个数
        particle.Imax = ones(1,particle.Dim);
        particle.Imin = zeros(1,particle.Dim); % 边界        
    case 'DTLZ2'
        particle.Obj = 3;    % 目标个数
        particle.Dim = 10;   % 变量个数
        particle.Imax = ones(1,particle.Dim);
        particle.Imin = zeros(1,particle.Dim); % 边界        
    case 'DTLZ3'
        particle.Obj = 3;    % 目标个数
        particle.Dim = 12;   % 变量个数
        particle.Imax = ones(1,particle.Dim);
        particle.Imin = zeros(1,particle.Dim); % 边界
        
    case 'DTLZ4'
        particle.Obj = 3;    % 目标个数
        particle.Dim = 12;   % 变量个数
        particle.Imax = ones(1,particle.Dim);
        particle.Imin = zeros(1,particle.Dim); % 边界
        
    case 'DTLZ5'
        particle.Obj = 3;  % 目标个数
        particle.Dim = 12;   % 变量个数
        particle.Imax = ones(1,particle.Dim);
        particle.Imin = zeros(1,particle.Dim); % 边界
        
    case 'DTLZ6'
        particle.Obj = 3;    % 目标个数
        particle.Dim = 10;   % 变量个数
        particle.Imax = ones(1,particle.Dim);
        particle.Imin = zeros(1,particle.Dim); % 边界
        
    case 'DTLZ7'
        particle.Obj = 3;    % 目标个数
        particle.Dim = 10;   % 变量个数
        particle.Imax = ones(1,particle.Dim);
        particle.Imin = zeros(1,particle.Dim); % 边界
        
    case 'UF1'
        particle.Obj = 2;    % 目标个数
        particle.Dim = 30;   % 变量个数
        particle.Imax = ones(1,particle.Dim);
        particle.Imin = [0,zeros(1,particle.Dim-1)-1]; % 边界
        
    case 'UF2'
        particle.Obj = 2;    % 目标个数
        particle.Dim = 30;   % 变量个数
        particle.Imax = ones(1,particle.Dim);
        particle.Imin = [0,zeros(1,particle.Dim-1)-1]; % 边界
        
    case 'UF3'
        particle.Obj = 2;    % 目标个数
        particle.Dim = 30;   % 变量个数
        particle.Imax = ones(1,particle.Dim);
        particle.Imin = zeros(1,particle.Dim);
        
    case 'UF4'
        particle.Obj = 2;    % 目标个数
        particle.Dim = 30;   % 变量个数
        particle.Imax = ones(1,particle.Dim);
        particle.Imin = zeros(1,particle.Dim);
        
    case 'UF5'
        particle.Obj = 2;    % 目标个数
        particle.Dim = 30;   % 变量个数
        particle.Imax = ones(1,particle.Dim);
        particle.Imin = zeros(1,particle.Dim);
        
    case 'UF6'
        particle.Obj = 2;  % 目标个数
        particle.Dim = 30;   % 变量个数
        particle.Imax = ones(1,particle.Dim);
        particle.Imin = zeros(1,particle.Dim);
        
    case 'UF7'
        particle.Obj = 2;    % 目标个数
        particle.Imax = ones(1,particle.Dim);
        particle.Imin = zeros(1,particle.Dim);
        
    case 'UF8'
        particle.Obj = 3;    % 目标个数
        particle.Dim = 30;   % 变量个数
        particle.Imax = ones(1,particle.Dim);
        particle.Imin = zeros(1,particle.Dim);
        
    case 'UF9'
        particle.Obj = 3;    % 目标个数
        particle.Dim = 30;   % 变量个数
        particle.Imax = ones(1,particle.Dim);
        particle.Imin = zeros(1,particle.Dim);
        
    case 'UF10'
        particle.Obj = 4;    % 目标个数 4  实际是3+1约束
        particle.Dim = 5;    % 变量个数
        particle.Imax = [7200 24 200 30 1000];  %太阳能板  风机  蓄电池  热泵 水箱
        particle.Imin = [1 1 1 5 1];
end