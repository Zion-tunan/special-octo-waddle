%% ����Ŀ�꺯������������ά����������
%% �����д��л־Զ

%% ��ͬ������о�

function particle = Test_fun(particle)

switch particle.Test_function
    case 'ZDT1'
        particle.Dim = 30;                                                %ά��
        particle.Obj = 2;                                                 %Ŀ������������Լ��������+1��
        particle.Imin = zeros(1,particle.Dim);                            %����λ������
        particle.Imax = ones(1,particle.Dim);                             %����λ������
    case 'ZDT2'
        particle.Dim = 30;                                                %ά��
        particle.Obj = 2;                                                 %Ŀ������������Լ��������+1��
        particle.Imin = zeros(1,particle.Dim);                            %����λ������
        particle.Imax = ones(1,particle.Dim);                             %����λ������
    case 'ZDT3'
        particle.Dim = 30;                                                %ά��
        particle.Obj = 2;                                                 %Ŀ������������Լ��������+1��
        particle.Imin = zeros(1,particle.Dim);                            %����λ������
        particle.Imax = ones(1,particle.Dim);                             %����λ������
    case 'ZDT4'
        particle.Dim = 10;                                                %ά��
        particle.Obj = 2;                                                 %Ŀ������������Լ��������+1��
        particle.Imin(1,1) = 0;                                           %����λ������
        particle.Imax(1,1) = 1;                                           %����λ������
        particle.Imin(1,2:particle.Dim) = -5*ones(1,particle.Dim-1);      %����λ������
        particle.Imax(1,2:particle.Dim) = 5*ones(1,particle.Dim-1);       %����λ������

    case 'ZDT6'
        particle.Dim = 10;                                                %ά��
        particle.Obj = 2;                                                 %Ŀ������������Լ��������+1��
        particle.Imin = zeros(1,particle.Dim);                            %����λ������
        particle.Imax = ones(1,particle.Dim);                             %����λ������
    case 'DTLZ1'
        particle.Obj = 3;    % Ŀ�����
        particle.Dim = 12;   % ��������
        particle.Imax = ones(1,particle.Dim);
        particle.Imin = zeros(1,particle.Dim); % �߽�        
    case 'DTLZ2'
        particle.Obj = 3;    % Ŀ�����
        particle.Dim = 10;   % ��������
        particle.Imax = ones(1,particle.Dim);
        particle.Imin = zeros(1,particle.Dim); % �߽�        
    case 'DTLZ3'
        particle.Obj = 3;    % Ŀ�����
        particle.Dim = 12;   % ��������
        particle.Imax = ones(1,particle.Dim);
        particle.Imin = zeros(1,particle.Dim); % �߽�
        
    case 'DTLZ4'
        particle.Obj = 3;    % Ŀ�����
        particle.Dim = 12;   % ��������
        particle.Imax = ones(1,particle.Dim);
        particle.Imin = zeros(1,particle.Dim); % �߽�
        
    case 'DTLZ5'
        particle.Obj = 3;  % Ŀ�����
        particle.Dim = 12;   % ��������
        particle.Imax = ones(1,particle.Dim);
        particle.Imin = zeros(1,particle.Dim); % �߽�
        
    case 'DTLZ6'
        particle.Obj = 3;    % Ŀ�����
        particle.Dim = 10;   % ��������
        particle.Imax = ones(1,particle.Dim);
        particle.Imin = zeros(1,particle.Dim); % �߽�
        
    case 'DTLZ7'
        particle.Obj = 3;    % Ŀ�����
        particle.Dim = 10;   % ��������
        particle.Imax = ones(1,particle.Dim);
        particle.Imin = zeros(1,particle.Dim); % �߽�
        
    case 'UF1'
        particle.Obj = 2;    % Ŀ�����
        particle.Dim = 30;   % ��������
        particle.Imax = ones(1,particle.Dim);
        particle.Imin = [0,zeros(1,particle.Dim-1)-1]; % �߽�
        
    case 'UF2'
        particle.Obj = 2;    % Ŀ�����
        particle.Dim = 30;   % ��������
        particle.Imax = ones(1,particle.Dim);
        particle.Imin = [0,zeros(1,particle.Dim-1)-1]; % �߽�
        
    case 'UF3'
        particle.Obj = 2;    % Ŀ�����
        particle.Dim = 30;   % ��������
        particle.Imax = ones(1,particle.Dim);
        particle.Imin = zeros(1,particle.Dim);
        
    case 'UF4'
        particle.Obj = 2;    % Ŀ�����
        particle.Dim = 30;   % ��������
        particle.Imax = ones(1,particle.Dim);
        particle.Imin = zeros(1,particle.Dim);
        
    case 'UF5'
        particle.Obj = 2;    % Ŀ�����
        particle.Dim = 30;   % ��������
        particle.Imax = ones(1,particle.Dim);
        particle.Imin = zeros(1,particle.Dim);
        
    case 'UF6'
        particle.Obj = 2;  % Ŀ�����
        particle.Dim = 30;   % ��������
        particle.Imax = ones(1,particle.Dim);
        particle.Imin = zeros(1,particle.Dim);
        
    case 'UF7'
        particle.Obj = 2;    % Ŀ�����
        particle.Imax = ones(1,particle.Dim);
        particle.Imin = zeros(1,particle.Dim);
        
    case 'UF8'
        particle.Obj = 3;    % Ŀ�����
        particle.Dim = 30;   % ��������
        particle.Imax = ones(1,particle.Dim);
        particle.Imin = zeros(1,particle.Dim);
        
    case 'UF9'
        particle.Obj = 3;    % Ŀ�����
        particle.Dim = 30;   % ��������
        particle.Imax = ones(1,particle.Dim);
        particle.Imin = zeros(1,particle.Dim);
        
    case 'UF10'
        particle.Obj = 4;    % Ŀ����� 4  ʵ����3+1Լ��
        particle.Dim = 5;    % ��������
        particle.Imax = [7200 24 200 30 1000];  %̫���ܰ�  ���  ����  �ȱ� ˮ��
        particle.Imin = [1 1 1 5 1];
end