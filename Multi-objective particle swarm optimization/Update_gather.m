%% ���·�֧�䵵��
%% �����д��л־Զ

%% ȫ������Ϊ��֧���

function pareto = Update_gather(pop,param,pareto)
%% ��ɸѡ��ǰ��Ⱥ�ķ�֧���
gather_Logic = ones(param.Size,1);                           %֧���ϵ�жϣ�0Ϊ�����ܣ�֧��

for u = 2:param.Size
    temp = repmat(pop.Value(u,:),u-1,1);                     %pop.Value(u,:)�ĵ�U�е�ȫ���У�����ֵ�����u-1��x3��
    apflag = Dominate(temp,pop.Value(1:u-1,:));              %��pop.Value(u,:)�ĵ�U�е�ȫ���е�ֵ��1:u-1�е�ȫ����ֵ���бȽϲ��죬����1�����u������֧�����������ݣ�����-1��������������֧�䡣����u-1 x 1����
    gather_Logic(u) = -1*(sum(apflag==-1)>0)+1;              %����1:u-1�е������б���u������֧��ʱ��gather_Logic(u)=0.
    temp1 = 1:u-1;                                           %����[1 2 3 4 5 ... u-1]���󣬼�1 x u-1
    temp2 = apflag==1;                                       %������u������֧���λ�ñ���1������λ��0������һ��u-1x1��ֻ��0��1�ľ���
    gather_Logic(temp1(temp2)) = 0;                          %  ?ʲô��˼  0��֧���   Ӧ���ǵ�u�б�֧��ʱ��Ϊ0
    %if pop.Value(u,4)>200
    %    gather_Logic(u)=0;
    %end
end
%����֧���gather_Logic=1

gather_Logic = logical(gather_Logic);                        %�ѷ�0ֵ��Ϊ1��0����0������֧����Ϊ1��0��֧���.  ���logical�Ǳ���ӵģ������±�ʽ�Ӳ�����

%temp_Value = pop.Value(gather_Logic,:,m);                   %��ǰ��Ⱥ�з�֧���
%temp_Index = pop.Index(gather_Logic,:,m);                   %��ǰ��Ⱥ�з�֧����Ӧ�ľ��߱���

temp_Value = [pop.Value(gather_Logic,:);pareto.Value];       %��pop.Value�еķ�֧���ȥ��
temp_Index = [pop.Index(gather_Logic,:);pareto.Index];

%% ��ǰ��Ⱥ�з�֧��� �� ԭ�е����ļ���

[temp_Num,~] = size(temp_Value);                             %�����н������

gather_Logic = ones(temp_Num,1);                             %����������temp_Num x 1�ľ���

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

pareto.Value = temp_Value(gather_Logic,:);                    %��ǰ��Ⱥ�з�֧���
pareto.Index = temp_Index(gather_Logic,:);                    %��ǰ��Ⱥ�з�֧����Ӧ�ľ��߱���

