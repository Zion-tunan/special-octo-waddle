%% �ж�֧���ϵ�����
%% �����д��л־Զ

%% temp_one = repmat(pop.Value(u,:),u-1,1)��temp_two = pop.Value(1:u-1,:)

function apflag = Dominate(temp_one,temp_two)

% temp_one = roundn(temp_one,-3);
% temp_two = roundn(temp_two,-3);

temp1 = all(temp_one <= temp_two,2);                         %all���������ж���ѡ��ĵ�Ԫ�ؼ��壻2Ϊ���з���ֵ.���ڼ���u�е�3�������ǲ���ȫ��1:u-1�е�����С������ȣ��ǵĻ�����1,1�����u������֧��1:u-1������.����u-1 x 1����
temp2 = all(temp_two <= temp_one,2);                         %���ڼ���u�е�3�������ǲ���ȫ��1:u-1�е����ݴ������ȣ��ǵĻ�����1.����u-1 x 1����
apflag = temp1.*temp2 + temp1-temp2;                         %1/-1/0.�����ʱ,����1.��u<1:u-1ʱ,temp1����1,temp2����0,apflag����1������֧��.��u>1:u-1ʱ,temp1����0,temp2����1,apflag����-1������֧��.����u-1 x 1����.0
%����0����֧��,����1�����u�����ݱ�ǰ�ߵ�����֧�������ȣ�����-1�����u������֧��ǰ�ߵ�����

end

