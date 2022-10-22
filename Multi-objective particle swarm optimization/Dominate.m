%% 判断支配关系并标记
%% 程序编写：谢志远

%% temp_one = repmat(pop.Value(u,:),u-1,1)，temp_two = pop.Value(1:u-1,:)

function apflag = Dominate(temp_one,temp_two)

% temp_one = roundn(temp_one,-3);
% temp_two = roundn(temp_two,-3);

temp1 = all(temp_one <= temp_two,2);                         %all函数用于判断所选择的的元素集体；2为按列返回值.用于检测第u行的3个数据是不是全比1:u-1行的数据小或者相等，是的话返回1,1代表第u行数据支配1:u-1的数据.生成u-1 x 1矩阵
temp2 = all(temp_two <= temp_one,2);                         %用于检测第u行的3个数据是不是全比1:u-1行的数据大或者相等，是的话返回1.生成u-1 x 1矩阵
apflag = temp1.*temp2 + temp1-temp2;                         %1/-1/0.当相等时,返回1.当u<1:u-1时,temp1返回1,temp2返回0,apflag返回1，代表支配.当u>1:u-1时,temp1返回0,temp2返回1,apflag返回-1，代表被支配.生成u-1 x 1矩阵.0
%返回0代表不支配,返回1代表第u行数据被前边的数据支配或者相等，返回-1代表第u行数据支配前边的数据

end

