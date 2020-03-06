close all
clear all
clc
addpath matlab_script/

%%
% this is the file with the right gll ordering
input_mesh='./fringe_m20.f00008';

% these are the file with the sts and their interval
input1="./sts_files/stsFST0.f00001";
t1 = [1.1384427 1.3376211];
input2="./sts_files/stsFST0.f00002";
t2 = [ 1.3376211 1.5865941];
input3="./sts_files/stsFST0.f00003";
t3 = [1.5865941  1.8347543];
input4="./sts_files2/stsFST0.f00001";
t4 = [1.931857  2.338979];
input5="./sts_files2/stsFST0.f00002";
t5 = [2.338979  2.7461015];
input6="./sts_files2/stsFST0.f00003";
t6 = [2.7461015  3.1486270];
input7="./sts_files2/stsFST0.f00004";
t7 = [3.1486270  3.5463085];
input8="./sts_files3/stsFST0.f00001";
t8 = [4.8238731 5.0282659];
input9="./sts_files3/stsFST0.f00002";
t9 = [5.0282659 5.2150766];
input10="./sts_files3/stsFST0.f00003";
t10 = [5.2150766 5.4018872];
input11="./sts_files3/stsFST0.f00004";
t11 = [5.4018872 5.5882415];
input12="./sts_files3/stsFST0.f00005";
t12 = [5.5882415 5.7750082];
input13="./sts_files3/stsFST0.f00006";
t13 = [5.7750082 5.9564271];
input14="./sts_files3/stsFST0.f00007";
t14 = [5.9564271 6.1417076];
input15="./sts_files3/stsFST0.f00008";
t15 = [6.1417076 6.3269880];

% Number of elements x and y direction
nelx = 165;
nely = 40;

% Choose the set of files to average and interpolate
set = 3;
output_sts = ['stsINT',num2str(set),'.mat'];
if set==1
    files = [input5, input6, input7];
    T = [t5; t6; t7];
elseif set==2
    files = [input9, input10, input11, input12, input13, input14, input15];
    T = [t9; t10; t11; t12; t13; t14; t15];
elseif set==3
    files = [input5, input6, input7, input9, input10, input11, input12, input13, input14, input15];
    T = [t5; t6; t7; t9; t10; t11; t12; t13; t14; t15];
end


[data_mesh,lr1,elmap,time,istep,fields,emode,wdsz,etag,header,status] = readnek(input_mesh);
[xx,yy] = reshapenek(data_mesh,nelx,nely);
[nx,ny] = size(xx);
N = length(files);
data_int = zeros(nx,ny,46);
tT = 0;
tic
for i=1:N
  i
  input_sts = files(i);
  ti = T(i,:);
  [data_sts,lr1,elmap,time,istep,fields,emode,wdsz,etag,header,status] = readnek(input_sts);
  data_i = interpolate_data(xx,yy,data_sts,ti);
  tT = tT +  ti(2) - ti(1);
  data_int = data_int + data_i;
end
data_int = data_int/tT;
data_int(:,:,1) = xx;
data_int(:,:,2) = yy;
toc
save(output_sts,'data_int')
