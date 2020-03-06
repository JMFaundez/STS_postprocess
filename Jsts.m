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

% Choose the files that will be interpolated
files = [input5, input6, input7];
T = [t5; t6; t7];

% this is the output
output_sts = 'stsINT0.f00001';

% Number of elements x and y direction
nelx = 165;
nely = 40;

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
save('data_int.mat','data_int')
