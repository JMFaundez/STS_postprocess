close all
clear all
clc
addpath matlab_script/

% Load an unroll data
L = load('stsINT3.mat');
data = L.data_int;

x = data(:,:,1);
y = data(:,:,2);
U = data(:,:,3);
V = data(:,:,4);
W = data(:,:,5);
uu = data(:,:,7);
vv = data(:,:,8);
ww = data(:,:,9);
uv = data(:,:,11);
[ny,nx] = size(x);

alpha = zeros(ny,nx);
dx = x(2,:) - x(1,:);
dy = y(2,:) - y(1,:);
L = sqrt(dx.^2 + dy.^2);
a = acos(dy./L);

for i=1:nx
  alpha(:,i) = a(i);
end

%figure()
%plot(x(1,:),alpha(1,:))

%figure()
%plot(x(1,:),dy./L)

%figure()
%hold on
%plot(x(1,:),abs(dy),'k.')
%plot(x(1,:),dx)
%figure()
%mesh(x,y,x*0)
%legend('dy','dx')

ind = 500;
ar = alpha(1,ind:end);
xr = x(1,ind:end);
p = polyfit(xr,ar,10);
x1 = [0.05 0.1 0.18 0.25 0.28, 0.3];
alpha_rad = polyval(p,x1) +pi/2
alpha_deg = alpha_rad*180/pi

cos_a = cos(alpha);
sin_a = sin(alpha);

% Compute tangential components
tt = uu.*cos_a.^2 + vv.*sin_a.^2 + 2*uv.*cos_a.*sin_a;
T = U.*cos_a + V.*sin_a;
pp = tt -T.^2;
pp(abs(pp)<1e-14)=0;
p = sqrt(pp);

% Create normal and tangential coordinates
S = zeros(ny,nx);
N = zeros(ny,nx);
for n=1:ny
  S(n,:) = x(1,:);
  dx = x(n,:) - x(1,:);
  dy = y(n,:) - y(1,:);
  dh = sqrt(dx.^2+dy.^2);
  N(n,:) = dh;
end

% Compute cf
p_max = zeros(nx,1);
n_max = zeros(nx,1);
for i=1:nx
  [p_max(i), jmax] = max(p(:,i));
  n_max(i) = N(jmax,i);
end

Re = 5.33333e5;
cf1 = zeros(nx,1);
cf2 = zeros(nx,1);
fin = 15;
for i=1:nx
  cf1(i) = (-p(1,i) + p(fin,i))/(N(fin,i) - N(1,i))/Re;
  cf2(i) = (-T(1,i) + T(fin,i))/(N(fin,i) - N(1,i))/Re;
end

%cf1 = smoothdata(cf1,'Gaussian',10);

in0 = 500;
figure()
subplot(411)
contourf(S(:,in0:end),N(:,in0:end),p(:,in0:end), 'LineStyle','none')
xlim([0,0.4])
ylim([0,0.01])
ylabel('Normal direction')
title("$\sqrt{\overline{u_{t}'u_{t}'}}$",'Interpreter','latex')
colorbar('manual','Position',[0.93, 0.73, 0.02, 0.20])

subplot(412)
plot(S(1,in0:end),p_max(in0:end),'linewidth',1.5)
xlim([0,0.4])
ylabel("$\left(\sqrt{\overline{u_{t}'u_{t}'}}\right)_{max}$",...
	   'Interpreter','latex', 'FontSize',20)

stp = 1;
subplot(413)
plot(S(1,in0:stp:end),cf1(in0:stp:end),'linewidth',1.5)
xlim([0,0.4])
ylabel("$\frac{1}{Re}\left(\frac{\partial u'_tu'_t}{\partial n}\right)_{n=0}$",...
	   'Interpreter','latex', 'FontSize',20)

subplot(414)
plot(S(1,in0:end),cf2(in0:end),'linewidth',1.5)
xlim([0,0.4])
ylabel('$\frac{1}{Re}\left(\frac{\partial U_T}{\partial n}\right)_{n=0}$',...
	   'Interpreter','latex', 'FontSize',20)
xlabel('Chord')
