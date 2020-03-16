close all
clear all
clc
addpath matlab_script/

% Load an unroll data
L = load('stsINT2.mat');
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
dx = x(10,:) - x(1,:);
dy = y(10,:) - y(1,:);
L = sqrt(dx.^2 + dy.^2);
a = acos(dy./L);
s = zeros(nx,1);
s(2:end,1) = cumsum(sqrt((x(1,2:end)-x(1,1:end-1)).^2+(y(1,2:end)-y(1,1:end-1)).^2));
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

% Compute tangential components
tt = uu.*cos_a.^2 + vv.*sin_a.^2 + 2*uv.*cos_a.*sin_a;
T = U.*cos_a + V.*sin_a;
pp = tt -T.^2;
pp(abs(pp)<1e-14)=0;
p = sqrt(pp);

% Compute normal component
nn = uu.*sin_a.^2 + vv.*cos_a.^2 - 2*uv.*cos_a.*sin_a;
Nor = -U.*sin_a + V.*cos_a;
np2 = nn -Nor.^2;
np2(abs(np2)<1e-14)=0;
np = sqrt(np2);

% Compute spanwise component
zz = ww - W.^2;
zz(abs(zz)<1e-14)=0;
zp = sqrt(zz);

u_pert = sqrt((pp + np2 + zz)/3);
U_av = sqrt(T.^2 + Nor.^2 + W.^2);


integrand = u_pert;%./U_av;

dat1 = base_case('fringe_m30.f00008',165,40);
y99 = dat1.y99;
TI = zeros(nx,1);
for i=1:nx
  for j=1:ny
    if y99(i)<=N(j,i)
    	int_i = integrand(j:end-1,i);
    	n_i = N(j:end-1,i);
    	TI(i) = trapz(n_i,int_i)/trapz(n_i,ones(length(n_i),1));
    	break
    end
  end
end

figure()
plot(S(1,ind:end),TI(ind:end))
xlim([0, 0.35])
xlabel('Chord')
ylabel('Ti')

% Compute cf
p_max = zeros(nx,1);
n_max = zeros(nx,1);
for i=1:nx
  [p_max(i), jmax] = max(p(:,i));
  n_max(i) = N(jmax,i);
end

in0 = 490;
x(1,in0)
Re = 5.33333e5;
s = s - s(in0);
s = s*Re;
cf1 = zeros(nx,1);
cf2 = zeros(nx,1);
fin = 15;
for i=1:nx
  cf1(i) = (-p(1,i) + p(fin,i))/(N(fin,i) - N(1,i))/Re;
  cf2(i) = (-T(1,i) + T(fin,i))/(N(fin,i) - N(1,i))/Re;
end
figure()
contourf(S(:,in0:end),N(:,in0:end),u_pert(:,in0:end), 'LineStyle','none')
xlim([0,0.4])
ylim([0,0.01])
ylabel('Normal direction')
xlabel('Chord')
colorbar()
%cf1 = smoothdata(cf1,'Gaussian',10);

figure()
subplot(222)
contourf(S(:,in0:end),N(:,in0:end),p(:,in0:end), 'LineStyle','none')
xlim([0,0.4])
ylim([0,0.01])
ylabel('Normal direction')
xlabel('Chord')
title("$\sqrt{\overline{u_{t}'u_{t}'}}$",'Interpreter','latex')
colorbar('manual','Position',[0.93, 0.73, 0.02, 0.20])

subplot(221)
plot(s(in0:end),p_max(in0:end),'linewidth',1.5)
%xlim([0,0.4])
xlim([0, max(s)])
grid on
ylabel("$\left(\sqrt{\overline{u_{t}'u_{t}'}}\right)_{max}$",...
	   'Interpreter','latex', 'FontSize',20)

stp = 1;
subplot(223)
plot(s(in0:stp:end),cf1(in0:stp:end),'linewidth',1.5)
xlabel("$Re_s$",'Interpreter','latex','FontSize',14)
%xlim([0,0.4])
xlim([0, max(s)])
grid on
ylabel("$\frac{1}{Re}\left(\frac{\partial\sqrt{ u'_tu'_t}}{\partial n}\right)_{n=0}$",...
	   'Interpreter','latex', 'FontSize',20)

subplot(224)
plot(s(in0:end),cf2(in0:end),'linewidth',1.5)
xlim([0, max(s)])
grid on
ylabel('$\frac{1}{Re}\left(\frac{\partial U_T}{\partial n}\right)_{n=0}$',...
	   'Interpreter','latex', 'FontSize',20)
xlabel("$Re_s$",'Interpreter','latex','FontSize',16)

figure()
contourf(S(:,in0:end),N(:,in0:end),np(:,in0:end), 'LineStyle','none')
xlim([0,0.4])
ylim([0,0.01])
ylabel('Normal direction')
xlabel('Chord')
title("$\sqrt{\overline{u_{n}'u_{n}'}}$",'Interpreter','latex')
colorbar()

figure()
plot(S(1,in0:end), s(in0:end))
[var ind] = min(abs(s-30300))

x(1,ind)
