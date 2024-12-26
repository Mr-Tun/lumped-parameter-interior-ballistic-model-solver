clear;close all;%清屏并关闭窗口
%% 输入已知数据
global n rho_p omega f k alpha u1 e1 d chi lambda mu chi_s lambda_s m diameter V0 lg p0  phi theta delta S l0 rho ZK
rho_p=1600;           %火药密度 kg/m^(-3)
omega=1.16;           %弹丸装药量   kg
f=950000;             %火药力 (J/kg)
k=1.25  ;             %比热比κ
alpha=0.001;          %余容  m^(-3)/kg
u1=5.127*10^(-8);     %装药烧速系数
n=0.765;              %装药的压力指数n，本处为了让曲线更接近实际的曲线修改了压力指数
e1=0.00055;           %装药火药药弧厚2e1
d=0.00055;            %火药的直径 m
chi=0.75;             %火药形状特征量
lambda=0.12;          %火药形状特征量
mu=0;                 %装药形状特征量μ
chi_s=1.696;          %装药分裂点形状特征量
lambda_s=-0.4104;     %装药分裂点形状特征量'
m=2.8;                %弹丸质量  kg
diameter=57*10^(-3);  %火炮口径 m
V0=0.00151;           %药室容积  m^3
lg=3.624;             %火炮身管的行程l m
phi=1.168;            %次要功系数φ
p0=3*10^7;            %弹丸起动压力  pa
%% 常量计算
theta=k-1;                                                           %θ=κ-1
delta=omega/V0;                                                      %Δ装填密度 kg/m^(-3)
S=3.1415926/4*diameter^2;                                            %火炮炮膛横断面积m^2
l0=V0/S;                                                             %等效圆柱体的长度     
rho=0.2956*(d/2+e1);
ZK=(e1+rho)/e1;
%% 经典内弹道初值计算
l_0=0;                                                               %炮弹的行程初值 
v_0=0;                                                               %炮弹的速度初值
psi0=((1/delta)-(1/rho_p))/((f/p0)+(alpha-(1/rho_p)));               %已燃烧体积百分比ψ的初值计算
lpsi0=l0*(1-delta/rho_p-delta*(alpha-1/rho_p)*psi0);                                                    %自由容积缩径长lψ的初值计算
Z0=((1+4*lambda*psi0/chi)^0.5-1)/2/lambda;                           %相对已燃烧厚度Z0初值计算
%% 经典内弹道循环计算
t0 = 0;  % 初始时间
y0 = [psi0,Z0,l_0,v_0,p0,lpsi0];  % 初始状态变量
h = 0.0001;  % 步长
tf = 0.1;  % 结束时间
% 调用龙格-库塔差分程序
[t1, y1] = run_rk4(@f, t0, y0, h, tf);
y1=y1';
t1=t1';
%% 寻找启动点
ni=1;
while y1(ni,3)<=lg
      ni=ni+1;
end
ni=ni-1;
lgi=ni;

%% 第一次绘图
subplot(2,2,1)%利用subplot功能划分四块区域进行绘图
plot(y1(1:lgi,3),y1(1:lgi,5)./1000000,'b');
hold on
legend(  'closed', 'Location', 'northeast' );%添加图例
% title('预测的57mm高射炮经典闭膛理论p-l曲线');
xlabel('l/m');
ylabel('p/MPa');

subplot(2,2,2)
plot(t1(1:lgi,1).*1000,y1(1:lgi,5)./1000000,'b');
box on;
% title('预测的57mm高射炮经典闭膛理论p-t曲线');
legend(  'closed', 'Location', 'northeast' );%添加图例
xlabel('t/ms');
ylabel('p/MPa');

subplot(2,2,3)
plot(t1(1:lgi,1).*1000,y1(1:lgi,4),'b');
box on;
%title('预测的57mm高射炮经典闭膛理论v-t曲线');
legend( 'closed', 'Location', 'northeast' );%添加图例
xlabel('t/s');
ylabel('v/m*s^-1')

subplot(2,2,4)
plot(y1(1:lgi,3),y1(1:lgi,4),'b')
%title('预测的57mm高射炮经典闭膛理论v-l曲线');
legend( 'closed', 'Location', 'northeast' );%添加图例
xlabel('l/m');
ylabel('v/m*s^-1');

