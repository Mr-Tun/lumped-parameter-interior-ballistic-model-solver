function run_rk4(t0, y0, h, tf)
    local t = {}
    local y = {}
    local n = math.floor((tf - t0) / h) + 1
    for i = 1, n do
        t[i] = t0 + (i - 1) * h
    end
    y[1]=y0
    for i = 1, n - 1 do
       
        local k1 = lumped(y[i])
        local k2 = lumped({y[i][1] + h * k1[1] / 2, y[i][2] + h * k1[2] / 2, y[i][3] + h * k1[3] / 2, y[i][4] + h * k1[4] / 2, y[i][5] + h * k1[5] / 2, y[i][6] + h * k1[6] / 2})
        local k3 = lumped({y[i][1] + h * k2[1] / 2, y[i][2] + h * k2[2] / 2, y[i][3] + h * k2[3] / 2, y[i][4] + h * k2[4] / 2, y[i][5] + h * k2[5] / 2, y[i][6] + h * k2[6] / 2})
        local k4 = lumped({y[i][1] + h * k3[1], y[i][2] + h * k3[2], y[i][3] + h * k3[3], y[i][4] + h * k3[4], y[i][5] + h * k3[5], y[i][6] + h * k3[6]})
        y[i + 1] = {}
        for j = 1, 6 do
            y[i + 1][j] = y[i][j] + (h / 6) * (k1[j] + 2 * k2[j] + 2 * k3[j] + k4[j])
        end
    end

    return t, y
end

function lumped(y)
    local dydt = {}
    dydt[2] = (u1 / e1 * (y[5]^n)) * (y[2] < ZK and 1 or 0)
    dydt[1] = (chi + 2 * lambda * chi * y[2] + 3 * mu * chi * y[2]^2) * dydt[2] * (y[2] < 1 and 1 or 0) + (chi_s / ZK * (1 + 2 * lambda_s / ZK * y[2])) * dydt[2] * (y[2] >= 1 and y[2] < ZK and 1 or 0)
    dydt[3] = y[4]
    dydt[4] = (S * y[5]) / (phi * m)
    dydt[6] = l0 * (-delta * (alpha - 1 / rho_p) * dydt[1])
    dydt[5] = 1 / (S * (y[3] + y[6])) * (f * omega * dydt[1] - theta * phi * m * y[3] * dydt[4] - S * y[5] * (dydt[3] + dydt[6]))
    return dydt
end
-- defined the global variable
rho_p = 1600             --the density of the propellants 火药密度 kg/m^(-3)
omega = 1.16             --mass of charge 装药量   kg
f = 950000               --propellants force 火药力 (J/kg)
k = 1.25                 --gamma specific heat ratio 比热比
alpha = 0.001            --gas specific volume 分子余容  m^(-3)/kg
u1 = 5.127e-8            --burning velocity coefficient 烧燃速系数
n = 0.765                --pressure coefficient 装药压力指数 
e1 = 0.00055             --the thickness of propellants 火药弧厚 m
d = 0.00055              --the diameter of propellant grains holes 火药柱孔径 m
chi = 0.75               --the shape parameters of the propellant grains χ
lambda = 0.12            --the shape parameters of the propellant grains λ 
mu = 0                   --the shape parameters of the propellant grains μ
chi_s = 1.696            --the shape parameters of the propellant grains χs
lambda_s = -0.4104       --the shape parameters of the propellant grains λs
m = 2.8                  --piston/projectile mass 活塞/弹丸质量
diameter = 57e-3         --the caliber of the gun or the pump tube 炮管口径m
V0 = 0.00151             --the volume of the combustion chamber 药室体积m^3
lg = 3.624               --length of the launch tube 发射管长度
phi = 1.168              --drag coefficient of piston/projectile 次要功系数 φ1
p0 = 3e7                 --start pressure of projectile/piston 弹丸/活塞启动压力 Pa

-------
theta = k - 1                     -- θ=γ-1
delta = omega / V0                -- charge density
S = math.pi / 4 * diameter^2      -- across area of gun
l0 = V0 / S                       -- reduction length of combustion chamber 
rho = 0.2956 * (d / 2 + e1)       -- the shape parameters of the propellant grains ρ
ZK = (e1 + rho) / e1              -- 多孔火药碎粒全部燃完时的火药相对燃烧厚度

--compute the initial condition for the lumped parameter interior ballistic equations
l_0 = 0  -- displacement of projectile
v_0 = 0  -- velocity of projectile
psi0 = ((1 / delta) - (1 / rho_p)) / ((f / p0) + (alpha - (1 / rho_p)))  -- percentage of relative burned volume ψ
lpsi0 = l0 * (1 - delta / rho_p - delta * (alpha - 1 / rho_p) * psi0)    -- free volume reduction length lψ
Z0 = ((1 + 4 * lambda * psi0 / chi)^0.5 - 1) / (2 * lambda)              -- relative burned thickness Z0
t0 = 0   --initial time       
y0 = {psi0, Z0, l_0, v_0, p0, lpsi0} --initial condtions of solution vectors
h = 0.0001 -- time marching step dt
tf = 0.1  --end time&conditions 

test={}
test=lumped(y0)

t1, y1 = run_rk4(t0, y0, h, tf)

-- find the end conditions
ni = 1
while y1[ni][3] <= lg do
    ni = ni + 1
end
ni = ni - 1
lgi = ni

local current_dir = io.popen("cd"):read("*l")
print("当前工作目录是: " .. current_dir)
-- 打开一个文件以写入数据
local file = io.open("interior_ballistic_parameter_data.txt", "w")
-- 检查文件是否成功打开
if file then
    -- 遍历 y 表中的每个时间点
    for i = 1, ni do
        -- 将每个状态变量写入文件，用空格分隔
        for j = 1, #y1[i] do
            file:write(y1[i][j])
            if j < #y1[i] then
                file:write(" ")  -- 添加空格分隔符
            end
        end
        file:write(" ")  -- 添加空格分隔符
        file:write(t1[i])
        file:write("\n")  -- 换行
    end

    -- 关闭文件
    file:close()
else
    print("无法打开文件进行写入")
end
