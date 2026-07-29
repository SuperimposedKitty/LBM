# 格子玻尔兹曼算法说明

本文档说明本项目中单相和两相格子玻尔兹曼方法的基本原理、代码实现方式和当前案例设置。

## 1. D2Q9 格子模型

项目采用二维九速度 D2Q9 模型。每个格点存储 9 个分布函数 `f_i`，方向包括 1 个静止方向、4 个轴向方向和 4 个对角方向。方向、权重和反向方向定义在 `include/lbm/lattice.hpp` 中。

宏观密度和速度由分布函数求矩得到：

```text
rho = sum_i f_i
rho * u = sum_i c_i * f_i
```

其中 `c_i` 是第 `i` 个离散速度。D2Q9 的格子声速平方为 `cs2 = 1/3`。

本次升级只改变碰撞与力项离散，不增加速度方向。模型仍是二维、等温、低马赫数方法；MRT 能改善稳定性和非流体模态误差，但不能把 D2Q9 变成高马赫数可压缩模型，也不能消除二值几何掩膜的阶梯边界误差。

## 2. MRT 与 BGK 碰撞模型

项目默认使用多松弛时间 MRT 碰撞，同时保留原单松弛时间 BGK 模型。两种模型共用 D2Q9 速度、权重、平衡分布和 9 个分布函数的内存布局。

BGK 模型为：

```text
f_i(x + c_i, t + 1) = f_i(x, t) - omega * (f_i - f_i_eq)
omega = 1 / tau
nu = cs2 * (tau - 0.5)
```

`tau` 控制运动黏度，必须大于 `0.5`，否则黏度非正，计算会失去物理意义。

MRT 先用固定矩阵 `M` 把分布函数转换到矩空间：

```text
m = M f
m* = m - S (m - m_eq)
f* = M^-1 m*
```

矩向量包含密度、能量、两个动量、两个能量通量和两个剪切应力矩。质量和单相动量模态不松弛；剪切应力模态仍使用 `1/tau`，因此运动黏度公式保持不变。能量、体积和能量通量模态使用 `MrtRelaxationRates` 中相互独立的松弛率，用于抑制不参与宏观流动的数值模态。

`SolverConfig::collision_model` 和 `TwoPhaseConfig::collision_model` 默认为 `CollisionModel::MRT`。设为 `CollisionModel::BGK` 可切回旧公式。MRT 的非守恒松弛率必须位于 `(0, 2)`。

每一步计算分为三段：

1. `compute_macroscopic()`：由分布函数求密度和速度。
2. `collide()`：按配置使用 MRT 或 BGK 把分布函数松弛到平衡态。
3. `stream_periodic()` 或 `stream_lid_driven_cavity()`：把碰撞后的分布函数沿格子方向迁移。

平衡分布 `f_i_eq` 使用低马赫数二阶展开：

```text
f_i_eq = w_i * rho * (1 + 3 c_i.u + 4.5 (c_i.u)^2 - 1.5 u.u)
```

该形式适合不可压或弱可压低速流动。实际使用时应保持速度远小于格子声速。

## 3. 单相案例

周期剪切波案例使用周期边界。初始速度场为正弦剪切波，随着黏性耗散逐渐衰减，可用于观察质量守恒和黏性扩散。

顶盖驱动方腔案例把四周设为固体壁面，顶壁以给定速度运动。固壁在迁移阶段使用反弹边界，顶盖额外加入运动壁面动量修正，从而形成经典方腔主涡结构。

两个单相案例最终只输出速度场 SVG 动图：

```text
result/periodic_shear_speed_animation.svg
result/lid_driven_cavity_speed_animation.svg
```

## 4. 两相 Shan-Chen 模型

两相求解器位于 `include/lbm/two_phase_solver.hpp` 和 `src/two_phase_solver.cpp`。当前实现使用两个分布函数集合：

```text
fa_i：红色注入相 A
fb_i：蓝色被驱替相 B
```

两个组分分别计算密度：

```text
rho_a = sum_i fa_i
rho_b = sum_i fb_i
rho = rho_a + rho_b
```

混合速度由总动量和半步力修正得到：

```text
u = (sum_i c_i * (fa_i + fb_i) + 0.5 * F_total) / rho
```

两相分离通过 Shan-Chen 伪势力实现。伪势函数为：

```text
psi(rho) = 1 - exp(-rho)
```

A 相受到邻近 B 相的伪势作用，B 相受到邻近 A 相的伪势作用。正的 `interaction_strength` 会促使两种组分分离并形成有限厚度界面。

带力碰撞使用 Guo 力项。BGK 模式用标量因子 `(1 - omega/2)` 修正；MRT 模式先把力源转换到矩空间，再按 `(I - S/2) M F` 修正每个模态，避免用同一个松弛频率处理所有物理矩。

两组分的 Shan-Chen 内力按无序格点链路成对计算。每次计算 A 格点受到的力时，同时给对应 B 格点施加等大反向力；反向组分组合也在同一链路内处理。因此计算域内部的流体相互作用全域合力在浮点误差范围内为零。开放入口和出口外侧使用零法向梯度的虚拟储液层补齐伪势邻居，避免边界格点因缺少半侧邻居产生各向异性加速；该部分属于外部储液层施加的边界作用，不计入内部力平衡。A/B 组分的动量矩按各自 `1/tau_a` 和 `1/tau_b` 向共同混合速度松弛，从而保留原多组分模型中的组分间动量交换。

## 5. 细管驱替和接触角

细管两相驱替案例位于 `examples/capillary_displacement.cpp`。当前设置为：

- 上下边界是固体细管壁面。
- 初始时管内流体全部为蓝色 B 相。
- 左边界持续注入红色 A 相。
- 右边界采用近似零梯度出口，复制近出口位置的密度和受限速度。
- 计算按设置的总步数执行，不再按出口红相占比提前结束。

接触角通过流体-固体黏附力近似体现。代码把 `contact_angle_degrees` 转换为 `cos(theta)`，再让固壁对 A/B 两相施加相反偏好：

```text
theta < 90 deg：A 相更润湿壁面
theta = 90 deg：中性润湿
theta > 90 deg：B 相更润湿壁面
```

这不是严格的几何接触角边界条件，而是 Shan-Chen 模型中常用的壁面润湿力近似。实际接触角还会受到网格分辨率、相互作用强度、密度比和壁面黏附强度共同影响。

## 6. 多孔介质和多尺度流动近似

当前细管案例在管道中间加入一段多孔介质。自由流动区域孔隙度为 `1.0`，多孔介质区域孔隙度为 `0.3`。多孔区不是固体墙，也不是把孔道逐个解析出来，而是在格点尺度上引入一个体积平均的阻力项，用来近似更小孔隙尺度对宏观流动的影响。

代码中的孔隙度场 `porosity_` 具有两层含义：

```text
epsilon = 1.0：自由流动区域
epsilon = 0.3：多孔介质代表体积单元
epsilon = 0.0：固体壁面
```

多孔介质等效渗透率采用 Kozeny-Carman 型关系：

```text
K = d_p^2 * epsilon^3 / (180 * (1 - epsilon)^2)
```

其中 `d_p` 是配置中的等效孔径，单位为格点长度。阻力项包含 Darcy 线性阻力和 Forchheimer 非线性惯性阻力：

```text
F_drag = -rho * (a_D * u_pore + a_F * |u_pore| * u_pore)
u_pore = u / epsilon
```

这里 `u` 是 LBM 求得的表观速度，`u_pore` 是孔隙内平均速度。阻力再按 A/B 两相密度占比分配到两个组分的 Guo 力项中。这样自由流区仍主要表现为两相 LBM，进入多孔区后则表现出 Darcy/Forchheimer 型宏观阻滞，是一种单网格上的多尺度近似。

## 7. 液滴撞击固体表面案例

液滴撞击案例位于 `examples/droplet_impact.cpp`，复用两相 Shan-Chen 求解器，但使用封闭固壁时间步 `step_closed()`。该时间步不施加入口/出口边界，分布函数撞到计算域边界和固体表面时执行反弹。

案例初始化方式如下：

- 红色 A 相被初始化为一个带有限厚度界面的圆形液滴。
- 蓝色 B 相填充液滴外部环境。
- 底部为固体撞击表面，侧边和顶部作为封闭反弹边界。
- 液滴初始给定向下速度；体力项可用于补充重力，本案例主要由初始下落速度驱动。
- 固体表面设置为较大接触角，使液滴在撞击后更容易收缩并出现回弹。

该案例输出：

```text
result/droplet_impact_animation.svg
```

控制台会输出液滴 A 相质心的初始高度、最低高度和最终高度。`final droplet center y - minimum droplet center y` 为正时，说明在计算后段出现了回弹趋势。

## 8. 两相可视化

两相动图输出的是相场：

```text
phi = (rho_a - rho_b) / rho
```

含义如下：

```text
phi ->  1：红色注入相占优
phi -> -1：蓝色被驱替相占优
phi ->  0：两相界面
```

之前界面看起来发白，是因为相场色带把 `phi = 0` 映射到接近白色的中间颜色。现在色带改为蓝-黄-红，界面附近用黄色显示，更容易从红蓝两相中辨认出来。

两相案例最终输出：

```text
result/capillary_phase_animation.svg
result/droplet_impact_animation.svg
result/geometry_displacement_animation.svg
```

细管中间的多孔介质区域会在动图上叠加斜线纹理，底色仍然保留相场颜色。

## 9. 几何掩膜导入与边界条件

通用几何读取器位于 `include/lbm/geometry_mask.hpp` 和 `src/geometry_mask.cpp`。它把 `.geom` 文本文件转换为与求解器网格一一对应的单元类型。文件第一行是物理上边界，而求解器内部 `y = 0` 位于底部，因此读取时会反转行坐标。

字符定义如下：

```text
.：自由流体单元
#：固体或障碍物单元
I：左边界入口单元
O：右边界出口单元
P：体积平均多孔流体单元
```

读取器会检查非空行等宽、字符合法性、入口/出口所在列以及是否存在流体单元。求解器初始化时还会检查至少存在一个入口和一个出口，并确认每个入口、出口与域内流体相连。错误信息包含行列位置，便于直接修改几何文件。

单相几何流使用 `initialize_masked_flow()` 和 `step_masked_flow()`。`#` 在迁移阶段执行半格反弹，`I` 每步重构为给定密度和速度的平衡分布，`O` 复制内侧相邻单元的分布函数，近似实现一阶零梯度出口。`P` 在单相求解器中与 `.` 相同。

两相几何驱替使用 `initialize_geometry_displacement()` 和普通 `step()`。入口持续注入红色 A 相，其他流体与出口初始填充蓝色 B 相；`P` 把 `porous_porosity` 写入孔隙度场，并复用第 6 节的 Darcy/Forchheimer 阻力。该模式的入口/出口只作用于文件中标记的格点，不改变细管驱替和液滴碰撞案例的默认边界。

默认几何 `geometry/channel_obstacle.geom` 同时包含菱形固体障碍和孔隙度 `0.3` 的多孔区域。对应输出为：

```text
result/geometry_flow_speed_animation.svg
result/geometry_displacement_animation.svg
```

当前边界是格点级阶梯近似，不包含 STL/CAD 曲面、插值反弹或子网格曲率修正。需要提高曲线精度时，应增加几何网格分辨率并重新标记单元。

## 10. 输出和构建约定

项目运行后只输出最终 SVG 动图到源码根目录下的 `result` 文件夹。编译产物按构建类型放置：

```text
build/RelWithDebInfo/bin/   exe, dll, pdb
build/RelWithDebInfo/lib/   lib, exp
```

默认构建类型是 `RelWithDebInfo`，推荐使用：

```powershell
.\scripts\build_vs2022.bat
```

## 11. 三维 D3Q19 与 D3Q27 扩展

三维模块与原有 D2Q9 模块并行存在，不改变 `D2Q9`、`Grid`、`Solver` 和
`TwoPhaseSolver` 的公开调用。新增接口为：

```text
Lattice3DModel::D3Q19 / D3Q27
Grid3D
Solver3D
TwoPhaseSolver3D
GeometryMask3D
```

D3Q19 包含静止方向、6 个面中心方向和 12 个棱方向，共 19 个分布函数，
默认用于三维单相流。D3Q27 再增加 8 个立方体角点方向，共 27 个分布函数，
对三维界面梯度和壁面方向具有更完整的立方对称性，因此两相求解器默认使用
D3Q27。两种格子的声速平方仍为 `1/3`，平衡分布仍使用低马赫数二阶展开。

三维 MRT 采用 Hermite 矩子空间分解。非平衡分布被分解为：

```text
密度模态
三个动量模态
体积应力模态
五个独立剪切应力模态
剩余高阶动理学模态
```

单相密度和动量保持守恒；剪切应力以 `1/tau` 松弛，体积应力和高阶模态分别
使用 `MrtRelaxationRates3D::bulk` 与 `kinetic`。两相 A/B 组分的动量以各自
`1/tau_a`、`1/tau_b` 向共同混合速度松弛。Guo 力项使用相同的矩分解，并对
每个子空间应用对应的 `(1 - s/2)` 修正。BGK 分支保留统一松弛公式。

三维两相求解器继续使用 Shan-Chen 多组分伪势。每条三维无序格点链路只计算
一次，并把作用力等大反向累加到两个端点。润湿力不再假定壁面位于底部，而是
根据当前流体格点周围的三维固体指示函数求局部壁面方向，因此体素障碍的侧面、
棱和曲面阶梯都能参与接触角作用。多孔单元继续使用 Kozeny-Carman 渗透率与
三维 Darcy/Forchheimer 阻力。

三维体素文件由 `GeometryMask3D` 读取，格式为：

```text
LBM_GEOMETRY_3D nx ny nz
SLICE 0
<ny 行等宽字符>
SLICE 1
<ny 行等宽字符>
...
```

字符仍为 `.`、`#`、`I`、`O`、`P`。切片按求解器 `z` 坐标递增排列，每个
切片的首行对应物理上边界。当前只支持 `x` 最小端入口和 `x` 最大端出口。
曲面仍是体素阶梯近似，不包含 STL/CAD 直接解析和插值反弹。

三维示例最终只生成一个动态 SVG，每帧同时显示三个正交中心切片：

```text
result/d3_lid_driven_cavity_speed_animation.svg
result/d3_geometry_flow_speed_animation.svg
result/d3_geometry_displacement_animation.svg
result/d3_droplet_impact_animation.svg
```

三维液滴示例用于验证球形界面、撞壁和铺展。当前小网格参数不把整体离壁回弹
作为默认验收条件；若研究三维回弹，需要进行网格收敛、表面张力标定和动态接触
角参数扫描，不能只提高人工壁面排斥力。

三维存储量显著增加。仅双缓冲分布函数，D3Q19 单相约为每格 `304 B`；
D3Q27 两相的四组分布函数约为每格 `864 B`，尚未计入密度、速度、力和固体场。
因此默认案例使用较小体素域。该扩展仍是等温、低马赫数、弱可压缩模型，不能
用于高马赫数可压缩流或直接替代三维热流模型。
