# 格子玻尔兹曼算法说明

本文档说明本项目中单相和两相格子玻尔兹曼方法的基本原理、代码实现方式和当前案例设置。

## 文件驱动的初始化与求解

通用入口为 `lbm_blockMesh`、`lbm_setFields`、`lbm_solver`。案例文件格式、配置参数和操作顺序见 README 的“从案例文件驱动通用 LBM 求解器”。
`src/foam_case.cpp` 按词法 token 解析 ASCII 文件，保留字典与列表层次，不执行代码或预处理指令。
导入器从 faces/owner/neighbour 恢复每个单元的八个顶点及六个面，检查轴对齐、统一边长、面朝向、相邻单元中心间距和 boundary 连续覆盖。
根据几何中心建立 OpenFOAM cell label 到 LBM scalar index 的映射，非均匀字段按该映射赋值，不能依赖文件恰好采用 x 优先编号。

二维保留单层 z 方向空边界，三维使用 D3Q19/D3Q27。物理固壁外填充一层固体格点，包围盒中缺失体素也视为固体；固壁位于流体格点与固体格点之间的半格处。
可选 `0/solid` 在输入网格的单元中心用 0/1 标记流体/固体；`setFields` 可以通过球形或盒形区域赋值，导出时固体从流体网格剔除。默认仍为半格反弹；单相可选解析距离和 BFL 曲面插值，见第 12 节。导出网格仍保留体素阶梯面。
OpenFOAM 几何坐标通过快照的 origin/spacing 保留。标准 LBM 的活动方向间距必须相等；未提供 `constant/lbmProperties` 时，仅改变几何尺度不会自动完成物理标定。

四类求解器新增 `initialize_fields`，按用户初值建立 `f_i=f_i^eq(rho,U)` 或 `fa_i=f_i^eq(rhoA,U)`、`fb_i=f_i^eq(rhoB,U)`，保留既有构造与 step 接口。
封闭域与导入几何均复用掩膜迁移和两相作用力。开放边界限制在左入口/右出口；单相新增未知分布重构的速度入口/密度出口选项，两相保留旧边界。
初始 0/ 值按用户文件保留，边界条件从第一个时间步持续施加。这里不读取非平衡分布函数，不能当作数值状态完全一致的重启。

新增 `lbm_case_input` 测试覆盖二维/三维单相及两相四条路径、反转 cell 编号、uniform/nonuniform 初值、球形相区、盒形多孔区、p 到 rho 转换、输入不被求解器改写，以及非法网格/字段/边界拒绝。

## OpenFOAM 场导出

四类求解器的 `openfoam_snapshot()` 复制当前宏观场；共用 `OpenFoamWriter` 将其写为 ASCII OpenFOAM 数据。
不会修改分布函数或执行额外碰撞。二维网格在 z 方向挤出一层，三维输出全部流体体素。
每个流体格点对应一个六面体，边长取输入网格的 `spacing`，面位于格点中心两侧半格处。
固体不作为流体单元导出，内部障碍的暴露面组成 `walls`；多孔区域仍是流体单元，用 `porosity` 表示。
内部面排在边界面之前，owner 编号小于 neighbour，顶点顺序保证法向由 owner 指向 neighbour。
共用顶点去重，边界按 patch 连续存储。二维前后面为 `empty`，其余场边界用邻接单元值作 `calculated` 后处理值。

未启用 SI 时，场采用无量纲格子单位，时间目录为整数步数；SI 单相模式按第 12 节转换物理单位和时间。
`U` 是求解器宏观混合速度，`rho` 为总密度，格子模式 `p = cs2 * rho = rho/3` 为理想部分压强。
`speed=|U|` 用于直接着色，`vorticity=∇×U`、`divergence=∇·U` 用相邻流体格点相对输入几何坐标的中心差分计算；紧邻固体或外边界时采用单边差分，缺少两侧相邻流体格点时该方向导数取 0。SI 模式单位为 1/s；格子模式仍相对输入坐标求导。障碍边缘的涡量仅用于定性观察，不替代壁面剪切应力计算。
两相另存 `pBulk = cs2 * (rho + G * psi(rhoA) * psi(rhoB))`，其中 `psi(rho)=1-exp(-rho)`。
该式来自对称双组分相互作用的均匀体相状态关系；它不包含界面梯度压力张量、壁面润湿应力或 Darcy 阻力造成的附加项，不能当作完整界面法向应力。

`alpha.A = rhoA/(rhoA+rhoB)`，`alpha.B = 1-alpha.A`，`phase = alpha.A-alpha.B`。
这里的 alpha 是密度分数，可用于当前近等密度案例的饱和度可视化；未定义各相本征密度时，不宣称它是严格的体积饱和度。
多孔介质占整个体素的 A 相体积分数代理为 `porosity * alpha.A`，而孔隙内饱和度代理为 `alpha.A`。
输出不把 `rho_low` 微量组分截断，保留原始计算值。

文件是后处理交换格式，不是将 LBM 转换为 OpenFOAM 有限体积求解器：周期耦合、移动壁面、接触角和驱动力仍由原 LBM 算法负责。
仅提供用于识别案例的 `system/controlDict`，没有续算用的 `fvSchemes`、`fvSolution` 和材料模型。
时间序列在每个 `output_interval` 写盘，固定保存到 `result/openfoam/<案例名>`。重跑时按 `.lbm-times` 清单清理已知场和旧时间目录，再覆盖网格，避免积累旧帧；全场 ASCII 数据比 SVG 切片占用更多磁盘。
`0/` 保存真实初始场。ParaView 的 `Skip Zero Time` 会隐藏该目录，需取消勾选后刷新读取器时间列表，必要时重新打开案例。
`lbm_openfoam_tests` 检查网格闭合、正体积、面和 patch 索引、二维 empty 面、四类求解器场映射、饱和度和与非法拓扑修改。

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

两个单相案例仅输出完整 OpenFOAM 场时间序列：

```text
result/openfoam/periodic_shear/lbm.foam
result/openfoam/lid_driven_cavity/lbm.foam
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
result/openfoam/droplet_impact/lbm.foam
```

控制台会输出液滴 A 相质心的初始高度、最低高度和最终高度。`final droplet center y - minimum droplet center y` 为正时，说明在计算后段出现了回弹趋势。

## 8. 两相可视化

OpenFOAM 的 `phase` 字段保存相场：

```text
phi = (rho_a - rho_b) / rho
```

含义如下：

```text
phi ->  1：红色注入相占优
phi -> -1：蓝色被驱替相占优
phi ->  0：两相界面
```

界面的显示颜色由 ParaView 色带决定。例如选择蓝-黄-红色带，可将 `phi = 0` 的界面显示为黄色。程序不再生成带预设色带的 SVG。

两相案例最终输出：

```text
result/openfoam/capillary_displacement/lbm.foam
result/openfoam/droplet_impact/lbm.foam
result/openfoam/geometry_displacement/lbm.foam
```

细管中间的多孔介质区域用 `porosity` 字段标记，在 ParaView 中可单独着色显示。

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
result/openfoam/geometry_flow/lbm.foam
result/openfoam/geometry_displacement/lbm.foam
```

纯 `.geom` 导入仍是格点级阶梯近似。通用单相文件入口可额外提供解析曲面距离做插值反弹；未引入 STL/CAD 解析或任意非结构格子。

## 10. 输出和构建约定

项目仅按输出间隔保存 OpenFOAM 网格与场数据到源码根目录下的 `result` 文件夹。可执行程序及 DLL/PDB 放到根目录 `bin`，LIB/EXP 放到同级 `lib`，均不加构建版本子目录：

```text
bin/                      exe, dll, pdb
lib/   lib, exp
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
纯 `.geom3d` 曲面仍是体素阶梯近似。`initialize_fields` 的单相路径可携带解析距离并启用 BFL；两相不使用该插值边界。

三维示例仅保存完整三维 OpenFOAM 网格与场，使用 ParaView 查看时间序列：

```text
result/openfoam/d3_lid_driven_cavity/lbm.foam
result/openfoam/d3_geometry_flow/lbm.foam
result/openfoam/d3_geometry_displacement/lbm.foam
result/openfoam/d3_droplet_impact/lbm.foam
```

三维液滴示例用于验证球形界面、撞壁和铺展。当前小网格参数不把整体离壁回弹
作为默认验收条件；若研究三维回弹，需要进行网格收敛、表面张力标定和动态接触
角参数扫描，不能只提高人工壁面排斥力。

三维存储量显著增加。仅双缓冲分布函数，D3Q19 单相约为每格 `304 B`；
D3Q27 两相的四组分布函数约为每格 `864 B`，尚未计入密度、速度、力和固体场。
因此默认案例使用较小体素域。该扩展仍是等温、低马赫数、弱可压缩模型，不能
用于高马赫数可压缩流或直接替代三维热流模型。

## 12. 单相物理尺度与精度路径

### SI 与格子变量

`include/lbm/units.hpp` 集中保存长度步长 dx、时间步长 dt、参考密度 rho0 和参考压强 p0。
仅存在 `constant/lbmProperties` 的单相文件案例启用 SI；不改变原有库构造或时间步签名。
内部始终在 dxLB=dtLB=rho0LB=1 的均匀格子上计算：

```text
uSI = uLB * dx/dt
rhoSI = rho0 * rhoLB
nuLB = nuSI * dt/dx²
tau = 0.5 + 3*nuLB
pSI = p0 + rho0*(dx/dt)²*(rhoLB - 1)/3
forceDensityLB = forceDensitySI * dt²/(rho0*dx)
accelerationLB = accelerationSI * dt²/dx
```

压强输入按上式反解密度，不能指定任意独立的 p 和 rho；同时提供二者时检查其内部与固定边界的一致性。
SI 模式禁止显式 `tau`，避免与物理黏度矛盾。输出 `U/rho/p` 分别为 m/s、kg/m³、Pa，时间目录为秒。
旧格子案例的压强仍为 rho/3，不减参考密度。二维 z 厚度独立，不参与 dx 和 tau 的推导。

`lbm_refineCase` 以原始均匀初场为模板执行扩散尺度细化：dx/r、dt/r²，因此 uLB/r、nuLB 和 tau 不变，Re 与物理时长不变。
初场几何区域在重新 blockMesh/setFields 后按细网格中心重新取样；它不是对已有非均匀场的插值重启。
目前仅支持单个均匀 hex 的 SI 单相模板；不做局部 AMR 或并行划分。

### 单相 Guo 力项

二维与三维单相碰撞共享总力密度 `F=force_density+rho*acceleration`。宏观速度为 `(sum(c*f)+F/2)/rho`。
原始源项为：

```text
Gi = wi * [3*(ci-u)·F + 9*(ci·u)*(ci·F)]
BGK: f* = f - omega*(f-feq) + (1-omega/2)*G
MRT: m* = m - S*(m-meq) + (I-S/2)*M*G
```

单相动量矩守恒，松弛率为 0，不能为了加力把它设置成剪切松弛率。
初始化使用 `f=feq-G/2`，使用户指定的初始宏观速度在半步力修正后仍精确成立。
测试覆盖 BGK/MRT 的 D2Q9、D3Q19、D3Q27 均匀加速，并用受力 Poiseuille 流检验黏度和速度剖面。
本轮未增加单相多孔阻力，也未把两相参数自动转换为 SI。

### 开放边界与 BFL 曲面

默认 `OpenBoundaryModel::Legacy` 保持旧单相调用的行为。`VelocityPressure` 仅重构迁移后未知的 x 向分布，保留已知分布。
入口指定 u 求 rho，出口指定 rho 求法向 u；使用反向非平衡关系和切向动量修正，使密度与三个动量矩一致。
加力时重构目标是原始动量 `rho*u-F/2`，体积力和加速度均参与入口密度/出口速度的代数求解。
文件入口仍仅支持 x 最小端入口、x 最大端出口；`inletProfile parabolic` 仅二维且入口值代表截面平均速度。
该压力出口不是无反射边界，也没有实现对流出口；不保证截短计算域后结果不变。

`constant/lbmGeometry` 的解析固体并集提供有符号距离 d，流体侧正、固体侧负。
流体 x 沿离散方向 ci 命中固体时，壁面链路比例 `q=d(x)/(d(x)-d(x+ci))`。
定义碰撞后分布 fout=f*i(x)、fopp=f*opp(x)、fbehind=f*i(x-ci)，壁面修正 C=6*wi*rho*ci·uw：

```text
q < 0.5: fopp(x,t+1) = 2q*fout + (1-2q)*fbehind - C
q >= 0.5: fopp(x,t+1) = [fout+(2q-1)*fopp-C]/(2q)
```

无有效解析距离或缺少流体后邻点时回退半格反弹。`interpolationFallbacks` 报告本步缺后邻点的回退数；未提供距离的掩膜壁面使用默认半格法。
两端距离的线性插值是曲面交点近似，不是解析射线求交。解析几何支持球、无限长圆柱、盒体，不支持 STL、贴体网格或曲面移动。
两相边界和专用案例不在本轮改造范围；单相周期基准可通过库配置的 periodic_x/periodic_z 配合固壁运行。

### 诊断、性能与验证

通用求解器在监测步打印格子质量、密度极值、Ma、相邻监测时刻速度相对变化 dU、密度相对变化 dRho、速度散度 RMS、格子动能和入口/出口近似质量流量。
通用单相 `maxMach` 在初始与每个内部步检查，超限终止；`densityWarning` 在监测步超限只打印告警。默认仍按 endTime 计算，只有显式启用 convergence 才可能提前结束。
入口/出口通量是格子截面宏观矩近似，不包含开放边界存储变化，不能直接与单步质量变化等同。
`massBalance` 则使用每步未知分布重构引入的精确质量增量以及 BFL 引入的壁面质量修正，核对完整的离散质量收支。
`wallMass` 显式暴露 BFL 的非严格局部质量守恒，不能用很小的 massBalance 掩盖它。

壁面力按出射与反射分布的动量交换累加；标记为 obstacle 的解析/掩膜内部固体另计力系数。
`Cd_ref/Cl_ref` 使用 `0.5*rho_ref*U_ref²*A_ref` 归一化；二维 A_ref=referenceLength（单位深度），三维 A_ref=referenceLength²。
三维不是自动识别球投影面积，多个障碍也不逐个分解；不能不检查参考面积便把 Cd_ref 当作标准物体阻力系数。
computeSeconds/outputSeconds 分开计时；fluidMLUPS 只按流体单元更新计数，不代表已经做了高 Re 或湍流 DNS。

新增 `lbm_accuracy`（解析解和矩约束）、`lbm_reference`（Ghia/DFG 与三维圆管）、`lbm_precision_input`（SI、细化和输入输出一致性）测试。
详细实测误差、网格阶数、参考来源及工程适用限制见 [ACCURACY.md](ACCURACY.md)。

## 13. 稳态监测与一致条件网格研究

`include/lbm/convergence.hpp` 定义联合阈值与连续窗口计数，`run_foam_case` 共用于文件求解和圆柱验证案例。
`Solver`/`Solver3D` 构造、初始化、step/run 接口不变；逐步 Ma 护栏属于通用运行器，不会静默改变专用老案例的控制逻辑。
不启用 convergence 时仍保留总时间结束条件，且输出步与最后一步均保存。

监测间隔可独立于输出间隔设置；间隔为 0 时回退至输出间隔。速度和密度指标为：

```text
dU   = sqrt(sum_fluid |u(t)-u(t_previous_monitor)|² / max(sum_fluid |u(t)|²,1e-30))
dRho = sqrt(sum_fluid (rho(t)-rho(t_previous_monitor))² / max(sum_fluid rho(t)²,1e-30))
forceRange = max(max_window Cd-min_window Cd, max_window Cl-min_window Cl)
```

每个单相内部步记录力的上下界，下一窗口也包含上次末值，避免只比较两个端点漏掉窗口内的力振荡。
至少连续两个、默认五个完整窗口通过全部阈值，且各窗口起点不早于 minTime 和入口 ramp 结束时刻，才标记 converged。
窗口失败或出现非有限指标时连续计数清零；最后不足完整间隔的采样只报告数据，不累积通过次数。
无可用非零参考速度/长度或未标记障碍时不使用力系数判据；不能以此认证不存在参考面积的障碍阻力。

稳态停止会额外保存真实终止步，即使该步不在写盘间隔上；未达到阈值则输出 termination=endTime。
maxMach 在每步直接检查宏观速度数组，避免复制快照；computeSeconds 包含该护栏和逐步力采样，不包含完整采样诊断与写盘。
密度告警仍不作为硬终止；此配置只用于单相稳态问题，不适用于周期涡脱落或两相界面演化。

`lbm_cylinder_validation` 在三个网格上保持直径、密度、运动黏度、入口速度、Re、壁面模型和压力出口完全一致，
同步使用 dt 正比 dx²；另只加长出口距离作敏感性对照。各档采用相同物理监测时长和阈值，而非把“相同步数”误认为相同物理时长。
三档为非等比细化，只在阻力差同号且可解正阶数时求解实际格距比对应的观察阶数；不盲用 log2，也不保证已进入渐近区。
该研究针对本项目 Re=18 的通道工况，不宣称等价于 DFG Re=20 的自然出口基准。

## 14. 固定坐标压力探针与细网格研究

`sample_pressure_probes` 位于 `include/lbm/pressure_probes.hpp` 和 `src/pressure_probes.cpp`，
使用快照的 origin/spacing 将指定坐标变为格点坐标，不依赖 OpenFOAM cell 编号。
单相 `FoamCase::pressure_probes` 与文件 `pressureProbes ((x y z) ...);` 按列表顺序一一对应。
该列表在 SI 细化时原样保留，因此比较的是同一物理点，而非同一个数组下标。

每个活动方向取相邻中心点，局部坐标为 f；低端权重 1-f、高端权重 f，角点权重为各方向权重的乘积。
二维四点、三维八点；恰在节点上的零权重项不参与，二维 z 固定在单层中心。
插值 `rho-1` 再加回基准密度，最后用与写盘相同的 `UnitScale::pressure` 转换：

```text
p_lattice = rho / 3
p_SI      = referencePressure + densityScale * (dx/dt)^2 * (rho - 1) / 3
deltaP01  = p[0] - p[1]
```

正权重角点必须全部是有限正密度流体；否则报错，不删除固体点后重新归一化，不夹取最近流体点。
这只定义单元中心包围区域的插值，不能采样固壁或边界面外侧。暂不提供两相探针，避免漏掉相互作用压力项。
通用运行器在创建输出 writer 前校验初始探针，非法点不会清除上次结果；0、每个监测步和真实最终步均采样。
只有控制台 PROBES 行和 `FoamRunSummary::pressures`，不改变默认 OpenFOAM-only 文件输出。
探针只用于诊断，自动停机仍使用第 13 节的联合判据，不宣称压差已单独满足某个收敛阈值。

`--fine-study` 增加 D=12/16/24 和 D=16、L=9D 的对照，物理条件和监测阈值沿用原研究。
前后流体探针分别为 (0.004,0.006,0.0005)、(0.012,0.006,0.0005) m；
`deltaCp=deltaP/(0.5*rho0*U0^2)` 仅使用参考动压无量纲化，不是圆柱壁面压力系数分布。
分别对 Cd 和 deltaP 求网格差与出口敏感性；各自差不单调或阶数不可解时输出 unresolved。
精细网格不会自动消除模型、边界和时间误差，实测值及限制见 ACCURACY.md。
