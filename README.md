# D2Q9 / D3Q19 / D3Q27 LBM C++ 示例项目

这是一个同时支持二维和三维的格子玻尔兹曼方法示例项目。原有 D2Q9
接口和案例保持不变，并新增 D3Q19 单相与 D3Q27 两相求解器。项目包含：

- 周期剪切波 `periodic shear wave`
- 顶盖驱动方腔流 `lid-driven cavity`
- 细管两相驱替 `capillary displacement`
- 液滴撞击固体表面并反弹 `droplet impact`
- 文件导入几何的单相障碍流 `imported geometry flow`
- 文件导入几何的两相驱替 `imported geometry displacement`
- D3Q19 三维顶盖驱动方腔流
- D3Q19 三维体素障碍流
- D3Q27 三维体素两相驱替
- D3Q27 三维液滴撞击与铺展

所有案例仅按输出间隔保存 OpenFOAM 网格与完整场时间序列，不再生成 SVG，统一写入项目根目录下的 `result` 文件夹。
`result` 与 `src` 同级；程序使用编译期写入的绝对路径，因此从 VSCode 运行或直接双击 exe 都会输出到同一个位置。

## 从案例文件驱动通用 LBM 求解器

现在可以修改案例目录中的网格、初始场和配置运行不同仿真，不必修改或重新编译案例 C++。
提供 Windows 原生命令，使用 OpenFOAM ASCII 字典和场文件的明确子集：

```powershell
.\scripts\build_vs2022.bat
.\bin\lbm_blockMesh.exe -case .\cases\channel2d
.\bin\lbm_setFields.exe -case .\cases\channel2d
.\bin\lbm_solver.exe -case .\cases\channel2d
```

`lbm_blockMesh` 只生成网格，不修改 `0/`；`lbm_setFields` 按区域修改已有初始场；`lbm_solver` 只读取输入并计算。
手工编辑 `0/` 后可直接运行求解器，**不要再执行 setFields，否则其字典会重新赋值**。
改变网格单元数后，应重建网格并把旧非均匀初始场恢复成 `uniform`，再执行 setFields；旧列表长度不匹配会报错。

```text
cases/channel2d/
  system/blockMeshDict    网格顶点、单元数、边界
  system/setFieldsDict    默认值和区域赋值
  system/lbmDict          格子、碰撞、松弛时间、两相及多孔参数
  system/controlDict     endTime、deltaT 和 writeInterval
  constant/lbmProperties 可选单相 SI 标定：密度、运动黏度、参考压强
  constant/lbmGeometry   可选解析固体，供单相曲面插值反弹使用
  constant/polyMesh/     points、faces、owner、neighbour、boundary
  0/U                   初始三分量速度及速度边界
  0/rho                 初始密度及密度边界，可用 0/p 替代
  0/alpha.A             两相时必须提供的组分密度分数
  0/porosity            可选材料孔隙度，默认 1
  0/solid               可选固体掩膜，0 为流体、1 为固体
```

结果写到 `result/openfoam/channel2d/lbm.foam`，不会覆盖输入案例的 `0/` 和网格。重复求解仍覆盖同名结果。
不同输入目录的末级名称应不同，否则它们对应相同结果目录。

提供五个模板：

| 案例目录 | 模型 | 可修改内容 |
| --- | --- | --- |
| `cases/channel2d` | D2Q9 单相圆柱绕流 | 192×96 网格，曲面插值反弹；格子单位，Re=18 |
| `cases/cylinder2d_si` | D2Q9 单相物理单位圆柱绕流 | 192×96 网格，0.024 m × 0.012 m，Re=18；SI 输入输出 |
| `cases/droplet2d` | D2Q9 封闭两相液滴 | sphereToCell 半径、相分布、湿润参数 |
| `cases/displacement3d` | D3Q27 三维驱替 | 入口相分数、boxToCell 多孔区 |
| `cases/closed3d` | D3Q19 三维封闭域 | 20³ 网格，正反向初始速度层；观察剪切流随时间衰减 |

运行 `lbm_blockMesh` 重建网格后，两个单相模板都**必须运行一次 `lbm_setFields`**：
`channel2d` 把圆柱写入 `0/solid`，`closed3d` 把初始速度层写入 `0/U`。
直接运行求解器只会读取当前 `0/`，不会自动执行区域赋值。
在 ParaView 中打开 `result/openfoam/channel2d/lbm.foam`，确认 Mesh Regions 中启用了 `internalMesh`，
点 `Apply` 和 `Reset Camera`，
跳到时间 6000 或 12000；用 `Surface` 表示并选择 `speed` 着色，固定色标范围约为 `0..0.035`，
即可看到圆柱两侧的加速区及后方低速尾流。切换到 `vorticity` 的 `Z` 分量，
使用发散色表并固定范围 `-0.2..0.2`，可区分圆柱上下符号相反的剪切层。
若想看流线，选中原始 `lbm.foam` 体网格，添加 `Stream Tracer`，矢量选 `U`；
若过滤器无法使用单元场，先添加 `Cell Data to Point Data` 再接 `Stream Tracer`。
将种子线放在圆柱上游、跨过通道高度。后期流场趋于稳定，时间动画变化较小；流线显示的是流动路径，
不是随时间移动的粒子。`closed3d` 是有限初始扰动的衰减实验，不是持续驱动流；用切片显示 `vorticity` 更合适。
`channel2d` 仍是格子单位案例，单独改变格点数不能视为网格收敛验证。SI 案例的等物理条件细化见下节；低雷诺数流动不保证出现交替脱落涡街。

### 单相精度升级与物理单位

新模板建议先复制到工作目录，避免 `setFields` 把版本控制中的均匀初值改成大列表：

```powershell
Copy-Item .\cases\cylinder2d_si .\build\cylinder2d_si -Recurse
.\bin\lbm_blockMesh.exe -case .\build\cylinder2d_si
.\bin\lbm_setFields.exe -case .\build\cylinder2d_si
.\bin\lbm_solver.exe -case .\build\cylinder2d_si
```

该模板运行 0.75 s（12000 格子步），入口在前 0.15 s 平滑启动。结果位于 `result/openfoam/cylinder2d_si/lbm.foam`；
用 `speed` 的 `0..0.07` m/s、`vorticity` 的 Z 分量约 `-400..400` 1/s 查看。它是启动瞬态，最后一帧不代表已经收敛。

存在 `constant/lbmProperties` 时启用 **仅单相** SI 模式，网格坐标为米、`deltaT` 为秒：

```text
referenceDensity 1000;
kinematicViscosity 6.666666666666667e-6;
referencePressure 0;
referenceLength 0.004;
```

此时 `0/U` 的 dimensions 为 `[0 1 -1 0 0 0 0]`，`0/rho` 为 `[1 -3 0 0 0 0 0]`，
`0/p` 为 `[1 -1 -2 0 0 0 0]`；`solid/porosity` 仍无量纲。`tau` 必须省略，由 `nu`、`deltaT` 和网格间距推导。
`p` 是相对 `referencePressure` 标定的等温弱可压压强，不是独立求解的不可压压力。
无 `lbmProperties` 的旧案例仍使用格子单位。SI 不适用于当前两相模型。

单相 `system/lbmDict` 可使用：

```text
openBoundary velocityPressure;
wallBoundary interpolated;
inletProfile uniform;
inletRampTime 0.15;
maxMach 0.1;
densityWarning 0.01;
```

`velocityPressure` 仅重构未知分布：入口 U 固定、rho/p 零梯度；出口 U 零梯度、rho/p 固定。
`legacy` 保留旧入口/出口算法。`inletProfile parabolic` 仅限二维 velocityPressure 模式，指定入口速度表示平均速度；入口平滑启动也仅用于单相。
单相 `bodyForce` 为单位体积力，`acceleration` 为加速度，两者相加；SI 模式分别用 N/m³、m/s²，否则为格子单位。
默认没有体力。`maxMach` 在初始及输出采样步检查，超限终止；密度波动超阈值仅告警，不会伪装成收敛。

`wallBoundary interpolated` 要求 `constant/lbmGeometry`，例如：

```text
solids
(
    cylinder { centre (0.008 0.006 0.0005); radius 0.002; axis (0 0 1); }
);
```

支持 `sphere`（centre/radius）、无限长 `cylinder`（centre/radius/axis）、`box`（min/max）的并集。
解析固体与 `0/solid` 合并；单相内部用有符号距离做 BFL 曲面插值反弹，额外掩膜固体没有距离时回退半格反弹。
**ParaView 的导出网格仍为体素阶梯面**，数值曲面插值不等于导出了贴体曲面网格。

保持物理几何、速度、黏度、总时间不变的细化命令（源模板的 `0/` 必须仍为 uniform）：

```powershell
.\bin\lbm_refineCase.exe -case .\cases\cylinder2d_si -output .\build\cylinder2d_si_r2 -factor 2
.\bin\lbm_blockMesh.exe -case .\build\cylinder2d_si_r2
.\bin\lbm_setFields.exe -case .\build\cylinder2d_si_r2
.\bin\lbm_solver.exe -case .\build\cylinder2d_si_r2
```

细化倍率为 2..8，活动方向格点数乘倍率，`deltaT` 除倍率平方，步数及输出间隔乘倍率平方，因此保持 `tau`、Re 和物理采样时间。
目标目录必须不存在；不自动覆盖输入，不重采样非均匀初场，不支持任意 polyMesh 或多块网格的细化。
二维加密两倍计算量约增至 16 倍，三维约 32 倍（未计 I/O）。解析解、参考数据和已验证误差见 [ACCURACY.md](ACCURACY.md)。

### 输入范围与边界

- 网格必须是轴对齐、活动方向等间距的正方形/立方体体素，可含挖空的固体区域；二维为一层 z 单元，前后面必须 `empty`，z 厚度可独立设置。
- 导入原生 `polyMesh` 会根据单元中心重新映射编号，支持外部 OpenFOAM `blockMesh` 生成的多块共形均匀网格；不支持非均匀加密、斜网格、曲面单元、内部零厚度挡板、周期 patch 或任意非结构网格。矩形包围盒含填充固体格点最多 200 万个。
- 自带 `lbm_blockMesh` 支持一个标准轴对齐 `hex`、`simpleGrading (1 1 1)`、空 `edges`/`mergePatchPairs`，以及命名为 `walls`、`inlet`、`outlet`、`frontAndBack` 的边界。它不是完整 OpenFOAM blockMesh 的替代品。
- 左侧入口：网格 `patch`，U 为 `fixedValue uniform (ux 0 0)`；默认 legacy 模式的 rho（或 p）和两相 alpha.A 也用统一 fixedValue，右侧出口的 U/rho/alpha.A 为 `zeroGradient`。单相 velocityPressure 模式改用上节的速度入口/压力出口。封闭域可以没有入口出口。外部网格的 patch 名称不限，类型和方向仍须符合上述条件。
- 静止固壁：网格 `wall`，U 为 `noSlip` 或零速度 `fixedValue`，标量为 `zeroGradient`。两相接触角通过 lbmDict 的 `contactAngle`/`wallAdhesion` 设置；尚不解析 OpenFOAM 接触角边界类型或移动壁面。原有顶盖、液滴碰撞专用 exe 保留原算法。
- 可选 `0/solid` 使用单元中心的 0/1 标量掩膜；`setFieldsDict` 可用 `boxToCell`、`sphereToCell` 赋值。标为 1 的格点从输出流体网格中去除，界面使用静止反弹边界；入口和出口不能标为固体。它不改变输入 `constant/polyMesh`，也不是亚网格曲面边界。
- `0/U`、`0/rho`（或 `0/p`）、`0/alpha.A` 支持 `uniform` 与 `nonuniform List<...>`，列表必须与输入网格单元数量一致；支持行注释、块注释、带引号的字典名。不支持 binary/gzip、宏引用、`#include`、`#codeStream` 或表达式，遇到会报错。
- 未启用 SI 时，字段必须声明无量纲 `dimensions [0 0 0 0 0 0 0]`；几何坐标保持输入数值，**仅缩放几何不会自动标定物理单位**。SI 单相字段维度按上节设置。
- 格子单位下 p 是理想部分压强，转换为 `rho = 3*p`；SI 换算见 LBM.md。rho 和 p 同时存在时内部与固定边界值都必须一致。两相由 rho 和 alpha.A 构造 rhoA/rhoB；入口 alpha.A 须严格位于 (0,1)，初始内部场允许 0 或 1。rhoA/rhoB、phase、pBulk 为派生输出，不是本入口的独立输入变量。
- `setFieldsDict` 支持 `defaultFieldValues`、`boxToCell`、`sphereToCell` 和 `volScalarFieldValue`/`volVectorFieldValue`，区域按单元中心选择，后面的区域覆盖前面的区域，保留字段的边界字典。

修改 `system/controlDict` 的 `endTime` 和 `writeInterval` 控制计算与输出；只支持 `startTime 0`、`writeControl timeStep`。格子模式要求 `deltaT 1`；SI 模式 `endTime` 必须为正 `deltaT` 的整数倍，`writeInterval` 始终为步数。
修改 `system/lbmDict` 可选择 singlePhase/twoPhase、D2Q9/D3Q19/D3Q27、BGK/MRT、tau/tauA/tauB。
两相另支持 interactionStrength、bodyForce、contactAngle、wallAdhesion、recoloring、poreDiameter、darcyDrag、forchheimerDrag；单相支持 Guo 体力/加速度，但暂不支持多孔阻力。
宏观初始场被重建为平衡分布（带力单相加半步源项修正），因此本入口是初始化求解，**不是保存全部分布函数的精确断点续算**。

安装了 OpenFOAM 时，也可以使用其 `blockMesh -case ...` 和 `setFields -case ...` 准备符合上述范围的 ASCII 文件，再运行 `lbm_solver`。
本项目不依赖 OpenFOAM 安装；所实现语法参考 [blockMesh 文档](https://doc.cfd.direct/openfoam/user-guide-v14/blockmesh) 和 [场文件说明](https://www.openfoam.com/documentation/user-guide/2-openfoam-cases/2-2-basic-inputoutput-file-format)。

## OpenFOAM 结果查看

十个案例均默认启用。运行方式和 `output_interval` 不变，初始步、采样步和最后一步均保存。
每次运行生成 `result/openfoam/<案例源文件名>/lbm.foam`，控制台会打印完整路径。
从 VSCode 或双击 exe 启动均使用同一项目结果根目录。同一案例重跑会覆盖网格并清理上次登记的时间步，因此缩短计算时间也不会混入旧帧。
不要同时运行同一个案例的多个进程。旧版本的运行编号子目录保留，但新结果不再写入其中。

```text
result/openfoam/geometry_displacement/
  lbm.foam
  constant/polyMesh/{points,faces,owner,neighbour,boundary}
  system/controlDict
  0/{U,speed,vorticity,p,rho,porosity,rhoA,rhoB,alpha.A,alpha.B,phase,pBulk}
  300/...
  600/...
```

在 ParaView 中打开 `lbm.foam`，点击 Apply，选择场并播放时间序列。三维输出包含全部体素，可在 ParaView 中查看切片或等值面。

**查看 0 时刻：**选中原始 `lbm.foam` 读取器，在 Properties 搜索 `Skip Zero Time`（必要时开启齿轮高级选项），取消勾选并 Apply。
该选项会忽略整个 `0/` 目录，并非求解器没有保存初始场。如果时间仍从第一个采样步开始，点击读取器的 Refresh/Refresh Times；仍不更新时删除该数据源并重新打开，首次 Apply 前取消勾选该选项。
更新后跳到第一帧，时间应为 0。输出文件不能控制 ParaView 客户端的 `Skip Zero Time` 偏好设置。
同一案例重新计算后，也需要刷新或重新打开读取器，避免显示旧网格或缓存场。

| 文件 | 含义 |
| --- | --- |
| `U` | 三分量混合物速度，二维第三分量为零；速度大小可在 ParaView 选择 Magnitude |
| `speed` | `U` 的模长，可直接用于着色 |
| `vorticity` | 相对输入几何坐标的离散旋度；紧邻固体和外边界采用单边差分，其余采用中心差分 |
| `divergence` | 相同差分规则计算的速度散度，SI 模式单位为 1/s |
| `rho`、`p` | 总密度和等温压强；格子模式 `p = rho / 3`，SI 模式按参考压强换算 |
| `porosity` | 求解器实际使用的孔隙度；单相 `P` 按普通流体处理，因此为 1 |
| `rhoA`、`rhoB` | 两相组分密度 |
| `alpha.A`、`alpha.B` | 密度分数 `rhoA/(rhoA+rhoB)` 和其补数，作为本模型的饱和度代理量 |
| `phase` | 两相序参量 `(rhoA-rhoB)/(rhoA+rhoB)` |
| `pBulk` | 两相 Shan–Chen 均匀体相压强近似，详见 LBM.md |

未启用 SI 的旧案例仍输出无量纲场，时间目录名为计算步数；SI 单相案例输出米、秒、m/s、kg/m³、Pa，涡量/散度为 1/s，时间目录名为物理秒数。
不同本征密度流体的真实体积饱和度不能直接等同于密度分数；纯相中保留的微量另一组分也不会被强行截成 0 或 1。

固体格点从流体网格中去除，流固界面形成 `walls`。几何入口/出口分别形成 `inlet`/`outlet`，其他外表面为 `outer`；二维前后面为 `empty`。
这些文件用于后处理，非空场边界采用 `calculated` 并保存邻接单元值，**不是 OpenFOAM 续算案例**，不将周期条件、移动壁面速度、润湿力等转换为有限体积求解边界条件。
案例不再生成 SVG，也不缓存用于 SVG 渲染的动画帧。显式 CSV/VTK 和 SVG 库接口仍保留供自定义调用。场文件在计算过程中写入，中途停止时已写完的时间步仍可读取；历史 SVG 文件不会被自动删除。

库调用：初始化后构造一次 writer，网格只写一次，后续快照只写场。改变网格或重跑时重新构造 writer。
导出器用 `.lbm-times` 登记自己生成的时间步，重跑仅逐文件清理这些时间目录中的已知场文件；不会递归删除案例目录。
案例根目录的用户文件保留。若旧时间目录含未知文件，程序报错而不删除；没有导出器清单的既有网格或时间目录也不会直接覆盖。

```cpp
lbm::OpenFoamWriter output(lbm::openfoam_result_path("my_case"), solver.openfoam_snapshot());
output.write(0, solver.openfoam_snapshot());
solver.step();
output.write(1, solver.openfoam_snapshot());
```

实现遵循 OpenFOAM 的 [polyMesh 说明](https://www.openfoam.com/documentation/user-guide/4-mesh-generation-and-conversion/4.1-mesh-description)
和 [场文件格式](https://www.openfoam.com/documentation/user-guide/2-openfoam-cases/2-2-basic-inputoutput-file-format)。

可选第三方兼容性检查：在安装了 `vtk` 的 Python 环境运行
`python scripts/validate_openfoam.py result/openfoam`。该脚本实际读取各个时间步，检查网格、场数组和两相分数之和；
正常编译与运行不需要 Python、VTK 或 OpenFOAM。C++ 自动测试另行验证正体积、面朝向、连通性及原始场值映射。

## 构建与测试

默认编译版本为 `RelWithDebInfo`。

```powershell
.\scripts\build_vs2022.bat
```

也可以手动运行：

```powershell
cmake -S . -B build -DCMAKE_BUILD_TYPE=RelWithDebInfo
cmake --build build --config RelWithDebInfo
```

编译产物目录统一为：

```text
bin/                      exe, dll, pdb
lib/   lib, exp
```

所有编译版本共用项目根目录的 `bin/` 和 `lib/`，不会创建版本子目录。切换编译版本时建议使用 `cmake --build build --config Debug --clean-first` 等方式完整重建，避免复用另一个版本的同名程序。

默认碰撞模型为 MRT。原有 D2Q9 构造函数和时间步接口不变；三维通过
`Solver3D`、`TwoPhaseSolver3D` 和独立配置新增。需要切回 BGK 时，可在相应配置中设置：

```cpp
lbm::SolverConfig config;
config.collision_model = lbm::CollisionModel::BGK;

lbm::TwoPhaseConfig two_phase_config;
two_phase_config.collision_model = lbm::CollisionModel::BGK;

lbm::Solver3DConfig config_3d;
config_3d.collision_model = lbm::CollisionModel::BGK;
```

数值测试不生成仿真结果；OpenFOAM 格式测试在 `build/openfoam_test_output` 下保留小网格验证文件：

```powershell
.\scripts\test_vs2022.bat
```

## 运行

```powershell
.\bin\lbm_periodic.exe
.\bin\lbm_cavity.exe
.\bin\lbm_capillary.exe
.\bin\lbm_droplet.exe
.\bin\lbm_geometry_flow.exe
.\bin\lbm_geometry_displacement.exe
.\bin\lbm_d3_cavity.exe
.\bin\lbm_d3_geometry_flow.exe
.\bin\lbm_d3_geometry_displacement.exe
.\bin\lbm_d3_droplet.exe
```

输出文件：

```text
result/openfoam/periodic_shear/lbm.foam
result/openfoam/lid_driven_cavity/lbm.foam
result/openfoam/capillary_displacement/lbm.foam
result/openfoam/droplet_impact/lbm.foam
result/openfoam/geometry_flow/lbm.foam
result/openfoam/geometry_displacement/lbm.foam
result/openfoam/d3_lid_driven_cavity/lbm.foam
result/openfoam/d3_geometry_flow/lbm.foam
result/openfoam/d3_geometry_displacement/lbm.foam
result/openfoam/d3_droplet_impact/lbm.foam
```

ParaView 显示建议（颜色由所选色带决定）：

- 单相案例：选择 `U` 的 Magnitude 查看流速。
- 细管两相驱替：选择 `alpha.A`，范围固定为 0 到 1。
- 液滴撞击：选择 `alpha.A` 或 `phase` 观察界面。
- 固体障碍：显示 `walls` 边界。
- 多孔介质：用 `porosity` 查看区域分布。
- 三维案例：在 ParaView 中添加 Slice 查看 XY、XZ 或 YZ 截面。

## 几何文件

几何案例默认读取 `geometry/channel_obstacle.geom`。也可以在命令行传入其他文件：

```powershell
.\bin\lbm_geometry_flow.exe .\geometry\channel_obstacle.geom
.\bin\lbm_geometry_displacement.exe .\geometry\channel_obstacle.geom
```

`.geom` 是等宽字符网格。文件第一行对应物理上边界，最后一行对应物理下边界；空行会被忽略。

```text
.  自由流体
#  固体或障碍物，使用静止反弹边界
I  左边界固定速度入口
O  右边界零梯度出口
P  多孔流体区；单相按普通流体，两相使用 Darcy/Forchheimer 阻力
```

每行字符数必须相同，`I` 只能出现在最左列，`O` 只能出现在最右列。障碍物曲线以格点阶梯边界近似，网格越密，形状分辨率越高。

三维几何使用 `.geom3d`。首行声明尺寸，随后按 `z` 从小到大列出切片；
每个切片中的第一行仍对应物理上边界：

```text
LBM_GEOMETRY_3D 24 12 8
SLICE 0
########################
...
SLICE 1
########################
I.............PPPPP....O
...
```

字符含义与 `.geom` 相同，入口仍限制在 `x = 0`，出口限制在
`x = nx - 1`。三维默认文件为 `geometry/channel_obstacle.geom3d`：

```powershell
.\bin\lbm_d3_geometry_flow.exe .\geometry\channel_obstacle.geom3d
.\bin\lbm_d3_geometry_displacement.exe .\geometry\channel_obstacle.geom3d
```

## VSCode

可直接使用任务：

1. `Terminal -> Run Build Task -> VS2022 Build`
2. `Terminal -> Run Task -> Run D2Q9 Example`
3. `Terminal -> Run Task -> Run Lid Driven Cavity`
4. `Terminal -> Run Task -> Run Capillary Displacement`
5. `Terminal -> Run Task -> Run Droplet Impact`
6. `Terminal -> Run Task -> Run Imported Geometry Flow`
7. `Terminal -> Run Task -> Run Imported Geometry Displacement`
8. `Terminal -> Run Task -> Run D3Q19 Cavity`
9. `Terminal -> Run Task -> Run D3 Geometry Flow`
10. `Terminal -> Run Task -> Run D3 Geometry Displacement`
11. `Terminal -> Run Task -> Run D3 Droplet Impact`
