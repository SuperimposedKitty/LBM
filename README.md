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
新增三个 Windows 原生命令，使用 OpenFOAM ASCII 字典和场文件的明确子集：

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
  system/controlDict     endTime 和 writeInterval（整数格子步数）
  constant/polyMesh/     points、faces、owner、neighbour、boundary
  0/U                   初始三分量速度及速度边界
  0/rho                 初始密度及密度边界，可用 0/p 替代
  0/alpha.A             两相时必须提供的组分密度分数
  0/porosity            可选材料孔隙度，默认 1
```

结果写到 `result/openfoam/channel2d/lbm.foam`，不会覆盖输入案例的 `0/` 和网格。重复求解仍覆盖同名结果。
不同输入目录的末级名称应不同，否则它们对应相同结果目录。

提供四个模板：

| 案例目录 | 模型 | 可修改内容 |
| --- | --- | --- |
| `cases/channel2d` | D2Q9 单相通道 | 入口速度、初始 U/rho、分辨率 |
| `cases/droplet2d` | D2Q9 封闭两相液滴 | sphereToCell 半径、相分布、湿润参数 |
| `cases/displacement3d` | D3Q27 三维驱替 | 入口相分数、boxToCell 多孔区 |
| `cases/closed3d` | D3Q19 三维封闭域 | 初始速度扰动、密度分布 |

### 输入范围与边界

- 网格必须是轴对齐、活动方向等间距的正方形/立方体体素，可含挖空的固体区域；二维为一层 z 单元，前后面必须 `empty`，z 厚度可独立设置。
- 导入原生 `polyMesh` 会根据单元中心重新映射编号，支持外部 OpenFOAM `blockMesh` 生成的多块共形均匀网格；不支持非均匀加密、斜网格、曲面单元、内部零厚度挡板、周期 patch 或任意非结构网格。矩形包围盒含填充固体格点最多 200 万个。
- 自带 `lbm_blockMesh` 支持一个标准轴对齐 `hex`、`simpleGrading (1 1 1)`、空 `edges`/`mergePatchPairs`，以及命名为 `walls`、`inlet`、`outlet`、`frontAndBack` 的边界。它不是完整 OpenFOAM blockMesh 的替代品。
- 左侧入口：网格 `patch`，U 为统一的 `fixedValue uniform (ux 0 0)`；rho（或 p）和两相 alpha.A 也用统一 fixedValue。右侧出口的 U/rho/alpha.A 为 `zeroGradient`。封闭域可以没有入口出口。外部网格的 patch 名称不限，类型和方向仍须符合上述条件。
- 静止固壁：网格 `wall`，U 为 `noSlip` 或零速度 `fixedValue`，标量为 `zeroGradient`。两相接触角通过 lbmDict 的 `contactAngle`/`wallAdhesion` 设置；尚不解析 OpenFOAM 接触角边界类型或移动壁面。原有顶盖、液滴碰撞专用 exe 保留原算法。
- `0/U`、`0/rho`（或 `0/p`）、`0/alpha.A` 支持 `uniform` 与 `nonuniform List<...>`，列表必须与输入网格单元数量一致；支持行注释、块注释、带引号的字典名。不支持 binary/gzip、宏引用、`#include`、`#codeStream` 或表达式，遇到会报错。
- 字段必须声明无量纲 `dimensions [0 0 0 0 0 0 0]`；速度、密度、时间和模型参数均用格子单位。几何坐标和间距保持输入数值，**网格尺寸缩放不自动完成物理单位换算**。
- p 是理想部分压强，转换为 `rho = 3*p`；rho 和 p 同时存在时内部值必须一致。两相由 rho 和 alpha.A 构造 rhoA/rhoB；入口 alpha.A 须严格位于 (0,1)，初始内部场允许 0 或 1。rhoA/rhoB、phase、pBulk 为派生输出，不是本入口的独立输入变量。
- `setFieldsDict` 支持 `defaultFieldValues`、`boxToCell`、`sphereToCell` 和 `volScalarFieldValue`/`volVectorFieldValue`，区域按单元中心选择，后面的区域覆盖前面的区域，保留字段的边界字典。

修改 `system/controlDict` 的 `endTime` 和 `writeInterval` 控制计算与输出；目前只支持 `startTime 0`、`deltaT 1`、`writeControl timeStep`。
修改 `system/lbmDict` 可选择 singlePhase/twoPhase、D2Q9/D3Q19/D3Q27、BGK/MRT、tau/tauA/tauB。
两相另支持 interactionStrength、bodyForce、contactAngle、wallAdhesion、recoloring、poreDiameter、darcyDrag、forchheimerDrag；单相暂不支持体力和多孔阻力。
宏观初始场被重建为平衡分布，因此本入口是初始化求解，**不是保存全部分布函数的精确断点续算**。

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
  0/{U,p,rho,porosity,rhoA,rhoB,alpha.A,alpha.B,phase,pBulk}
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
| `rho`、`p` | 总密度和理想格子压强 `p = rho / 3` |
| `porosity` | 求解器实际使用的孔隙度；单相 `P` 按普通流体处理，因此为 1 |
| `rhoA`、`rhoB` | 两相组分密度 |
| `alpha.A`、`alpha.B` | 密度分数 `rhoA/(rhoA+rhoB)` 和其补数，作为本模型的饱和度代理量 |
| `phase` | 两相序参量 `(rhoA-rhoB)/(rhoA+rhoB)` |
| `pBulk` | 两相 Shan–Chen 均匀体相压强近似，详见 LBM.md |

**全部场、坐标和时间均为无量纲格子单位**，文件 `dimensions` 为 `[0 0 0 0 0 0 0]`，时间目录名为计算步数。
没有隐含的米、秒或帕换算。不同本征密度流体的真实体积饱和度不能直接等同于密度分数；纯相中保留的微量另一组分也不会被强行截成 0 或 1。

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
