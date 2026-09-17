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

所有案例默认按输出间隔保存 OpenFOAM 网格与完整场时间序列，并在计算完成后生成 SVG 动图，统一写入项目根目录下的 `result` 文件夹。
`result` 与 `src` 同级；程序使用编译期写入的绝对路径，因此从 VSCode 运行或直接双击 exe 都会输出到同一个位置。

## OpenFOAM 格式输出

十个案例均默认启用。运行方式和 `output_interval` 不变，初始步、采样步和最后一步均保存。
每次运行生成 `result/openfoam/<案例源文件名>/<运行编号>/lbm.foam`；运行编号为微秒时间戳加冲突序号，控制台会打印完整路径。
从 VSCode 或双击 exe 启动均使用同一项目结果根目录。旧运行保留，不会混入新时间序列。

```text
result/openfoam/geometry_displacement/<运行编号>/
  lbm.foam
  constant/polyMesh/{points,faces,owner,neighbour,boundary}
  system/controlDict
  0/{U,p,rho,porosity,rhoA,rhoB,alpha.A,alpha.B,phase,pBulk}
  300/...
  600/...
```

在 ParaView 中打开 `lbm.foam`，点击 Apply，选择场并播放时间序列。三维输出包含全部体素，不限于 SVG 的中心切片。

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
原有 SVG 和显式 CSV/VTK 导出接口仍可使用。场文件在计算过程中写入，故此版本不再仅在结束时生成一个动图；中途停止时已写完的时间步仍可读取。

库调用：初始化后构造一次 writer，网格只写一次，后续快照只写场。目录必须不存在或为空；改变网格需创建新 writer。

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
build/RelWithDebInfo/bin/   exe, dll, pdb
build/RelWithDebInfo/lib/   lib, exp
```

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
.\build\RelWithDebInfo\bin\lbm_periodic.exe
.\build\RelWithDebInfo\bin\lbm_cavity.exe
.\build\RelWithDebInfo\bin\lbm_capillary.exe
.\build\RelWithDebInfo\bin\lbm_droplet.exe
.\build\RelWithDebInfo\bin\lbm_geometry_flow.exe
.\build\RelWithDebInfo\bin\lbm_geometry_displacement.exe
.\build\RelWithDebInfo\bin\lbm_d3_cavity.exe
.\build\RelWithDebInfo\bin\lbm_d3_geometry_flow.exe
.\build\RelWithDebInfo\bin\lbm_d3_geometry_displacement.exe
.\build\RelWithDebInfo\bin\lbm_d3_droplet.exe
```

输出文件：

```text
result/periodic_shear_speed_animation.svg
result/lid_driven_cavity_speed_animation.svg
result/capillary_phase_animation.svg
result/droplet_impact_animation.svg
result/geometry_flow_speed_animation.svg
result/geometry_displacement_animation.svg
result/d3_lid_driven_cavity_speed_animation.svg
result/d3_geometry_flow_speed_animation.svg
result/d3_geometry_displacement_animation.svg
result/d3_droplet_impact_animation.svg
```

颜色说明：

- 周期剪切波和顶盖驱动方腔：颜色表示速度大小，深灰表示固壁。
- 细管两相驱替：红色为注入相，蓝色为被驱替相，深灰为固壁。
- 液滴撞击：红色为液滴相，蓝色为环境相，黄色为两相界面，深灰为固壁。
- 导入几何单相流：颜色表示速度，深灰表示文件中定义的障碍物。
- 导入几何两相驱替：红色为注入相，蓝色为被驱替相，斜线区域为多孔介质。
- 三维案例在一个动态 SVG 中同时显示 `XY`、`XZ`、`YZ` 三个中心切片。

## 几何文件

几何案例默认读取 `geometry/channel_obstacle.geom`。也可以在命令行传入其他文件：

```powershell
.\build\RelWithDebInfo\bin\lbm_geometry_flow.exe .\geometry\channel_obstacle.geom
.\build\RelWithDebInfo\bin\lbm_geometry_displacement.exe .\geometry\channel_obstacle.geom
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
.\build\RelWithDebInfo\bin\lbm_d3_geometry_flow.exe .\geometry\channel_obstacle.geom3d
.\build\RelWithDebInfo\bin\lbm_d3_geometry_displacement.exe .\geometry\channel_obstacle.geom3d
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
