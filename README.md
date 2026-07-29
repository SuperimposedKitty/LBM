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

所有案例运行完成后都只输出一个 SVG 动态图片，统一写入项目根目录下的 `result` 文件夹。
`result` 与 `src` 同级；程序使用编译期写入的绝对路径，因此从 VSCode 运行或直接双击 exe 都会输出到同一个位置。

## 构建

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

数值测试不会生成结果文件：

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
