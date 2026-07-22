# D2Q9 LBM C++ 示例项目

这是一个二维 D2Q9 格子玻尔兹曼方法示例项目，包含：

- 周期剪切波 `periodic shear wave`
- 顶盖驱动方腔流 `lid-driven cavity`
- 细管两相驱替 `capillary displacement`
- 液滴撞击固体表面并反弹 `droplet impact`
- 文件导入几何的单相障碍流 `imported geometry flow`
- 文件导入几何的两相驱替 `imported geometry displacement`

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

## 运行

```powershell
.\build\RelWithDebInfo\bin\lbm_periodic.exe
.\build\RelWithDebInfo\bin\lbm_cavity.exe
.\build\RelWithDebInfo\bin\lbm_capillary.exe
.\build\RelWithDebInfo\bin\lbm_droplet.exe
.\build\RelWithDebInfo\bin\lbm_geometry_flow.exe
.\build\RelWithDebInfo\bin\lbm_geometry_displacement.exe
```

输出文件：

```text
result/periodic_shear_speed_animation.svg
result/lid_driven_cavity_speed_animation.svg
result/capillary_phase_animation.svg
result/droplet_impact_animation.svg
result/geometry_flow_speed_animation.svg
result/geometry_displacement_animation.svg
```

颜色说明：

- 周期剪切波和顶盖驱动方腔：颜色表示速度大小，深灰表示固壁。
- 细管两相驱替：红色为注入相，蓝色为被驱替相，深灰为固壁。
- 液滴撞击：红色为液滴相，蓝色为环境相，黄色为两相界面，深灰为固壁。
- 导入几何单相流：颜色表示速度，深灰表示文件中定义的障碍物。
- 导入几何两相驱替：红色为注入相，蓝色为被驱替相，斜线区域为多孔介质。

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

## VSCode

可直接使用任务：

1. `Terminal -> Run Build Task -> VS2022 Build`
2. `Terminal -> Run Task -> Run D2Q9 Example`
3. `Terminal -> Run Task -> Run Lid Driven Cavity`
4. `Terminal -> Run Task -> Run Capillary Displacement`
5. `Terminal -> Run Task -> Run Droplet Impact`
6. `Terminal -> Run Task -> Run Imported Geometry Flow`
7. `Terminal -> Run Task -> Run Imported Geometry Displacement`
