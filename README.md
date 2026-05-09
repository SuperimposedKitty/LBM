# D2Q9 LBM C++ 示例项目

这是一个二维 D2Q9 格子玻尔兹曼方法示例项目，包含：

- 周期剪切波 `periodic shear wave`
- 顶盖驱动方腔流 `lid-driven cavity`
- 细管两相驱替 `capillary displacement`

三个案例运行完成后都只输出一个 SVG 动态图片，统一写入项目根目录下的 `result` 文件夹。
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
```

输出文件：

```text
result/periodic_shear_speed_animation.svg
result/lid_driven_cavity_speed_animation.svg
result/capillary_phase_animation.svg
```

颜色说明：

- 周期剪切波和顶盖驱动方腔：颜色表示速度大小，深灰表示固壁。
- 细管两相驱替：红色为注入相，蓝色为被驱替相，深灰为固壁。

## VSCode

可直接使用任务：

1. `Terminal -> Run Build Task -> VS2022 Build`
2. `Terminal -> Run Task -> Run D2Q9 Example`
3. `Terminal -> Run Task -> Run Lid Driven Cavity`
4. `Terminal -> Run Task -> Run Capillary Displacement`
