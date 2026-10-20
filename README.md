# tilt：NiO₆ 八面体倾转分析（Windows）

输入 **POSCAR / CONTCAR**，自动生成每个 NiO₆ 的倾角、带正负号的投影角、三维转动分量、Ni—O—Ni 桥角、逐层关系以及 PNG/SVG 图。

这是 Python 源码运行包 **v1.0.1**，不是免安装 exe。支持 Windows 10/11 x64 + Python 3.11；不需要 VASP、Linux、服务器连接或 Git。首次安装依赖需要网络，安装后结构分析可离线进行。

## 最简单的使用方法

1. 下载 [Windows 源码运行包](dist/tilt-windows-v1.0.1.zip)，点击文件页的 **Download raw file**。也可点仓库 **Code → Download ZIP**。
2. **完整解压**到桌面等有写权限的位置，不要直接在压缩包里运行。
3. 安装 [Python 3.11 的 Windows 64-bit 安装包](https://www.python.org/downloads/release/python-3119/)，安装时包含 **Python Launcher**，并勾选 **Add python.exe to PATH**。已有 Python 3.11 x64 可跳过。此包固定依赖版本，不建议直接换用 Python 3.13/3.14。
4. 双击 **`setup_windows.bat`**。它在本文件夹创建 `.venv`，安装依赖并运行测试。看到 `Installation and tests completed` 即完成。
5. 双击 **`demo_windows.bat`** 跑人工样例；结果放在 `results` 下带时间戳的新文件夹。
6. 处理自己的结构：把 **POSCAR/CONTCAR 拖到 `run_windows.bat`**；或双击它，在文件选择框里选结构。中文路径及空格路径受支持。

每次运行创建新结果目录，结束后窗口会保留输出。不会覆盖输入或既有结果。只选择普通文本 POSCAR/CONTCAR，不要选择 ZIP、VESTA 文件或 OUTCAR。

默认入口采用**自动面内参考方向**。比较同一组构型时应固定参考轴；常规 Amam 晶胞通常用45°，√2超胞或其他晶胞不要机械照搬。

## 推荐的科研用法：固定参考方向

在解压目录打开 PowerShell，运行：

```powershell
.\.venv\Scripts\python.exe -X utf8 .\ni_octahedra.py "D:\结构\POSCAR" --azimuth 45 --out ".\results\case01"
```

将 `D:\结构\POSCAR` 改成自己的文件路径。`case01` 必须是新目录；第二个结构可使用 `case02`。也可用桌面入口：

```powershell
.\run_windows.bat "D:\结构\CONTCAR" --azimuth 45
```

自动推断参考方向：

```powershell
.\.venv\Scripts\python.exe -X utf8 .\ni_octahedra.py "D:\结构\POSCAR" --out ".\results\auto01"
```

**不需要激活虚拟环境，也不需要修改 PowerShell 执行策略**；直接调用 `.venv\Scripts\python.exe` 即可，见 [Python 官方 venv 文档](https://docs.python.org/3.11/library/venv.html)。

## 结果怎样看

| 文件 | 内容 |
|---|---|
| `octahedra.csv` | 逐 Ni 的倾角、x/y/z 转动分量、拟合残差和警告；可用 Excel 打开 |
| `tilt_yz.png` / `.svg` | 侧视投影：+y向右、+z向上，红色为顺时针、蓝色为逆时针，标注POSCAR中的Ni编号 |
| `layers.csv` | 分层、层内均值/RMS；有正负交替时不要只看均值 |
| `bridges.csv` | 真正共享氧的 Ni—O—Ni 三维角 |
| `interlayer_pairs.csv` | 相邻层最近配准候选，同/反向、近零及配准多解信息 |
| `analysis.json` | 全部数据，包括源哈希、参考坐标、周期近邻和旋转矩阵 |

`apical_tilt_deg` 是两个顶点 O 连线相对层法向的夹角。`omega_x/y/z_deg` 是六键最佳刚体拟合的旋转向量，**不是三个 Euler 角**，也不一定等于顶点轴倾角。两者差异可以来自八面体形变。

默认假设 **ab 为层面、c 为堆垛方向**。五/七配位不会强行拟合成 NiO₆；强畸变和邻居歧义会标记。层间关系不是自动的 Glazer `+`/`−`，程序不强行认证空间群或给所有缺陷超胞贴 Glazer 标签。自动面内参考会吸收统一绕z的转动，定量比较请固定 `--azimuth`。

详细公式、论文关系和限制：[方法说明](docs/METHOD.md)。

## 人工样例与测试

`examples/POSCAR_rx8` 是人工构建的孤立理想 NiO₆：在10 Å立方盒中绕+x右手旋转8°，不是La₃Ni₂O₇的DFT结果。使用 `--azimuth 0` 应得到：

- 六键旋转向量约 `(8, 0, 0)` 度；
- 顶点轴倾角 `8°`；
- yz投影顺时针角 `−8°`；
- 六键拟合残差接近零。

手动运行测试：

```powershell
.\.venv\Scripts\python.exe -X utf8 -m unittest discover -s . -p "test_*.py" -v
```

测试包括已知旋转、周期像、非正交晶胞、缺氧/过配位、中文与空格路径、UTF-8 BOM、PNG/SVG生成、重复输出拒绝和入口错误码。实际 Windows 自动测试见仓库 [Actions](https://github.com/LHCstat/tilt/actions/workflows/windows-tests.yml)；绿色通过表示该提交的测试通过，不能代替你自己结构的物理验收。

## 常见问题

- **找不到 Python**：安装 Python 3.11 x64（含Launcher）。确认 `py -3.11 --version` 或 `python --version` 显示3.11。
- **下载依赖失败**：保留安装窗口里的报错，检查网络，恢复后重新运行setup；不用重新下载程序。依赖固定为NumPy1.26.4、SciPy1.13.1、ASE3.26.0、Matplotlib3.9.4。
- **找不到文件**：路径要加双引号。PowerShell的当前目录须是解压后的程序目录。
- **输出目录已经存在**：换一个 `--out` 名称，不要覆盖以前结果。
- **没有完整NiO₆ / 角度为空**：检查邻居和结构；`--cutoff 2.6` 可调，但不能为了“凑6个”而掩盖真实配位变化。
- **没有文件选择窗口**：直接把文件拖到bat，或使用上述命令行。
- **搬动了整个程序目录后不能运行**：虚拟环境不适合直接搬动；在新位置重新解压程序并运行setup创建环境。

源码、测试和人工样例可直接查看。包内不包含Python本体、第三方依赖安装包、论文PDF或原始DFT计算目录；无需POTCAR/WAVECAR/CHGCAR。论文依据：Chen等 [Adv. Mater. 2025](https://doi.org/10.1002/adma.202504238)、Dong等 [Nat. Mater. 2025](https://doi.org/10.1038/s41563-025-02351-2)。
