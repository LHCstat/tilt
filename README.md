# tilt：NiO₆ 八面体倾转分析（Windows）

输入 **POSCAR / CONTCAR**，自动生成每个 NiO₆ 的倾角、带正负号的投影角、三维转动分量、Ni—O—Ni 桥角、逐层关系以及 PNG/SVG 图。

这是 Python 源码运行包 **v1.1.0**，不是免安装 exe。支持 Windows 10/11 x64 + Python 3.11；不需要 VASP、Linux、服务器连接或 Git。首次安装依赖需要网络，安装后结构分析可离线进行。

## 最简单的使用方法

1. 下载 [Windows 源码运行包](dist/tilt-windows-v1.1.0.zip)，点击文件页的 **Download raw file**。也可点仓库 **Code → Download ZIP**。
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
| `tilt_xz.png` / `.svg` | xz投影：+x向右、+z向上；与yz图同时生成 |
| `tilt_yz_slices.png` / `.svg`、`tilt_xz_slices.png` / `.svg` | 检出投影重合时额外生成的前后两片并排图 |
| `projection_slices.json` | 分片范围、每片 Ni 编号及分片前后仍重合的 Ni 对 |
| `layers.csv` | 分层、层内均值/RMS；有正负交替时不要只看均值 |
| `bridges.csv` | 真正共享氧的 Ni—O—Ni 三维角 |
| `interlayer_pairs.csv` | 相邻层最近配准候选，同/反向、近零及配准多解信息 |
| `analysis.json` | 全部数据，包括源哈希、参考坐标、周期近邻和旋转矩阵 |

`apical_tilt_deg` 是两个顶点 O 连线相对层法向的夹角。`omega_x/y/z_deg` 是六键最佳刚体拟合的旋转向量，**不是三个 Euler 角**，也不一定等于顶点轴倾角。两者差异可以来自八面体形变。

两张图使用实际六个氧投影的凸包作为半透明轮廓，并叠加真实 Ni—O 十字键线。每个 Ni 旁保留原子编号和一个角度数值，例如 `−6.24°`；该数值为当前 yz / xz 投影相对 +z 的有符号角，上端向右为正。图中不再绘制 z 参考虚线、圆弧或三维总倾角文字。完整精度和三维总倾角仍保存在 CSV/JSON 中。xz 投影为零不代表三维总倾角为零；畸变结构的轮廓不一定是理想菱形。

默认假设 **ab 为层面、c 为堆垛方向**。五/七配位不会强行拟合成 NiO₆；强畸变和邻居歧义会标记。层间关系不是自动的 Glazer `+`/`−`，程序不强行认证空间群或给所有缺陷超胞贴 Glazer 标签。自动面内参考会吸收统一绕z的转动，定量比较请固定 `--azimuth`。

对于 √2 超胞，同一投影中的 Ni 中心距离小于 0.6 Å 时，会额外输出两片并排图。yz 图沿分数坐标 a 分片，xz 图沿分数坐标 b 分片；周期性地划为两个半胞，自动选取边界以减少片内重合，并尽量远离 Ni。每个独立 Ni 只归属一片，所属八面体完整绘出；启用边界补齐时可额外绘制它的周期像。图上写明分数坐标范围，跨越 1→0 的范围属于同一片；非正交晶胞中这些是晶体学分片，不完全等同于笛卡尔等深度分片。各结构独立选择边界，Slice 1/2 不表示跨结构固定的原子对应关系。

两片共享坐标范围、比例和角度配色，完整叠加图同时保留。0.6 Å 只用于识别 Ni 中心的近重合，不是键长或物理判据；复杂超胞可能仍有轮廓相交，或需要更多切片。若两片中仍有 Ni 中心近重合，图顶和 `projection_slices.json` 会提示。

### sqrt2 的统一方向、边界补齐和间隙氧

命令行和桌面入口默认 `--display auto`：对 La24Ni16O(56+δ) 且含四个 Ni 层的结构，统一以底部 Ni 的三原子边为 x、两原子边为 y。v0 对应 x 平行输入 b；右手坐标系的 y 指向输入 −a（非正交时取层面内与 x 垂直方向）。其他结构按底层排列识别 a/b；如果原点截断双层，会先将 z 原点移到双层之间的宽间隙。`--display native` 可使用原坐标轴与原点；`--display lno` 强制上述规则，排列不明确时会报错。

显示范围对应截图的分数坐标：x/y 为 −0.05～1.05，z 为 −0.005～1.005。分片图在屏幕横向和 z 方向补入落在范围内的 Ni 周期像，其完整八面体随中心绘出；观察深度方向每个独立 Ni 只选一个代表，因此不增加第三片。补出的像沿用编号并加 `′`。底部在 xz 图中为两行三列，上部在 yz 图中为两行三列；上下两组面内排列错开，正交投影中另一组为两行两列，不人为复制非边界位置的原子。

统一方向是对分析副本做晶格基矢重排和原点平移，不改变周期结构，也不覆盖输入 POSCAR。变换矩阵和原点记录在 `analysis.json` 的 `display` 中；显示的每个原子像及整数平移记录在 `projection_slices.json` 中。投影角与旋转分量会随新坐标轴重新计算，三维顶点轴倾角不变。`--azimuth` 此时相对于统一后的 x 轴，不能直接沿用旧坐标的数值。

间隙氧用绿色菱形和原始 O 编号标记：sqrt2 要求 16 个 Ni 均为完整六配位、骨架共使用 56 个不同 O；La12Ni8 单胞则要求 8 个完整六配位 Ni、28 个骨架 O。只有满足对应骨架判据，才把其余过量 O 识别为间隙氧，不依赖“排在文件最后”的假设。单胞沿用输入晶轴，也会绘制间隙氧标记。不满足骨架判据时记录 `unresolved`，不将所有未配位 O 强行认定为间隙氧。周期像不计入间隙氧数量；v0 无间隙氧。

v1.1.0 保留 v1.0.1 的近邻识别、六配位判断、六键拟合和倾角计算算法；新增坐标约定和展示处理。`run_windows.bat` 输出结构 PNG/SVG 及 CSV/JSON；仅输入 POSCAR 不能得到总能量、形成能或相对 v0 的晶格变化，桌面入口不自动生成带这些信息的 PDF 图集。版本记录见 [CHANGELOG.md](CHANGELOG.md)。

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
