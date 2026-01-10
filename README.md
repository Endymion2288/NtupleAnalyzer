# Ntuple Angular Correlation Analysis

Ntuple数据的角度关联和运动学分析工具，用于分析重建后的 J/ψ+J/ψ+φ (JJP) 和 J/ψ+Υ+φ (JUP) 过程。

## 文件结构

```
NtupleAnalyzer/
├── analyze_ntuple_JJP.py    # JJP分析主程序
├── analyze_ntuple_JUP.py    # JUP分析主程序
├── plot_ntuple_results.py   # 绘图脚本
├── run_jjp_analysis.sh      # JJP运行脚本
├── run_jup_analysis.sh      # JUP运行脚本
├── run_all.sh               # 同时运行JJP和JUP
├── output/                  # 输出目录
│   ├── jjp_ntuple_correlations.root
│   ├── jup_ntuple_correlations.root
│   ├── plots_JJP/
│   └── plots_JUP/
└── README.md
```

## 数据源

### JJP (J/ψ + J/ψ + φ)
- **路径**: `/eos/user/x/xcheng/JpsiJpsiPhi_muon_pt_cut/merged_rootNtuple/`
- **Tree**: `mkcands/X_data`
- **粒子**: J/ψ₁, J/ψ₂, φ

### JUP (J/ψ + Υ + φ)
- **路径**: `/eos/user/x/xcheng/JpsiUpsPhi/merged_rootNtuple/`
- **Tree**: `mkcands/X_data`
- **粒子**: J/ψ, Υ, φ

## 事件选择 (Cuts)

### JJP 选择条件
| 变量 | 条件 |
|------|------|
| J/ψ 质量 | 2.9 - 3.3 GeV |
| φ 质量 | 0.99 - 1.10 GeV |
| J/ψ pT | > 3.0 GeV |
| φ pT | > 2.0 GeV |
| J/ψ 顶点概率 | > 0.05 |
| φ 顶点概率 | > 0.05 |
| φ K介子 pT | > 2.0 GeV |
| J/ψ Muon ID | soft (默认) |

### JUP 选择条件
| 变量 | 条件 |
|------|------|
| J/ψ 质量 | 2.9 - 3.3 GeV |
| Υ 质量 | 8.5 - 11.4 GeV |
| φ 质量 | 0.99 - 1.10 GeV |
| J/ψ pT | > 3.0 GeV |
| Υ pT | > 4.0 GeV |
| φ pT | > 2.0 GeV |
| J/ψ 顶点概率 | > 0.05 |
| Υ 顶点概率 | > 0.10 |
| φ 顶点概率 | > 0.05 |
| φ K介子 pT | > 2.0 GeV |
| J/ψ Muon ID | soft (默认) |
| Υ Muon ID | tight (默认) |

### Multi-candidate 选择
对于每个事件中的多个候选，选择 pT score 最高的候选：
```
score = sqrt(pT₁² + pT₂² + pT_φ²)
```

## 使用方法

### 快速开始

```bash
cd /eos/user/x/xcheng/CMSSW_14_0_18/src/NtupleAnalyzer

# 运行JJP分析
./run_jjp_analysis.sh

# 运行JUP分析
./run_jup_analysis.sh

# 同时运行两个分析
./run_all.sh
```

### JJP 分析

```bash
# 基本用法
python3 analyze_ntuple_JJP.py

# 限制处理事件数（用于测试）
python3 analyze_ntuple_JJP.py -n 10000

# 指定输出文件
python3 analyze_ntuple_JJP.py -o my_output.root

# 修改Muon ID要求
python3 analyze_ntuple_JJP.py --muon-id tight
python3 analyze_ntuple_JJP.py --muon-id none  # 不做Muon ID要求
```

### JUP 分析

```bash
# 基本用法
python3 analyze_ntuple_JUP.py

# 限制处理事件数
python3 analyze_ntuple_JUP.py -n 10000

# 自定义Muon ID
python3 analyze_ntuple_JUP.py --jpsi-muon-id soft --ups-muon-id tight
python3 analyze_ntuple_JUP.py --jpsi-muon-id none --ups-muon-id none
```

### 绘图

```bash
# JJP绘图
python3 plot_ntuple_results.py -i output/jjp_ntuple_correlations.root -o plots_JJP -p JJP

# JUP绘图
python3 plot_ntuple_results.py -i output/jup_ntuple_correlations.root -o plots_JUP -p JUP
```

## 输出直方图

### 角度关联
- `h_dy_*`: Δy (快度差) 分布
- `h_dphi_*`: Δφ (方位角差) 分布
- `h2_dy_dphi_*`: 2D Δy vs Δφ 关联图

### 运动学分布
- `h_*_pt`: pT 分布
- `h_*_eta`: η (伪快度) 分布
- `h_*_y`: y (快度) 分布
- `h_*_phi`: φ (方位角) 分布

### 不变质量
- `h_mass_*`: 两体/三体不变质量分布

## 输出图像

每个过程生成以下图像：
- `delta_y_comparison_*.pdf/png`: Δy 比较图
- `delta_phi_comparison_*.pdf/png`: Δφ 比较图
- `correlation_2d_*_*.pdf/png`: 单独的2D关联图
- `correlation_2d_all_*.pdf/png`: 组合2D关联图
- `pt_distributions_*.pdf/png`: pT 分布
- `eta_distributions_*.pdf/png`: η 分布
- `rapidity_distributions_*.pdf/png`: y 分布
- `phi_distributions_*.pdf/png`: φ 分布
- `invariant_mass_*.pdf/png`: 不变质量分布

## 参考

- 事件选择参考: `NPSExtraction/sPlotFit/run_mass_fit_JJP.py`, `run_mass_fit_JUP.py`
- 角度关联计算参考: `JJPMCAnalyzer/analyze_gen_correlations.py`
- 绘图风格参考: `JJPMCAnalyzer/plot_results.py`
