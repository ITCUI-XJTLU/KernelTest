# KernelTest Package Documentation Index

Complete guide to all documentation files in the KernelTest package.

## 🚀 Quick Start

**New to the package?** Start here:
1. Read [README.md](README.md) - Package overview and quick start
2. Run [FINAL_DEMO.R](FINAL_DEMO.R) - See the package in action
3. Check [HOW_TO_USE.md](HOW_TO_USE.md) if you encounter issues

## 📚 Documentation Files

### User Documentation

| File | Purpose | Language | Audience |
|------|---------|----------|----------|
| [README.md](README.md) | Package overview, installation, quick start | English | All users |
| [FINAL_DEMO.R](FINAL_DEMO.R) | Complete working example with visualization | R code | All users |
| [HOW_TO_USE.md](HOW_TO_USE.md) | Detailed usage guide and troubleshooting | 中文 | All users |
| [SOLUTION_SUMMARY.md](SOLUTION_SUMMARY.md) | Problem diagnosis and solutions | 中文 | Troubleshooting |
| [CHANGELOG.md](CHANGELOG.md) | Version history and updates | English | All users |

### Technical Documentation

| File | Purpose | Audience |
|------|---------|----------|
| [CLAUDE.md](CLAUDE.md) | Code architecture, development guide | Developers |
| [src/kernel_core.cpp](src/kernel_core.cpp) | C++ implementation source code | Developers |
| [R/Functions.R](R/Functions.R) | R wrapper functions source code | Developers |

### R Package Documentation

| Command | Documentation |
|---------|---------------|
| `?TS_twosample` | Two-sample kernel test function |
| `?est.c` | Bias term estimation |
| `?TS_kernel` | Single-sample kernel test |
| `?NormTransformation` | Variance stabilization |

## 🔍 Find What You Need

### I want to...

#### Install the package
→ See [README.md - Installation](README.md#installation)

#### Run a basic analysis
→ Run [FINAL_DEMO.R](FINAL_DEMO.R) or see [README.md - Quick Start](README.md#quick-start)

#### Fix "unused argument (plot = TRUE)" error
→ See [SOLUTION_SUMMARY.md](SOLUTION_SUMMARY.md) or [HOW_TO_USE.md](HOW_TO_USE.md)

#### Generate plots
→ See [README.md - With Automatic Visualization](README.md#with-automatic-visualization-)

#### Customize plots with gene names
→ See [README.md - Advanced Visualization Options](README.md#advanced-visualization-options)

#### Understand the algorithm
→ See [CLAUDE.md - Core Implementation](CLAUDE.md#core-implementation)

#### Improve performance
→ See [README.md - Performance](README.md#performance)

#### Contribute to development
→ See [CLAUDE.md - Development Commands](CLAUDE.md#development-commands)

#### See what's new
→ See [CHANGELOG.md](CHANGELOG.md)

## 📖 Documentation by Topic

### Installation & Setup
- [README.md - Installation](README.md#installation)
- [README.md - Troubleshooting Installation](README.md#troubleshooting-installation)
- [SOLUTION_SUMMARY.md - 解决方案](SOLUTION_SUMMARY.md)

### Basic Usage
- [README.md - Quick Start](README.md#quick-start)
- [FINAL_DEMO.R](FINAL_DEMO.R)
- [HOW_TO_USE.md - 快速开始](HOW_TO_USE.md)

### Visualization
- [README.md - Visualization Features](README.md#visualization-features)
- [FINAL_DEMO.R - Step 2 & 4](FINAL_DEMO.R)
- [README.md - Advanced Visualization Options](README.md#advanced-visualization-options)

### Troubleshooting
- [SOLUTION_SUMMARY.md](SOLUTION_SUMMARY.md)
- [HOW_TO_USE.md - 故障排除](HOW_TO_USE.md)
- [README.md - Troubleshooting Installation](README.md#troubleshooting-installation)

### Development
- [CLAUDE.md](CLAUDE.md)
- [CHANGELOG.md](CHANGELOG.md)
- [src/kernel_core.cpp](src/kernel_core.cpp)

## 🎯 Common Tasks

### Task: Run a complete analysis

```r
# Step 1: Load the package
devtools::load_all()

# Step 2: Load data
data(data1)
data(data4)

# Step 3: Estimate bias term
tao <- est.c(data1, data4, max1 = 5, max4 = 5)

# Step 4: Run analysis with visualization
results <- TS_twosample(
  data1 = data1,
  data4 = data4,
  tao = tao,
  band = 180,
  quant = c(0.01, 0.01, 0.01),
  plot = TRUE,
  alpha = 0.05
)

# Step 5: Examine results
print(results$significant_genes)
print(results$p_values)
```

**Learn more:** [FINAL_DEMO.R](FINAL_DEMO.R)

### Task: Troubleshoot installation issues

1. Check if `R/Functions_old.R` exists → should be renamed to `.bak`
2. Run `devtools::document()`
3. Run `devtools::load_all()`
4. Test: `args(TS_twosample)` should show 10 parameters

**Learn more:** [SOLUTION_SUMMARY.md](SOLUTION_SUMMARY.md)

### Task: Customize visualizations

```r
results <- TS_twosample(
  data1, data4, tao, band = 180, quant = c(0.01, 0.01, 0.01),
  plot = TRUE,
  plot_genes = c(1, 3, 5),  # Which genes to plot
  gene_names = c("BRCA1", "TP53", "EGFR", "MYC", "KRAS"),
  condition_names = c("Normal", "Tumor"),
  alpha = 0.01  # Stricter significance level
)
```

**Learn more:** [README.md - Advanced Visualization Options](README.md#advanced-visualization-options)

## 📝 File Descriptions

### README.md
Main documentation file with:
- Package overview and features
- Installation instructions
- Quick start guide
- Basic and advanced usage examples
- Performance benchmarks
- Visualization guide
- Recent updates and future plans

### FINAL_DEMO.R
Executable R script demonstrating:
- Complete workflow from data loading to visualization
- Basic and advanced visualization options
- Custom gene names and condition labels
- Examining results and statistics

### HOW_TO_USE.md (中文)
Chinese language guide covering:
- Problem diagnosis (问题诊断)
- Solution steps (解决方案)
- Quick start (快速开始)
- Advanced usage (高级用法)
- Troubleshooting (故障排除)

### SOLUTION_SUMMARY.md (中文)
Detailed technical solution guide:
- Root cause analysis (根本原因)
- Step-by-step solution (解决方案)
- Verification steps (验证修复)
- File changes (文件变更)
- Prevention measures (预防措施)

### CHANGELOG.md
Version history documenting:
- New features and enhancements
- Performance improvements
- Bug fixes
- API changes
- Migration guide

### CLAUDE.md
Developer documentation with:
- Package architecture
- Development workflow
- Code organization
- Build and test commands
- Implementation details

### src/kernel_core.cpp
C++ source code containing:
- Gaussian kernel implementation
- Nadaraya-Watson smoother
- Hat matrix computation
- Variance estimation
- Eigenvalue statistics

## 🆘 Getting Help

### Quick Help
```r
# View function help
?TS_twosample

# Check function parameters
args(TS_twosample)

# List all available functions
help(package = "KernelTest")
```

### File-based Help
- Installation issues → [SOLUTION_SUMMARY.md](SOLUTION_SUMMARY.md)
- Usage questions → [HOW_TO_USE.md](HOW_TO_USE.md)
- Examples → [FINAL_DEMO.R](FINAL_DEMO.R)
- API reference → [README.md](README.md)

### Development Help
- Architecture questions → [CLAUDE.md](CLAUDE.md)
- Build issues → [CLAUDE.md - Development Commands](CLAUDE.md#development-commands)
- Performance → [README.md - Performance](README.md#performance)

## 📊 Documentation Statistics

- Total documentation files: 8
- Lines of documentation: ~1500+
- Languages: English, 中文, R
- Code examples: 20+
- Screenshots/visualizations: Built-in plots

## 🔄 Last Updated

- README.md: 2025-01-11
- FINAL_DEMO.R: 2025-01-11
- HOW_TO_USE.md: 2025-01-11
- SOLUTION_SUMMARY.md: 2025-01-11
- CHANGELOG.md: 2025-01-11
- DOCUMENTATION_INDEX.md: 2025-01-11

---

**Tip:** Bookmark this file for quick access to all documentation resources!
