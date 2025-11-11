# 问题解决总结

## 问题描述

用户在尝试使用 `TS_twosample` 函数时遇到以下错误：

```r
Error in TS_twosample(data1 = data1, data4 = data4, tao = tao, band = 180,  : 
  unused argument (plot = TRUE)
```

## 根本原因

包中存在**两个文件**定义了同一个函数 `TS_twosample`：

1. **R/Functions.R** (第134行) - 新版本
   - 包含完整的参数列表：`alpha`, `plot`, `plot_genes`, `gene_names`, `condition_names`
   - 包含可视化功能
   - 使用 Rcpp 加速的 C++ 函数

2. **R/Functions_old.R** (第138行) - 旧版本  
   - 只有基本参数：`data1`, `data4`, `tao`, `band`, `quant`
   - 没有可视化功能
   - 使用纯 R 实现

当 R 包加载时，会按照某种顺序加载 `R/` 目录下的所有 `.R` 文件。**由于两个文件都包含 `@export` 标签**，旧版本的函数定义覆盖了新版本，导致用户无法使用新功能。

## 解决方案

### 步骤 1: 重命名旧文件

```bash
mv R/Functions_old.R R/Functions_old.R.bak
```

将 `Functions_old.R` 重命名为 `Functions_old.R.bak`，这样 R 就不会自动加载它（R 只加载 `.R` 结尾的文件）。

### 步骤 2: 重新生成文档

```r
devtools::document()
```

这会：
- 重新生成 NAMESPACE 文件（只包含新版本的导出）
- 更新函数文档 (.Rd 文件)
- 清除旧的函数定义

### 步骤 3: 重新编译和加载

```r
devtools::clean_dll()      # 清理旧的编译文件
pkgbuild::compile_dll()    # 重新编译 C++ 代码
devtools::load_all()       # 加载更新后的包
```

## 验证修复

运行以下代码验证问题已解决：

```r
devtools::load_all()
data(data1)
data(data4)

tao <- est.c(data1, data4, max1 = 5, max4 = 5)

# 这应该能正常工作了
results <- TS_twosample(
  data1 = data1,
  data4 = data4,
  tao = tao,
  band = 180,
  quant = c(0.01, 0.01, 0.01),
  plot = TRUE  # ✓ 不再报错
)
```

## 测试结果

```
✓ Function executed successfully!
Number of significant genes: 3 
Significant gene indices: 1 2 5
```

## 文件变更

| 文件 | 变更 | 说明 |
|------|------|------|
| `R/Functions_old.R` | 重命名为 `.bak` | 防止旧版本被加载 |
| `FINAL_DEMO.R` | 已更新 | 使用新的参数和功能的完整示例 |
| `HOW_TO_USE.md` | 已更新 | 添加问题诊断和解决方案说明 |
| `NAMESPACE` | 已重新生成 | 只导出新版本函数 |
| `man/TS_twosample.Rd` | 已重新生成 | 更新函数文档 |

## 预防措施

为了避免将来出现类似问题：

1. **不要在多个文件中定义同名函数**
2. **旧代码应该：**
   - 移到 `inst/` 目录
   - 或者移除 `@export` 标签
   - 或者重命名为非 `.R` 扩展名（如 `.R.bak` 或 `.R.old`）

3. **开发时定期检查：**
   ```r
   # 检查函数来源
   environment(TS_twosample)
   
   # 检查函数参数
   args(TS_twosample)
   ```

## 后续建议

1. **完全移除旧文件**（如果不再需要）：
   ```bash
   rm R/Functions_old.R.bak
   ```

2. **或者移到归档目录**：
   ```bash
   mkdir -p inst/archive
   mv R/Functions_old.R.bak inst/archive/Functions_old.R
   ```

3. **更新 README.md** 说明包现在使用 Rcpp 加速

4. **运行完整的包检查**：
   ```r
   devtools::check()
   ```

## 技术细节

### R 包加载顺序

R 包会按照以下顺序处理函数：
1. 读取 `R/` 目录下所有 `.R` 文件（按文件名字母顺序）
2. 执行每个文件中的代码
3. 如果有重复的函数名，**后加载的会覆盖先加载的**
4. `@export` 标签决定哪些函数会被导出到 NAMESPACE

### 为什么 `devtools::load_all()` 也失败了？

`devtools::load_all()` 会：
1. 重新加载所有 R 代码
2. 但**不会**重新生成 NAMESPACE（除非你运行 `document()`）
3. 仍然会按顺序加载所有 `.R` 文件

所以即使使用 `load_all()`，旧版本仍会覆盖新版本。

## 总结

✅ **问题已完全解决**

用户现在可以：
- 使用 `plot = TRUE` 参数生成可视化
- 使用 `alpha` 参数设置显著性水平
- 使用 `gene_names` 和 `condition_names` 自定义标签
- 使用 `plot_genes` 指定要详细绘制的基因
- 访问新的返回值：`p_values`, `significant_genes`, `z_critical`

所有功能已测试并正常工作！
