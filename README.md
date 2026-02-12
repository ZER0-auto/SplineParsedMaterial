# SplineParsedMaterial: 基于样条插值的MOOSE相场材料类

**该说明由deepseek根据代码自动生成。**

## 项目背景

在相场模拟中，自由能函数`f(c)`的精确表示至关重要。对于复杂自由能数据，高阶多项式拟合方法存在严重问题：
- **无法捕捉微观特征**：在浓度范围内存在多个微小的GP区势阱
- **拟合失真**：高阶多项式会扭曲关键的物理特征

## 解决方案

`SplineParsedMaterial` 将样条插值的精确性与MOOSE框架的自动导数系统相结合。
1. **样条插值保真度**：使用三次样条精确插值自由能数据
2. **自动导数计算**：继承`DerivativeMaterialInterface`，提供连续的一阶和二阶导数
3. **无缝集成**：与MOOSE现有的相场内核完全兼容

## 技术实现

### 架构设计

```cpp
SplineParsedMaterial
├── 继承: DerivativeMaterialInterface<Material>
├── 核心: SplineInterpolation (样条插值)
├── 功能: 计算 f(c), df/dc, d²f/dc²
└── 兼容: 与SplitCHParsed等内核直接集成
```

### 关键特性

1. **精确插值**：在输入数据点上完全匹配自由能值
2. **连续导数**：提供C¹连续的一阶导数和C⁰连续的二阶导数
3. **边界处理**：支持自然边界条件或指定边界导数
4. **定义域保护**：自动处理超出定义域的浓度值

## 输入参数

| 参数名 | 类型 | 必需 | 默认值 | 描述 |
|--------|------|------|---------|------|
| `x` | `std::vector<Real>` | 是 | - | 样条插值的横坐标值（浓度） |
| `y` | `std::vector<Real>` | 是 | - | 样条插值的纵坐标值（自由能） |
| `yp1` | `Real` | 否 | `1e30` | 左边界一阶导数（自然样条） |
| `ypn` | `Real` | 否 | `1e30` | 右边界一阶导数（自然样条） |
| `spline_variable` | `std::string` | 是 | - | 样条函数的变量名（如"c"） |
| `coupled_variables` | `std::vector<VariableName>` | 是 | - | 耦合的变量列表 |
| `property_name` | `std::string` | 是 | - | 材料属性名称 |
| `derivative_order` | `unsigned int` | 否 | `2` | 计算的导数阶数（最大2） |
| `enable_jit` | `bool` | 否 | `false` | JIT编译（忽略，仅为兼容性） |

## 使用示例

### 基本用法

```python
[Materials]
  [free_energy]
    type = SplineParsedMaterial
    # 样条插值数据
    x = '0.0 0.025 0.05 0.075 0.1 0.125 0.15 0.175 0.2 0.225 0.25'
    y = '0.0 0.0732865 0.120379 0.141422 0.144002 0.143766 0.146562 0.136021 0.107556 0.0643669 -0.00117124'
    # 边界导数条件
    yp1 = 0.1   # 在c=0处的导数
    ypn = 1.0   # 在c=0.25处的导数
    # 变量和属性设置
    spline_variable = c
    coupled_variables = 'c'
    property_name = F_total
    derivative_order = 2
  []
[]
```

## 与现有类型的连接

### 直接兼容的内核

1. **`SplitCHParsed`**：自动访问`F_total_c`和`F_total_c_c`导数属性
2. **`SplitCHWRes`**：使用化学势计算
3. **`CoupledTimeDerivative`**：标准时间导数项

### 辅助对象

1. **`TotalFreeEnergy`**：计算系统总自由能
2. **`DerivativeParsedMaterial`**：相似的API设计，便于迁移

## 工作机理详解

### 1. 样条插值初始化

```cpp
// 使用SplineInterpolation类
_spline.setData(x_values, y_values, yp1, ypn);
```

- **三次样条**：每段为三次多项式
- **连续性**：节点处C²连续（函数值、一阶导、二阶导连续）
- **边界条件**：可指定边界导数或使用自然样条

### 2. 导数属性声明

```cpp
// 通过DerivativeMaterialInterface注册
_dF_dc = &declarePropertyDerivative<Real>(property_name, var_name);
_d2F_dc2 = &declarePropertyDerivative<Real>(property_name, var_name, var_name);
```

这种声明方式确保MOOSE的导数系统能正确识别和访问这些属性。

### 3. 实时计算

```cpp
void SplineParsedMaterial::computeQpProperties()
{
  Real c_val = _c_val[_qp];           // 当前浓度
  _f[_qp] = _spline.sample(c_val);    // 自由能值
  (*_dF_dc)[_qp] = _spline.sampleDerivative(c_val);      // 一阶导数
  (*_d2F_dc2)[_qp] = _spline.sample2ndDerivative(c_val); // 二阶导数
}
```

### 4. 与相场内核算子的交互

```
SplitCHParsed 内核
    ↓ 请求导数
getMaterialPropertyDerivative("F_total", "c")
    ↓ 通过MOOSE导数系统
SplineParsedMaterial::_dF_dc
    ↓ 返回样条计算的导数
实时计算 df/dc
```

## 验证和测试

### 数学验证

1. **插值精度**：在输入节点处误差 < 10⁻¹⁰
2. **导数连续性**：C¹连续的一阶导数
3. **数值验证**：与有限差分法比较，误差 < 10⁻⁵

### 代码结构

```
SplineParsedMaterial/
├── SplineParsedMaterial.h    # 头文件
├── SplineParsedMaterial.C    # 源文件
├── README.md                 # 本文档
```
