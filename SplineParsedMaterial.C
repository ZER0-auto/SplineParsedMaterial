//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "SplineParsedMaterial.h"

// 请注意替换为你的项目名称+App
registerMooseObject("testApp", SplineParsedMaterial);

InputParameters
SplineParsedMaterial::validParams()
{
  InputParameters params = DerivativeMaterialInterface<Material>::validParams();

  // 添加样条特定参数
  params.addRequiredParam<std::vector<Real>>("x", "Abscissa values for spline interpolation");
  params.addRequiredParam<std::vector<Real>>("y", "Ordinate values for free energy f(c)");
  params.addParam<Real>("yp1", 1e30,
    "First derivative at left boundary (natural spline if not specified)");
  params.addParam<Real>("ypn", 1e30,
    "First derivative at right boundary (natural spline if not specified)");

  // 匹配你的输入文件参数
  params.addRequiredParam<std::string>("spline_variable", "The variable for spline interpolation");
  params.addRequiredCoupledVar("coupled_variables", "The coupled variables");

  // 属性名称参数 - 匹配你的输入文件
  params.addRequiredParam<std::string>("property_name", "Name of the material property");

  // 导数阶数 - 匹配你的输入文件
  params.addParam<unsigned int>("derivative_order", 2, "Maximum order of derivatives to compute");

  // enable_jit参数（暂时不实现，先忽略）
  params.addParam<bool>("enable_jit", false, "Enable JIT compilation (not implemented yet)");

  params.addClassDescription("Material that defines free energy using spline interpolation");

  return params;
}

SplineParsedMaterial::SplineParsedMaterial(const InputParameters & parameters)
  : DerivativeMaterialInterface<Material>(parameters),
    _x_values(getParam<std::vector<Real>>("x")),
    _y_values(getParam<std::vector<Real>>("y")),
    _c_val(coupledValue("coupled_variables")),
    _property_name(getParam<std::string>("property_name")),
    _var_name(coupledName("coupled_variables", 0)),
    _f(declareProperty<Real>(_property_name)),
    _derivative_order(getParam<unsigned int>("derivative_order")),
    _dF_dc(nullptr),
    _d2F_dc2(nullptr),
    _x_min(_x_values.front()),
    _x_max(_x_values.back())
{
  // 获取边界条件
  Real yp1 = getParam<Real>("yp1");
  Real ypn = getParam<Real>("ypn");

  // 验证输入数据
  if (_x_values.size() != _y_values.size())
    paramError("y", "x and y arrays must have the same size");

  if (_x_values.size() < 2)
    paramError("x", "At least two data points are required for spline interpolation");

  // 检查单调性
  for (size_t i = 1; i < _x_values.size(); ++i)
  {
    if (_x_values[i] <= _x_values[i-1])
      paramError("x", "x values must be strictly increasing");
  }

  // 设置样条数据
  _spline.setData(_x_values, _y_values, yp1, ypn);

  // 检查spline_variable参数是否与coupled_variables匹配
  std::string spline_var_name = getParam<std::string>("spline_variable");

  if (spline_var_name != _var_name)
  {
    mooseWarning("spline_variable ('", spline_var_name,
                 "') does not match the first coupled_variable ('",
                 _var_name, "'). Using the coupled variable.");
  }

  // 声明导数属性 - 使用DerivativeMaterialInterface的declarePropertyDerivative
  if (_derivative_order >= 1)
  {
    _dF_dc = &declarePropertyDerivative<Real>(_property_name, _var_name);
  }

  if (_derivative_order >= 2)
  {
    _d2F_dc2 = &declarePropertyDerivative<Real>(_property_name, _var_name, _var_name);
  }

  // 打印样条信息用于验证
  Moose::out << "SplineParsedMaterial initialized:" << std::endl;
  Moose::out << "  Property name: " << _property_name << std::endl;
  Moose::out << "  Spline variable: " << spline_var_name << std::endl;
  Moose::out << "  Coupled variable: " << _var_name << std::endl;
  Moose::out << "  Domain: [" << _x_min << ", " << _x_max << "]" << std::endl;
  Moose::out << "  Number of data points: " << _x_values.size() << std::endl;
  Moose::out << "  Derivative order: " << _derivative_order << std::endl;

  // 打印导数属性名
  if (_derivative_order >= 1 && _dF_dc)
  {
    Moose::out << "  First derivative property declared via DerivativeMaterialInterface" << std::endl;
  }

  if (_derivative_order >= 2 && _d2F_dc2)
  {
    Moose::out << "  Second derivative property declared via DerivativeMaterialInterface" << std::endl;
  }

  // 打印样条数据用于调试
  if (_x_values.size() <= 20)
  {
    Moose::out << "  X values: ";
    for (size_t i = 0; i < _x_values.size(); ++i)
    {
      Moose::out << _x_values[i];
      if (i < _x_values.size() - 1) Moose::out << ", ";
    }
    Moose::out << std::endl;

    Moose::out << "  Y values: ";
    for (size_t i = 0; i < _y_values.size(); ++i)
    {
      Moose::out << _y_values[i];
      if (i < _y_values.size() - 1) Moose::out << ", ";
    }
    Moose::out << std::endl;
  }
}

Real
SplineParsedMaterial::computeValue(Real c) const
{
  // 确保值在样条定义域内
  if (c < _x_min || c > _x_max)
  {
    // 只在第一次出现时警告，避免大量输出
    static bool warned = false;
    if (!warned && _t_step == 0)
    {
      mooseWarning("Value ", c, " outside spline domain [",
                   _x_min, ", ", _x_max,
                   "]. Clamping to domain boundaries.");
      warned = true;
    }
    c = std::max(_x_min, std::min(_x_max, c));
  }

  return _spline.sample(c);
}

Real
SplineParsedMaterial::computeDerivative(Real c, unsigned int order) const
{
  // 确保值在样条定义域内
  if (c < _x_min || c > _x_max)
  {
    c = std::max(_x_min, std::min(_x_max, c));
  }

  switch (order)
  {
    case 0:
      return _spline.sample(c);
    case 1:
      return _spline.sampleDerivative(c);
    case 2:
      return _spline.sample2ndDerivative(c);
    default:
      // 对于三次样条，三阶及以上导数为0
      return 0.0;
  }
}

void
SplineParsedMaterial::computeQpProperties()
{
  // 获取当前积分点的变量值
  Real c_val = _c_val[_qp];

  // 计算函数值
  Real f_val = computeValue(c_val);
  _f[_qp] = f_val;

  // 计算并存储导数
  if (_derivative_order >= 1 && _dF_dc)
  {
    Real df_dc = computeDerivative(c_val, 1);
    (*_dF_dc)[_qp] = df_dc;
  }

  if (_derivative_order >= 2 && _d2F_dc2)
  {
    Real d2f_dc2 = computeDerivative(c_val, 2);
    (*_d2F_dc2)[_qp] = d2f_dc2;
  }

  // 验证输出（只在第一个时间步的第一个积分点）
  if (_qp == 0 && _t_step == 0)
  {
    Moose::out << "=== SPLINE CALCULATION ===" << std::endl;
    Moose::out << "At QP " << _qp << ":" << std::endl;
    Moose::out << "  c = " << c_val << std::endl;
    Moose::out << "  f(c) = " << f_val << std::endl;

    if (_dF_dc)
    {
      Real df_dc = (*_dF_dc)[_qp];
      Moose::out << "  df/dc = " << df_dc << std::endl;

      // 数值导数验证
      Real eps = 1e-6;
      Real f_plus = computeValue(c_val + eps);
      Real f_minus = computeValue(c_val - eps);
      Real num_df = (f_plus - f_minus) / (2 * eps);
      Moose::out << "  Numerical df/dc = " << num_df << std::endl;
      Moose::out << "  Difference = " << std::abs(df_dc - num_df) << std::endl;
    }

    if (_d2F_dc2)
    {
      Real d2f_dc2 = (*_d2F_dc2)[_qp];
      Moose::out << "  d2f/dc2 = " << d2f_dc2 << std::endl;

      // 数值二阶导数验证
      Real eps = 1e-6;
      Real f_plus = computeValue(c_val + eps);
      Real f_minus = computeValue(c_val - eps);
      Real f_center = f_val;
      Real num_d2f = (f_plus - 2 * f_center + f_minus) / (eps * eps);
      Moose::out << "  Numerical d2f/dc2 = " << num_d2f << std::endl;
      Moose::out << "  Difference = " << std::abs(d2f_dc2 - num_d2f) << std::endl;
    }
  }

  // 额外的调试信息：检查SplitCHParsed能否正确访问导数
  if (_qp == 0 && _t_step == 0)
  {
    Moose::out << "=== DERIVATIVE PROPERTY CHECK ===" << std::endl;
    Moose::out << "  These properties should be accessible by SplitCHParsed:" << std::endl;
    Moose::out << "  - " << _property_name << " (free energy)" << std::endl;
    if (_dF_dc)
      Moose::out << "  - Derivative of " << _property_name << " w.r.t. " << _var_name << std::endl;
    if (_d2F_dc2)
      Moose::out << "  - Second derivative of " << _property_name << " w.r.t. " << _var_name << std::endl;
  }
}
