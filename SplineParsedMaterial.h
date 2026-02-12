//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "DerivativeMaterialInterface.h"
#include "Material.h"
#include "SplineInterpolation.h"

/**
 * Material that uses spline interpolation for free energy function
 */
class SplineParsedMaterial : public DerivativeMaterialInterface<Material>
{
public:
  static InputParameters validParams();
  SplineParsedMaterial(const InputParameters & parameters);

protected:
  virtual void computeQpProperties() override;

  // 使用样条计算函数值
  virtual Real computeValue(Real c) const;

  // 使用样条计算导数
  virtual Real computeDerivative(Real c, unsigned int order) const;

private:
  // 样条插值对象
  SplineInterpolation _spline;

  // 存储插值数据
  std::vector<Real> _x_values;
  std::vector<Real> _y_values;

  // 变量值
  const VariableValue & _c_val;

  // 属性名称
  std::string _property_name;

  // 变量名
  std::string _var_name;

  // 材料属性
  MaterialProperty<Real> & _f;

  // 导数阶数
  unsigned int _derivative_order;

  // 导数属性
  MaterialProperty<Real> * _dF_dc;
  MaterialProperty<Real> * _d2F_dc2;

  // 定义域边界
  Real _x_min;
  Real _x_max;
};
