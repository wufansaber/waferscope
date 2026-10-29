# SCATMECH 与 FEM 的复场相位参考

日期：2026-10-26，接触修订版 r2。

## 两个不同的问题

几何接触由 `delta` 或 `gap_nm` 定义。本版平整晶圆默认二者均为 0，球心高度 h=d/2，球底始终在晶圆面 z=0。

复场还需要约定相位原点。Bobbert–Vlieger 球谐展开围绕球心定义；FEM 的平面波和体源辐射使用晶圆表面坐标。
直接把原生 JonesDSC 输入成像会混用参考原点，即使角向强度和收集功率相同，倾斜通道的像面位置及峰值像素也可能不同。

## 转换

本软件采用 exp(-iωt) 约定。球心高度 h=radius+delta，入射单位方向为 (0,sinθi,-cosθi)，观测方向为 ns。
成像前使用几何平移相位：

```text
A_wafer(ns) = A_BV(ns) * exp[-i * (2*pi/lambda) * h * (cos(theta_i) + ns_z)]
```

这是入射平面波在球心处的相位与出射远场平移相位之积，不是拟合参数，不会移动球体，不会改变 NA 或材料，也不会缩放能量。
相位因子的模恒为 1，因此不改变 DSC、总收集功率和 s/p 权重；只修正复场参考。
有镀膜时 z=0 指最上层界面，不穿过膜层移动到底层基底。

## 验证

- 单元测试核对几何相位公式及单位模。
- 独立自由空间 Rayleigh 点偶极子检验：使用已知位置的辅助小球，比较 NIST 复远场与解析偶极子平移场，只消除一个不可观测的全局常相位，不拟合强度、位置或角相关相位。
- 10 nm 接触球与 FEM 的复光瞳对照中，边路误差由约 1.84% 降至约 0.014%；这不表示全部尺寸的 FEM 均已收敛。
- 积分功率保持不变；100 nm 剩余功率差异仍需要网格、阶数、计算域和 PML 收敛检查。

## 数据保存

`jones_dsc_um` 保留库原生、球心参考的 JonesDSC，不覆盖原始数据。
`sphere_to_wafer_phase` 单独保存实际施加的相位数组。
`pupil_modes`、`image_e_modes`、`image_hz0_modes` 使用晶圆顶面参考；结果记录 `field_reference` 和球心高度。
旧结果不自动转换或改写；需要用新版本重新计算 SCATMECH 相干成像。

同时优化了 FEM 的严格零源项情况：仅当组装后的右端向量每个元素都恰为零时跳过矩阵分解，返回零散射场；非零源项仍执行原有求解和残差检查。它不省略背景镜面场，也不关闭 haze。

参考：[NIST 模型参数](https://pages.nist.gov/SCATMECH/docs/bobvlieg.htm)、[官方球心展开及 qq=k(radius+delta) 源码](https://github.com/usnistgov/SCATMECH/blob/master/code/bobvlieg1.cpp)、[官方入射/出射相位实现](https://github.com/usnistgov/SCATMECH/blob/master/code/bobvlieg2.cpp)。
