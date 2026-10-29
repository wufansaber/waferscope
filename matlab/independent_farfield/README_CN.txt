独立 MATLAB 近场转远场程序

打开 abgle.m，修改开头的 USER INPUT，再点击 Run。
不依赖本项目 UI、Python、SCATMECH 或 Signal Processing Toolbox。
建议使用 MATLAB R2016b 或更新版本。此次按要求只修改，未运行验证。

最重要的输入：
1. file_name：近场 TXT/CSV 文件；找不到时会弹出文件选择框。
2. coordinate_unit：坐标的实际导出单位 m/mm/um/nm。电场须为 V/m。
3. lam：真空波长，单位 m；n_out：采样面上方均匀介质的实折射率。
4. input_exp_plus_iwt：COMSOL 数据通常设 true。
5. 8 列数据必须填写 z_plane_m；9 列则读取并检查文件中的 z。
6. reference_file：可选无颗粒参考文件；设置后先做复场的样品减参考。
   两文件必须具有相同的坐标、采样高度、单位、频率和相位参考。
   留空时不会自动去掉入射场或晶圆背景。

输入列格式（仅允许以下格式，每个文件一个频率）：
8 列：x y Ex_real Ex_imag Ey_real Ey_imag Ez_real Ez_imag
9 列：x y z Ex_real Ex_imag Ey_real Ey_imag Ez_real Ez_imag
允许带 % 或 # 表头注释，允许空格、制表符或逗号分隔实数。
只接受完整、等间距的矩形网格；重复点、缺失点、NaN/Inf 均停止处理。

transform_method：
  'fft'：二维 FFT 加线性插值，默认较快；pad_factor 仅增加频谱采样。
  'direct'：在每个方向上直接求离散平面傅里叶积分，较慢，无频谱插值。
  两种方式都是有限窗口后处理，均不会消除窗口截断或周期复制误差。

enable_window 默认 false；打开后表示对输入场加 Tukey 窗。
已移除不成立的自动能量修正；加窗后的绝对功率不是未加窗场的恢复值。
theta_max_deg 默认 85，严格小于 90。该设置不保证高角度精度。

输出：
MATLAB 工作区中的 result 以及近场文件旁 farfield_output/时间戳/ 中的文件。
farfield_result.mat：完整复 Ftheta/Fphi、角度、功率、设置及诊断。
farfield_angles.csv：逐方向的复远场振幅、归一化强度、角功率密度及 DSC。
figure_*.png：角谱、二维/三维方向图和截线。
方位角周期端点不重复；仅绘图时闭合接缝。

单位与相位：
Ftheta/Fphi：V；E_far(r) = exp(i*k*r)/r * F，工作相量为 exp(-i*omega*t)。
角功率密度：(|Ftheta|^2+|Fphi|^2)/(2*eta)，单位 W/sr。
积分功率只覆盖 theta <= theta_max_deg 的角帽，不能叫全空间总功率。
incident_intensity_w_m2 填写真实正值后，才计算 DSC [m^2/sr]；留 NaN 不推定入射强度。
observation_distance_m 留 NaN 则只输出 F；填写真实远场距离后另输出该距离处 E [V/m]。
复场的全局相位原点为坐标 (0,0,0)，程序补偿平面 xy 坐标和 z 高度相位。
若只知道相对坐标，结果也只能使用这个相对坐标系作相位参考。

HPBW 使用正负角组成的完整 XZ/YZ 切面，以该切面自身最强瓣为参考，
计算两侧半功率交点的角度差；缺少交点时报告 NaN，不外推宽度。

适用前提：采样平面在全部散射体上方，上方均匀、无源、无损，且待变换场
只向 +z 出射。下方存在基底本身不妨碍这种变换，但源场必须已正确包含其作用。
角谱圆外 K>1 是倏逝分量，不能作为同一介质中的独立远场辐射相加。
