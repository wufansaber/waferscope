# 暗场颗粒近场、远场与三路成像的理论模型

版本日期：2026-09-27。对应 `fullwave_prototype`、`optical_propagation.py` 和
源码版 `Full-wave FEM (experimental)` 页面。

本文区分三个层面：连续电磁方程、实际离散算法、已完成的数值验证。
“满足 Maxwell 方程的建模形式”不等于“任意参数下已经得到精确结果”。
完整的一般颗粒加粗糙面问题通常没有有限项闭式解；本程序采用解析平面背景、
有限元全波散射解、解析背景 Green 函数的远场提取和理想矢量成像传播。

## 1. 范围与几何

当前 FEM 支持一个均匀、各向同性、非磁性球形颗粒，空气上方、均匀半无限晶圆下方，
以及平面上的有限三角粗糙面。颗粒可以与平面分离；非零粗糙面时必须高于最高表面点。
复杂材料的折射率由用户给定，不自动查询材料色散表。

坐标与单位：
- 晶圆参考面为 `z=0`，空气为 `z>0`。
- 激光从负 KY 一侧入射，传播方向具有正 Y 和负 Z 分量。
- 入射角 `theta_i` 相对于法线，不是相对于晶圆平面。
- 计算长度用 um；界面上的波长、粒径、间隙以 nm 输入并转换。
- 三个物镜沿绕 Y 轴旋转的三个方向分别收集，默认倾角为 +45、-45、0 度。
- 默认物镜焦距 20 mm、倍率 20，对应管镜焦距 400 mm；三个像面分别计算。

默认示例的颗粒 n=1.59、k=0 只是模型输入，并非覆盖所有波长的 PSL 实测材料模型。
特别是改成 266 nm 时，真实 PSL、晶圆和膜层的色散及吸收需要单独提供。

## 2. 时间、介质与场的约定

采用时间因子

$$
\mathbf E(\mathbf r,t)=\operatorname{Re}\{\mathbf E(\mathbf r)e^{-i\omega t}\}.
$$

被动介质写成 $\tilde n=n+i\kappa$，$\kappa\ge0$，相对介电常数

$$
\epsilon_r=(n+i\kappa)^2,
\qquad k_0=2\pi/\lambda_0,
\qquad \mu_r=1.
$$

FEM 输出使用归一化场

$$
\mathbf e=\mathbf E/E_0,\qquad
\mathbf h=Z_0\mathbf H/E_0.
$$

这里的磁场数组名称为 `hz0`，不是只有 Z 分量，而是整个矢量 $Z_0\mathbf H/E_0$。
入射空气平面波的照度为 $I_0=|E_0|^2/(2Z_0)$。
因而归一化 Poynting 通量乘以 $I_0$ 就得到物理照度：

$$
\langle\mathbf S\rangle/I_0=\operatorname{Re}(\mathbf e\times\mathbf h^*).
$$

不能把近场 $|\mathbf e|^2$ 无条件解释为净传输能量；驻波和倏逝场中两者并不等价。

## 3. 解析平面背景

令 $s=\sin\theta_i$、$c=\cos\theta_i$，则

$$
\hat{\mathbf k}_{i}=(0,s,-c),\qquad
\hat{\mathbf k}_{r}=(0,s,c),\qquad
\mathbf k_t/k_0=(0,s,-q),
\quad q=\sqrt{\epsilon_s-s^2}.
$$

平方根选择向下传播/衰减的被动分支：$\operatorname{Im}q\ge0$，无损时取非负实部。
Fresnel 系数为

$$
r_s=\frac{c-q}{c+q},\qquad
r_p=\frac{\epsilon_s c-q}{\epsilon_s c+q}.
$$

s 基为 $(-1,0,0)$；p 基按 $\hat{\mathbf p}=\hat{\mathbf k}\times\hat{\mathbf s}$ 定义。
这使反射 p 基与入射 p 基的切向分量符号相反。比较不同程序的 Jones 系数时必须
先统一基矢，不能看到 $r_p$ 符号不同就直接判断物理错误。

空气侧背景是入射场与反射场之和，晶圆内是透射场。背景在平面上满足切向 E/H 连续。
`planar.py` 提供其数值表达，`solver.py` 构造相同的有限元系数函数。

## 4. 全波散射方程与有限元

令 $\epsilon_b(z)$ 是平面参考介质，实际介质为 $\epsilon(\mathbf r)$，
$\Delta\epsilon=\epsilon-\epsilon_b$。分解总场为背景与散射场：

$$
\mathbf e=\mathbf e_b+\mathbf e_s,
\qquad
\nabla\times\nabla\times\mathbf e_s-k_0^2\epsilon\mathbf e_s
=k_0^2\Delta\epsilon\mathbf e_b.
$$

这是散射场形式的完整线性 Maxwell 方程，**不是把总场近似成入射场的一阶 Born 近似**。
实际系数 $\epsilon$ 出现在左端，求解包含颗粒与表面的多次相互作用。

连续弱形式可写为

$$
\int_\Omega [ (\nabla\times\mathbf e_s)\cdot(\nabla\times\overline{\mathbf v})
-k_0^2\epsilon\mathbf e_s\cdot\overline{\mathbf v}]\,dV
=k_0^2\int_\Omega\Delta\epsilon\mathbf e_b\cdot\overline{\mathbf v}\,dV.
$$

NGSolve 使用实形函数展开复系数，代码以复对称双线性矩阵组装；它不是 Hermitian
正定问题，不能因求解器名称带 `cholesky` 就按实正定能量形式理解。
采用切向连续的 Nedelec $H(\mathrm{curl})$ 元，允许介质界面上法向 E 跳变。
球面采用曲边网格，场的阶数与几何曲率阶数分别处理。[NGSolve H(curl) 文档](https://docu.ngsolve.org/v6.2.1708/i-tutorials/unit-2.3-hcurlhdiv/hcurlhdiv.html)

### 4.1 静态凝聚

将内部自由度 L 与共享自由度 E 分开：

$$
S=A_{EE}-A_{EL}A_{LL}^{-1}A_{LE},\qquad
b'_E=b_E-A_{EL}A_{LL}^{-1}b_L.
$$

先解 $Sx_E=b'_E$，再恢复 $x_L=A_{LL}^{-1}(b_L-A_{LE}x_E)$。
这改变计算成本，不改变离散方程。代码恢复内部场后才提取近场与复极化源。
[NGSolve 静态凝聚说明](https://ngsolve.org/ngsolve/docs/i-tutorials/unit-1.4-staticcond/staticcond.html)

### 4.2 开放边界与误差

有限区域外使用径向复坐标拉伸 PML：

$$
\tilde r=r+i\alpha T[(r-R)/T]^p,\qquad R<r<R+T.
$$

切向 Jacobian 特征值为 $\tilde r/r$，径向特征值为 $d\tilde r/dr$。
PML 通过复坐标变换同时作用于 curl 与质量项，而不是仅添加一个经验损耗常数。
最外边界将散射场切向分量置零。径向变换保持 `z=0` 参考界面。
[NGSolve PML 文档](https://docu.ngsolve.org/latest/how_to/pml.html)

连续 PML 的理想无反射性质不消除离散误差。拉伸越强也不保证越准确：过快衰减可能
无法由现有单元解析；太弱则会产生有限厚度截断反射。应联合检查网格、阶数、厚度与强度。
`integration_bonus` 用于独立改变组装积分阶数；没有实测改善时不将其当作精度升级。
小代数残差只表示离散线性系统求解一致，不能代替物理误差评估。
改变 PML 的目标单元尺寸会重新生成整张网格，可能连带改变物理区域的网格拓扑；
不能把这种对比描述为“只替换 PML 单元、其他单元完全不变”。

## 5. 颗粒、凹凸粗糙面与相干背景

有限三角表面 $z=h(x,y)$ 之外保持平面。高度按实际三角面的**投影面积**调整至
零均值和给定 RMS，而不是仅对节点值做标准差归一化。

参考介质上的三类体积对比度为

$$
\Delta\epsilon=
\begin{cases}
\epsilon_p-1 & \text{空气中的球},\\
\epsilon_s-1 & \text{参考平面上方凸起},\\
1-\epsilon_s & \text{参考平面下方凹坑}.
\end{cases}
$$

该构造可描述当前三角高度场的凹凸形状，但不是任意凹形缺陷建模器，也不包含悬垂结构。
同一张网格分别求有颗粒 $E_1$ 和无颗粒 $E_0$。定义 $\Delta E=E_1-E_0$ 后有

$$
|E_1|^2=|E_0|^2+|\Delta E|^2+2\operatorname{Re}(\Delta E E_0^*).
$$

像面的严格通量分解使用 E 与 H 的交叉项。颗粒导致的场增量功率不等于
“有颗粒图减无颗粒图”，后者包含可正可负的相干干涉项。
有限斑块、单个随机种子不代表无限晶圆粗糙度的统计 BRDF；需另做面积与样本数收敛。

## 6. 从近场到远场

从曲边有限元单元积分提取体极化源

$$
\mathbf P_c(\mathbf r')=\Delta\epsilon(\mathbf r')\mathbf e(\mathbf r'),
\qquad \mathbf p_j=\mathbf P_c(\mathbf r_j)\,w_j.
$$

这里 $\mathbf P_c$ 已除去 $\epsilon_0E_0$，$w_j$ 包含实际单元 Jacobian。
使用总场而不是背景场，保证多次散射已包含在源项中。

在均匀空气中，远场系数为

$$
\mathbf A(\hat{\mathbf k})=\frac{k_0^2}{4\pi}
(I-\hat{\mathbf k}\hat{\mathbf k}^{T})
\int\mathbf P_c(\mathbf r')e^{-ik_0\hat{\mathbf k}\cdot\mathbf r'}dV',
\quad
\mathbf e_s\sim\frac{e^{ik_0r}}r\mathbf A.
$$

半空间中，空气侧源还包含一次由**解析背景 Green 函数**产生的反射传播项。
令 $\hat k_d=(k_x,k_y,-k_z)$，$s=(-k_y,k_x,0)/\sqrt{k_x^2+k_y^2}$，
$p_u=\hat k\times s$、$p_d=\hat k_d\times s$，则上方源的两个振幅为

$$
A_s=\frac{k_0^2}{4\pi}s\cdot\int_{z'>0}\mathbf P_c
[e^{-ik_0\hat k\cdot r'}+r_s e^{-ik_0\hat k_d\cdot r'}]dV',
$$

$$
A_p=\frac{k_0^2}{4\pi}\int_{z'>0}
[p_u e^{-ik_0\hat k\cdot r'}+r_p p_d e^{-ik_0\hat k_d\cdot r'}]\cdot\mathbf P_c\,dV'.
$$

下方凹坑源通过互易的透射背景波提取：相位为
$\exp[-ik_0(k_xx'+k_yy'+qz')]$，s 测试矢量为 $(1+r_s)s$；
p 测试矢量的横向分量为 $(1-r_p)p_{u,xy}$，法向分量为
$-2(k_xp_{u,x}+k_yp_{u,y})/(\epsilon_s k_z+q)$。
这些矢量与复源做不取复共轭的互易内积，不能误写成 Hermitian 内积。
`far_field.py` 实现该表达；上、下方源在振幅层面相加。

这个“传播时的一次反射”不是只保留颗粒与表面一次相互作用：反复相互作用已由全波总场求出。
远场 DSC 为

$$
\frac{d\sigma}{d\Omega}=|A_s|^2+|A_p|^2\quad[\mathrm{um}^2/\mathrm{sr}],
\qquad P_{\mathcal C}=I_0\int_{\mathcal C}\frac{d\sigma}{d\Omega}d\Omega.
$$

显示的半球图是 DSC 关于方向的图，不等于均匀平面面积上的功率密度。
若投影到 $u=k_x,v=k_y$ 平面，则 $d\Omega=du\,dv/k_z$，积分时必须保留 Jacobian。
当前散射图不包含无限平面背景的镜面反射光束；物镜若收入镜面方向，也不会自动加入该光束。

## 7. Jones 与 Mueller 的位置

对于确定性线性散射系统，Jones 矩阵映射两个入射偏振复振幅至两个出射复振幅：
$\mathbf A=J\mathbf a$。本 FEM 一次只解一个 s 或 p 输入，不把一次解冒充完整 Jones 矩阵。
若要构成完整 J，需要两次独立入射偏振解，并保持共同相位约定。

非相干输入用相干矩阵 $C=\langle\mathbf a\mathbf a^\dagger\rangle$ 描述，
$C_{out}=JCJ^\dagger$。选定 Stokes/Pauli 约定后

$$
M_{ij}=\tfrac12\operatorname{Tr}(\sigma_iJ\sigma_jJ^\dagger).
$$

Mueller 强度信息不能一般地恢复空间相位；不能从 $\sqrt{\mathrm{DSC}}$ 随意补一个
零相位再声称得到正确近场或成像。当前 FEM 链路始终保留复源和复光瞳。

## 8. 三路收集、几何率与能量率

第 j 路的轴 $a_j=(\sin\beta_j,0,\cos\beta_j)$，空气中半角
$\alpha_j=\arcsin(\mathrm{NA}_j)$，上半球收集区域为

$$
\mathcal C_j=\{\hat k:k_z>0,\ \hat k\cdot a_j\ge\cos\alpha_j\}.
$$

NA=0 表示禁用。能量占比与几何占比分别是

$$
\eta_E=\frac{\int_{\cup_j\mathcal C_j}\mathrm{DSC}\,d\Omega}
{\int_{k_z>0}\mathrm{DSC}\,d\Omega},
\qquad
\eta_G=\frac{\int_{\cup_j\mathcal C_j}d\Omega}{2\pi}.
$$

并集只计一次；但三个独立物镜的像面不能据此相加成一幅相干图。
采用共同球半径 R 截断的锥体体积 $V=R^3\Omega/3$，其与半球体积的比值才等于上述几何率。
若用不同高度的直圆锥体积，含义不同。未被地平线裁剪的单圆锥有
$\Omega=2\pi(1-\cos\alpha)$，裁剪或重叠时需做区域积分。
本节给出区分原则，不表示 FEM 实验页已增加新的三路并集统计控件。

## 9. 从复光瞳到三路像面

每路采用局部方向 $(u,v,w)$，$w=\sqrt{1-u^2-v^2}$。先将远场 s/p 复振幅转为
笛卡尔矢量，再用旋转将该方向的横向场准直到局部光轴，保持偏振和矢量长度。
定义平面瞳函数

$$
\mathbf P(u,v)=\frac{R_{\hat k\to\hat z}\mathbf A(\hat k)}{\sqrt w},
\qquad\int|\mathbf P|^2du\,dv=\int|\mathbf A|^2d\Omega.
$$

检偏器单位矢量 a 作用为 $P\mapsto a(a^TP)$；功率透过率 T 作用为 $P\mapsto\sqrt T P$。
倍率 M 下，理想像方方向为

$$
\hat k' =(-u/M,-v/M,w'),\qquad w'=\sqrt{1-(u^2+v^2)/M^2}.
$$

像方局部矢量由准直横向瞳场旋转回球面方向，并除以 $\sqrt{w'}$ 保持轴向通量。
略去统一传播相位，像面场为

$$
\mathbf e_{img}(x',y',z')=\frac1{\lambda_0 M}
\int \frac{R_{\hat z\to\hat k'}\mathbf P}{\sqrt{w'}}
e^{-i\frac{2\pi}{\lambda_0 M}(ux'+vy')}
e^{ik_0z'(w'-1)}du\,dv.
$$

磁场将积分中的矢量替换为 $\hat k'\times\mathbf e_{angular}$，
像面照度计算为 $I_0\operatorname{Re}(e_xh_y^*-e_yh_x^*)$。
不能在高 NA 情形简单以 $|E_x|^2+|E_y|^2$ 替代轴向能流。
这种理想消球差、正弦条件的矢量衍射结构对应 Richards–Wolf 框架；
程序不包含真实透镜表面、像差或有限镜筒遮挡。

FFT 是上述角谱积分的离散实现。中点光瞳网格需要半个采样间距的相位修正。
零填充改善像面采样间隔，但不增加原始光瞳信息，不消除周期视野/边缘截断误差。
高 NA 方形网格四角即使场为零，也必须先给出有限的占位方向，避免 `0 * NaN` 污染 FFT。

### 9.1 不通过 FFT 的独立成像基准

均匀 x 偏振圆瞳可先解析积分方位角，令 $r\in[0,a]$、$s=r/M$、
$w=\sqrt{1-s^2}$、$b=2\pi xs/\lambda$、$W=2\pi r e^{ikz(w-1)}/(\lambda M\sqrt w)$：

$$
e_x(x,0,z)=\int_0^a W[J_0(b)-\tfrac{s^2}{2(1+w)}(J_0(b)-J_2(b))]dr,
$$
$$
e_z(x,0,z)=-i\int_0^a WsJ_1(b)dr,
\quad
h_y(x,0,z)=\int_0^a W[J_0(b)-\tfrac{s^2}{2(1+w)}(J_0(b)+J_2(b))]dr.
$$

其他有关分量在该截线上由对称性为零。该一维 Bessel 积分不用 FFT、不调用旋转函数，
可以分别检验复 E、复 H 和强度；近轴极限再退化到标量 Airy 图样。

## 10. 独立球体理论、守恒与镜像解

### 10.1 Mie 系数与截面

均匀球在均匀空气中，令 $x=k_0a$、$m=n+i\kappa$、
$\psi_l(z)=zj_l(z)$、$\xi_l(z)=zh_l^{(1)}(z)$：

$$
a_l=\frac{m\psi_l(mx)\psi_l'(x)-\psi_l(x)\psi_l'(mx)}
{m\psi_l(mx)\xi_l'(x)-\xi_l(x)\psi_l'(mx)},
\quad
b_l=\frac{\psi_l(mx)\psi_l'(x)-m\psi_l(x)\psi_l'(mx)}
{\psi_l(mx)\xi_l'(x)-m\xi_l(x)\psi_l'(mx)}.
$$
$$
C_{sca}=\frac{2\pi}{k_0^2}\sum_{l=1}^\infty(2l+1)(|a_l|^2+|b_l|^2),
\quad
C_{ext}=\frac{2\pi}{k_0^2}\sum_{l=1}^\infty(2l+1)\operatorname{Re}(a_l+b_l).
$$

代码中独立使用 SciPy 球 Bessel 函数计算截面，并与 `miepython` 的递推算法交叉检查。
此直接 Bessel 小型验证器限于本文的小尺寸参数，不替代稳定的大尺寸 Mie 专用算法。
球内/球外的复 E/H 则由独立 `miepython.field` 的矢量球谐展开提供。
转换入射坐标和球心相位是已知几何变换，不拟合幅度或相位。
报告中的近场相对 L2 是各球壳采样点上不加权的离散范数，不是体积积分范数，
也不是界面处最大误差；远场的球面积分另使用 Gauss-Legendre 权重。
该库吸收折射率输入采用 n-ik，适配层明确转换，避免把吸收写成增益。

小球 Rayleigh 极限为

$$
\alpha_v=4\pi a^3\frac{\epsilon_p-1}{\epsilon_p+2},
\qquad C_{sca}\sim\frac{k_0^4}{6\pi}|\alpha_v|^2.
$$

其适用条件包含 $x\ll1$ 及材料内部尺度要求；不能把本式强行用于 100 nm/266 nm 后
把差异全部归咎于全波程序。[Mie 算法公式](https://miepython.readthedocs.io/en/latest/07_algorithm.html)

### 10.2 光学定理与能量守恒

在均匀无损背景下，按本程序的振幅约定有

$$
C_{ext}=\frac{4\pi}{k_0}\operatorname{Im}[e_i^*\cdot A(\hat k_i)],
\quad
C_{abs}=k_0\int\operatorname{Im}\epsilon_p\,|e|^2dV,
\quad C_{ext}=C_{sca}+C_{abs}.
$$

另外可以从极化源得到 $C_{ext}=k_0\operatorname{Im}\int e_i^*\cdot\Delta\epsilon e\,dV$。
远场前向值与该体积分共享同一组源，彼此一致主要检验归一化/提取实现，
不能把机器精度的一致性当作 FEM 物理误差已经达到机器精度。
真实独立检查还包括理论 Mie 值、近场解析展开、Poynting 球面通量和网格收敛。
有损半空间不能直接套用这个自由空间截面守恒式，遗漏进入晶圆的能量会得到错误结论。

### 10.3 理想导体镜面极限

电偶极子 $p=(p_x,p_y,p_z)$ 位于 $(x,y,h)$ 时，理想导体面上的镜像为
$p'=(-p_x,-p_y,p_z)$，位置 $(x,y,-h)$。自由空间偶极子与镜像的辐射之和
可直接给出上半球复远场，不调用 Fresnel 传播核。
程序将介质吸收参数逐渐增大，检验半空间核是否趋向这一镜像解，分别检查水平和垂直偶极子。
这是理想导体极限检查，不表示实际硅片是理想导体，也不是对任意有限介电常数的完整证明。

### 10.4 闭合面等效电磁流：第二条近远场路径

对于包围源、且完全位于均匀空气中的闭合面 S，散射场给出等效电流与磁流
$j_s=n\times h_s$、$m_s=-n\times e_s$。令
$\widetilde j=\int_S j_s e^{-ik_0\hat k\cdot r'}dS'$、
$\widetilde m=\int_S m_s e^{-ik_0\hat k\cdot r'}dS'$，则

$$
A(\hat k)=\frac{ik_0}{4\pi}
[(I-\hat k\hat k^T)\widetilde j-\hat k\times\widetilde m].
$$

`huygens.py` 使用这一表面积分，不使用体积介电对比度，也不调用原来的体源辐射核。
将解析偶极子 E/H 输入该路径，可独立检验符号和绝对幅度；将 FEM 球面 E/H 输入，
可以与体源提取结果比较。该实现仅用于均匀空气球体基准，不能直接拿一个跨越有损晶圆的
球面套用自由空间 Green 函数。跨介质界面需相应半空间 Green 函数或分域处理。
依据：[Purdue 电磁场等效原理讲义](https://engineering.purdue.edu/wcchew/ece604f23/EMFTEDX012324r.pdf)。

## 11. 论文依据与验证证据的层级

1. Bobbert, P. A. and Vlieger, J. (1986), *Light scattering by a sphere on a substrate*,
   Physica A 137, 209–242, DOI: [10.1016/0378-4371(86)90072-5](https://doi.org/10.1016/0378-4371(86)90072-5)。
   本项目通过 [NIST 官方模型说明](https://pages.nist.gov/SCATMECH/docs/bobvlieg.htm)
   和 [官方实现](https://pages.nist.gov/SCATMECH/code/bobvlieg1.cpp) 使用其平面球体参照，
   并需检查多极展开截断；不宣称复现论文全部图表或实验。
2. Richards, B. and Wolf, E. (1959), *Electromagnetic diffraction in optical systems II*,
   Proc. R. Soc. A 253, 358–379, DOI: [10.1098/rspa.1959.0200](https://doi.org/10.1098/rspa.1959.0200)。
   本轮检验的是本文明确列出的同类矢量积分及其 Bessel 化简，不是对付费全文每个公式逐页复核。
3. [miepython 官方场 API](https://miepython.readthedocs.io/en/latest/miepython.html)：
   独立球内外复场实现；[归一化说明](https://miepython.readthedocs.io/en/latest/03a_normalization.html)
   用于避免将单位积分相函数错当 DSC。
4. [Purdue 电磁学讲义：Rayleigh/Mie 与光学定理](https://engineering.purdue.edu/wcchew/ece604f19/Lecture%20Notes/Lect34.pdf)。
5. [Princeton 偶极子镜像问题](https://kirkmcd.princeton.edu/examples/image_antenna.pdf)，
   [UT Austin 偶极子辐射推导](https://farside.ph.utexas.edu/teaching/em/lectures/node95.html)。
   前一链接本轮检索可见但全文抓取失败；镜像测试公式在本文独立列出，可直接复算。

公开文档、独立代码、解析极限和数值收敛提供不同种类的证据；它们都不能代替真实样品测量。
材料色散、颗粒形貌、膜层、污染、照明光束和真实镜头像差的不确定性可能超过数值误差。

本轮实测数字、未通过试验和最终采用的设置见 [理论交叉校验结果](THEORY_VALIDATION_CN.md)。
关于二维采样峰值、界面距离和局部网格加密的区别，见 [近场收敛诊断](NEAR_CONVERGENCE_CN.md)。

## 12. 可选的分区多项式阶数

为降低直接求解的内存成本，可仅降低外部吸收层的有限元阶数。设物理半径为 R，
PML 厚度为 T，保护比例为 f，单元顶点平均位置为 c_K，则待降阶单元满足
`material(K) 属于 PML` 且 `|c_K| > R + f*T`。物理区和其他吸收层单元
保持原阶数；与保护单元相连的边、面也保持原阶数。这是相容 H(curl) 空间的
数值离散调整，不是材料、边界方程或几何的修改。

默认关闭降阶。独立 Mie 测试表明，降低吸收层阶数会增加反射/截断相关误差；
局部物理区加密不能自动补偿这部分误差。保护比例、网格尺寸、PML 阶数必须分别做
敏感性检查，不能仅凭残差小或像面功率守恒判断物理解准确。
实际内存、Mie 误差及未通过算例见 [低内存加密报告](LOW_MEMORY_REFINEMENT_CN.md)。

## 13. 后续扩展

2026-09-28 新增可选矢量有限角谱照明、出瞳 Zernike 相位、相机响应和复出瞳缓存。
默认关闭光束与像差扩展，保留本文件原有平面波、理想物镜基准。
新增模型的公式、参数单位、独立校验和未完成研究项见 [模型扩展说明](MODEL_UPGRADE_CN.md)。
