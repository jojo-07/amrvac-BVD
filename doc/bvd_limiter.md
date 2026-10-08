# BVD limiter 实现文档

本文档介绍 MPI-AMRVAC 中新增的 **BVD（Boundary Variation Diminishing，边界变差递减）限制器**的算法原理、代码实现、以及 AMRVAC 宏（VAC 预处理器）系统的使用方法。面向不熟悉 AMRVAC 宏系统的二次开发者。

---

## 目录

1. [BVD 算法概述](#1-bvd-算法概述)
2. [与 AMRVAC 的索引约定映射](#2-与-amrvac-的索引约定映射)
3. [AMRVAC 宏系统（VAC 预处理器）详解](#3-amrvac-宏系统vac-预处理器详解)
4. [代码走读：mod_bvd.t](#4-代码走读mod_bvdt)
5. [集成点：所有修改的文件](#5-集成点所有修改的文件)
6. [使用方法](#6-使用方法)
7. [编译与测试](#7-编译与测试)
8. [验证结果](#8-验证结果)
9. [二次开发指南](#9-二次开发指南)

---

## 1. BVD 算法概述

### 1.1 核心思想

BVD 是一种**单元剖面择优**的重构方法。对每个单元 i，构造一条完整剖面 q_i(x)，剖面同时给出单元**两端**的重构值。

**命名约定**：为避免与界面的"左/右"混淆，单元端用 **1/2** 编号（与原程序 `poly_q(i,1:2)` 一致）：

- **2 端**（右端）q_2(i)：面 i+1/2 的左侧值
- **1 端**（左端）q_1(i)：面 i-1/2 的右侧值

注意视角翻转：单元的"1 端"落在界面的**右**侧、单元的"2 端"落在界面的**左**侧。

然后按**单元总边界变差 TBV** 为该单元选择剖面——被选剖面的两端**一起采用**，保证同一单元的两个界面值来自同一候选（剖面自洽）。

```
TBV(q, i) = |q_2(i) - q_1(i+1)| + |q_2(i-1) - q_1(i)|
```

即单元 i 两个界面处的**跨面跳变**之和（经典 BVD 定义，Sun et al. 2016）。TBV 越小表示剖面越"平缓"。

### 1.2 候选剖面

本实现有两个（可选三个）候选：

**候选 1：二次多项式**（3 点二次插值，对应 `second_degree_polynomial`）

```
2 端 = (-1·U_{i-1} + 5·U_i + 2·U_{i+1}) / 6
1 端 = ( 2·U_{i-1} + 5·U_i - 1·U_{i+1}) / 6
```

**候选 2：THINC**（Tangent of Hyperbola for Interface Capturing，双曲正切界面捕捉）

单元内剖面是一条**移位 tanh 曲线**：

```
q(x) = umin + 0.5·umax·(1 + r·tanh(βξ + atanh(A)))
```

其中：

- `umin = min(U_{i-1}, U_{i+1})`，`umax = |U_{i+1} - U_{i-1}|`（左右邻居定义允许区间）
- `r = sign(U_{i+1} - U_{i-1})`（方向符号）
- `cc = (U_i - umin) / umax`（归一化单元平均）
- **A 由单元平均守恒条件解出**（保证剖面平均值严格等于 U_i）：

```
B = exp(r·β·(2·cc - 1))
A = (B / cosh(β) - 1) / tanh(β)
1 端 = umin + 0.5·umax·(1 + r·A)
2 端 = umin + 0.5·umax·(1 + r·(tanh(β) + A)/(1 + A·tanh(β)))
```

THINC 的关键性质：

1. **保单调**：两端值恒落在 `[umin, umin+umax]` 内
2. **守恒精确**：剖面单元平均严格等于 U_i
3. **极值回退**：非单调模板（局部极值，即 `(U_i-U_{i-1})(U_{i+1}-U_i) ≤ 0`）时两端退回一阶值 U_i

**候选 3：THINC2**（仅两级模式使用）：同样的 THINC 公式，陡度参数换为 `bvd_beta2`。

### 1.3 择优流程

**第一级**（`select_function_1st` 对应）——"先默认、后覆盖"三步：

1. **全范围默认**：候选覆盖范围内**所有**单元先取高阶候选（对应原程序
   `us = poly` 的语义）；覆盖范围含择优范围两侧各一个**外延单元**——它们
   会被第二级 TBV 读取，必须一并初始化（见 §8.5）
2. **命中覆盖**：逐单元比较 TBV，`TBV_thinc < TBV_poly`（严格小于）的单元
   自身换成 THINC（平局取高阶候选）
3. **邻域传播**：命中单元的左右邻居（i-1、i+1）也换成 THINC（传播一层、
   不级联——flag 先整批固定，再统一改邻居）

**第二级**（`bvd_twostage=.true.` 时启用，`select_function_2st` 对应）：

1. 中间剖面与 THINC2(beta2) 再比一次 TBV（同样的跨面跳变定义）
2. 严格比较，THINC2 更小时覆盖
3. **无邻域传播**

---

## 2. 与 AMRVAC 的索引约定映射

这是理解代码的关键。AMRVAC 的 limiter 接口（`reconstruct_LR` 中调用）输出两个数组 `wLp`/`wRp`，它们**按面编号存储**，与"单元两端"的对应关系如下：

| AMRVAC 数组 | 含义 | 单元剖面视角 |
|-------------|------|-------------|
| `wLp(j)` | 面 j+1/2 的**左侧值** | 单元 **j 的 2 端（右端）** |
| `wRp(j)` | 面 j+1/2 的**右侧值** | 单元 **j+1 的 1 端（左端）** |

即：单元 j 的两端分别落在 `wLp(j)`（2 端）和 `wRp(j-1)`（1 端）两个位置。注意"单元 1 端 → 界面右侧、单元 2 端 → 界面左侧"的翻转。

与参考程序（标量版）的对应：

| 参考程序 | 含义 | AMRVAC |
|----------|------|--------|
| `poly_q(i,1)` | 单元 i 左端 | `wRp(i-1)` |
| `poly_q(i,2)` | 单元 i 右端 | `wLp(i)` |
| `fl(i)=poly_q(i-1,2)` | 面 i 左侧值 | `wLp(i-1)` |
| `fr(i)=poly_q(i,1)` | 面 i 右侧值 | `wRp(i-1)` |

### 2.1 BVDlimiter 内部的单元范围

```
候选剖面范围 ic : [iLmin-1, iLmax+2]   （候选两端值在这里计算）
择优范围     js : [iLmin,   iLmax+1]   （TBV 比较在这里进行）
传播目标范围     : js-1 与 js+1        （即 [iLmin-1, iLmax+2] 内部）
输出右值范围 iLp : [iLmin+1, iLmax+1]  （供 wRC 使用）
默认覆盖范围 iw3 : [iLmin-1, iLmax+1]  （uend2 默认初始化；WENO3 候选同范围）
默认覆盖范围 iw3p: [iLmin,   iLmax+2]  （uend1 默认初始化）
```

为什么 ic 要比 iL 左右各扩一个单元？因为 TBV 的跨面跳变定义需要**邻居单元的剖面端值**，而邻域传播还要**改写邻居的剖面**。第一级默认覆盖（iw3/iw3p）则必须覆盖**全部会被读取的单元**（含两侧外延单元），否则边界单元会读到未初始化值——排查与修复见 §8.5。

---

## 3. AMRVAC 宏系统（VAC 预处理器）详解

### 3.1 什么是 .t 文件

AMRVAC 源码文件扩展名为 `.t`（template）。它们**不是标准 Fortran**——在标准 Fortran 90 之上叠加了一套宏语法，由 Perl 脚本 `src/vacpp.pl` 在编译前展开成标准 Fortran（`.f` 文件）。

**为什么需要宏？** Fortran 90 没有模板/泛型机制。同一份求解器代码要覆盖 1D/2D/3D，还要物理模块无关。宏让开发者写一份代码，编译时按维度自动展开。

### 3.2 构建链条

```make
# arch/rules.make 中的规则
%.f: %.t
	$(VACPP) -d=$(NDIM) -v=$(NDIR) $< > $(@)
```

即：

```
.t 文件 --vacpp.pl 按维度展开--> .f 文件 --mpif90 编译--> .o --> libamrvac.a
```

`setup.pl -d=2` 就是通过 `NDIM=2` 让 vacpp 生成 2D 版本。展开产物在 `lib/<ndim>d<ndir>_<arch>/` 目录（如 `lib/1d1_default/mod_bvd.f`），**调试宏展开问题时直接看这个目录下的 .f 文件**。

### 3.3 宏速查表

| 宏 | 含义 | 1D 展开示例 |
|----|------|------------|
| `^D` | 维度索引 | `1` |
| `^ND` | 维度总数 | `1` |
| `^LIM` | 索引范围 `^D:^D` | `1:1` |
| `^L`（变量名后缀） | 声明/实参用的范围列表 | `ic^L` → `icmin1,icmax1` |
| `^S`（变量名后缀） | 数组下标用的范围 | `ic^S` → `icmin1:icmax1` |
| `^IFTWOD` / `^IFTHREED` | 仅 2D/3D 存在的代码块 | 1D 时整块删除 |
| `{...\}` | 宏语句块（`\` 是续行标记） | 展开为显式语句 |
| `^C` | 分量循环展开 | |
| `{^NOONED ...}` | 1D 下注释掉的块 | 1D 时删除 |

### 3.4 展开实例（来自 mod_bvd.t 真实代码）

**源码（.t）：**

```fortran
subroutine BVDlimiter(ixI^L,iL^L,idims,dxdim,w,wLC,wRC)
  integer          :: ic^L, icm^L, icp^L
  icmin^D=iLmin^D-kr(idims,^D); icmax^D=iLmax^D+2*kr(idims,^D);
  icm^L=ic^L-kr(idims,^D);
  icp^L=ic^L+kr(idims,^D);
```

**1D 展开（vacpp.pl -d=1，真实产物见 lib/1d1_default/mod_bvd.f）：**

```fortran
subroutine BVDlimiter(ixImin1,ixImax1,iLmin1,iLmax1,idims,dxdim,w,wLC,wRC)
  integer :: icmin1,icmax1, icmmin1,icmmax1, icpmin1,icpmax1
  icmin1=iLmin1-kr(idims,1); icmax1=iLmax1+2*kr(idims,1);
  icmmin1=icmin1-kr(idims,1);icmmax1=icmax1-kr(idims,1);
  icpmin1=icmin1+kr(idims,1);icpmax1=icmax1+kr(idims,1);
```

注意 vacpp 的一个细节：`icm^L=ic^L-kr(...)` 这样的语句展开后，**min 和 max 两条赋值用 `;` 连在同一行**（`icmmin1=...;icmmax1=...;`）。而源码里本来就带空格分号的（`icmin^D=...; icmax^D=...`）则保留空格分号分开。

**2D 时同一份代码会展开出（真实产物）：**

```fortran
subroutine BVDlimiter(ixImin1,ixImin2,ixImax1,ixImax2,iLmin1,iLmin2,iLmax1,&
   iLmax2,idims,dxdim,w,wLC,wRC)
  integer :: icmin1,icmin2,icmax1,icmax2, icmmin1,icmmin2,icmmax1,icmmax2, ...
  icmin1=iLmin1-kr(idims,1);icmin2=iLmin2-kr(idims,2)
  icmax1=iLmax1+2*kr(idims,1);icmax2=iLmax2+2*kr(idims,2);
```

注意：`^L` 的展开顺序是**先所有维度的 min，再所有维度的 max**（`icmin1,icmin2,icmax1,icmax2`），不是 min/max 交错。

### 3.5 各宏的展开规律

**`^L` 后缀（用于声明和实参列表）**：**先所有维度的 min，再所有维度的 max**

```fortran
! 源码                          ! 2D 展开
integer :: ic^L                  integer :: icmin1,icmin2,icmax1,icmax2
call foo(ixI^L, ic^L)            call foo(ixImin1,ixImin2,ixImax1,ixImax2, &
                                          icmin1,icmin2,icmax1,icmax2)
```

**`^S` 后缀（用于数组下标）**：每个维度的范围用冒号连接

```fortran
! 源码                          ! 2D 展开
w(ic^S)                          w(icmin1:icmax1,icmin2:icmax2)
```

**`^D`（维度索引）**：当前语句按维度重复（同维度多条语句以 `;` 连在同一行）

```fortran
! 源码                          ! 2D 展开
icmin^D=iLmin^D-kr(idims,^D);    icmin1=iLmin1-kr(idims,1);icmin2=iLmin2-kr(idims,2)
```

**`kr(idims,^D)`（单位偏移）**：`kr` 是 Kronecker delta 向量，`kr(idims,^D)` 在第 idims 维为 1、其余为 0。用它做索引偏移，保证只在目标方向上平移一格。

### 3.6 重要注意事项

1. **行内中文注释会被破坏**：vacpp 对代码行做"引号掩码"处理（把引号内字符转到高 ASCII 区间），行尾 `!` 之后的中文 UTF-8 字节会被破坏成乱码。**中文注释必须独立成行**（行首 `!`），独立注释行是原样透传的。
2. **代码行长度限制**：vacpp 对超长且无法断行的语句会报 `Couldn't break line`。长字符串（如 `mpistop` 消息）要保持短。
3. **`&` 续行**：数组表达式写太长时用 Fortran 标准 `&` 续行，vacpp 会自动处理。
4. **`amrvac.h`**：`setup.pl` 创建的空头文件，供 gfortran 的 `-cpp` 标准预处理器使用——与 vacpp 是**两套东西**，不要混淆。

---

## 4. 代码走读：mod_bvd.t

完整源码位于 `src/limiter/mod_bvd.t`（306 行）。结构如下：

```
module mod_bvd
  ├── 参数：bvd_eps = 1.0d-20（THINC 归一化小量）
  ├── 公开接口：BVDlimiter
  ├── 助手子程序：bvd_thinc（完整 THINC 剖面）
  └── 主入口：BVDlimiter（候选计算 → 择优 → 传播 → 写回）
```

### 4.1 助手子程序 bvd_thinc

对应参考程序的 `thinc(q,thinc_q,bt)`，逐位移植。输入单元范围 `ic`、陡度 `beta`、单元中心值 `w`，输出剖面两端 `tend1/tend2`（1 端=左端，2 端=右端）。

```fortran
subroutine bvd_thinc(ixI^L,ic^L,idims,beta,w,tend1,tend2)
  ...
  icm^L=ic^L-kr(idims,^D);   ! 左邻居 U_{i-1}
  icp^L=ic^L+kr(idims,^D);   ! 右邻居 U_{i+1}

  do iw=1,nwflux
     ! 允许区间：以左右邻居为界
     uminC(ic^S) = min(w(icm^S,iw), w(icp^S,iw))
     urangeC(ic^S) = max(w(icm^S,iw), w(icp^S,iw)) - uminC(ic^S)
     ! 方向符号
     rsignC(ic^S) = sign(one, w(icp^S,iw)-w(icm^S,iw))
     ! 归一化单元平均
     ccC(ic^S) = (w(ic^S,iw)-uminC(ic^S)+bvd_eps)/(urangeC(ic^S)+bvd_eps)
     ! 守恒约束解出移位参数 A
     bigB(ic^S) = exp(rsignC(ic^S)*beta*(2.0d0*ccC(ic^S)-1.0d0))
     bigA(ic^S) = (bigB(ic^S)/cbeta-1.0d0)/tbeta
     ! 两端值（1 端 = 左端，2 端 = 右端）
     tend1(ic^S,iw) = uminC(ic^S)+half*urangeC(ic^S)*(1.0d0+rsignC(ic^S)*bigA(ic^S))
     tend2(ic^S,iw) = uminC(ic^S)+half*urangeC(ic^S) &
          *(1.0d0+rsignC(ic^S)*(tbeta+bigA(ic^S))/(1.0d0+bigA(ic^S)*tbeta))
     ! 局部极值 → 退回一阶
     where ((w(ic^S,iw)-w(icm^S,iw))*(w(icp^S,iw)-w(ic^S,iw)) <= 0.0d0)
        tend1(ic^S,iw) = w(ic^S,iw)
        tend2(ic^S,iw) = w(ic^S,iw)
     end where
  end do
end subroutine bvd_thinc
```

要点：

- 所有临时数组（`uminC/urangeC/ccC/...`）都是 `ixI^S` 大小，只在 `ic^S` 范围内赋值——这是 AMRVAC 的惯例（大数组一次分配，逐调用复用范围）
- `do iw=1,nwflux` 循环逐变量处理（原始变量逐分量择优）

### 4.2 主入口 BVDlimiter

流程分五步：

**第 1 步：定义单元范围**

```fortran
icmin^D=iLmin^D-kr(idims,^D); icmax^D=iLmax^D+2*kr(idims,^D);  ! 候选范围
jsmin^D=iLmin^D; jsmax^D=iLmax^D+kr(idims,^D);                  ! 择优范围
icm^L=ic^L-kr(idims,^D);  icp^L=ic^L+kr(idims,^D);
jsm^L=js^L-kr(idims,^D);  jsp^L=js^L+kr(idims,^D);
iLp^L=iL^L+kr(idims,^D);
iw3min^D=iLmin^D-kr(idims,^D); iw3max^D=iLmax^D+kr(idims,^D);  ! uend2 默认初始化范围
iw3p^L=iw3^L+kr(idims,^D);                                      ! uend1 默认初始化范围
```

**第 2 步：计算候选剖面两端**

```fortran
! 二次多项式（对应 second_degree_polynomial）
pend2(ic^S,1:nwflux) = (-1.0d0*w(icm^S,1:nwflux)+5.0d0*w(ic^S,1:nwflux) &
     +2.0d0*w(icp^S,1:nwflux))/6.0d0
pend1(ic^S,1:nwflux) = (2.0d0*w(icm^S,1:nwflux)+5.0d0*w(ic^S,1:nwflux) &
     -1.0d0*w(icp^S,1:nwflux))/6.0d0
! 完整 THINC（对应 thinc）
call bvd_thinc(ixI^L,ic^L,idims,thinc_beta1,w,tend1,tend2)
```

**第 3 步：第一级择优（全范围默认 → 命中覆盖 → 邻域传播）**

```fortran
! ① 全范围默认：所有会被读取的单元先取高阶候选（对应原程序 us = poly）
uend2(iw3^S,1:nwflux)=pend2(iw3^S,1:nwflux)
uend1(iw3p^S,1:nwflux)=pend1(iw3p^S,1:nwflux)
! TBV：单元两个界面的跨面跳变之和
tbvP(js^S,1:nwflux) = abs(pend2(js^S,1:nwflux)-pend1(jsp^S,1:nwflux)) &
     +abs(pend2(jsm^S,1:nwflux)-pend1(js^S,1:nwflux))
tbvT(js^S,1:nwflux) = abs(tend2(js^S,1:nwflux)-tend1(jsp^S,1:nwflux)) &
     +abs(tend2(jsm^S,1:nwflux)-tend1(js^S,1:nwflux))

do iw=1,nwflux
   ! ② 命中覆盖：严格小于才换 THINC（未命中保持默认高阶候选）
   flagj(js^S)=.false.
   where (tbvT(js^S,iw) < tbvP(js^S,iw))
      uend2(js^S,iw)=tend2(js^S,iw)
      uend1(js^S,iw)=tend1(js^S,iw)
      flagj(js^S)=.true.
   end where
   ! ③ 邻域传播：命中单元的左右邻居（i-1、i+1）也换 THINC（一层，不级联）
   where (flagj(js^S))
      uend2(jsm^S,iw)=tend2(jsm^S,iw)
      uend1(jsm^S,iw)=tend1(jsm^S,iw)
      uend2(jsp^S,iw)=tend2(jsp^S,iw)
      uend1(jsp^S,iw)=tend1(jsp^S,iw)
   end where
end do
```

**第 4 步：第二级择优（可选，无传播）**

```fortran
if (bvd_twostage) then
   call bvd_thinc(ixI^L,ic^L,idims,thinc_beta2,w,t2end1,t2end2)
   tbvU(js^S,1:nwflux) = abs(uend2(js^S,1:nwflux)-uend1(jsp^S,1:nwflux)) &
        +abs(uend2(jsm^S,1:nwflux)-uend1(js^S,1:nwflux))
   tbvT2(js^S,1:nwflux) = abs(t2end2(js^S,1:nwflux)-t2end1(jsp^S,1:nwflux)) &
        +abs(t2end2(jsm^S,1:nwflux)-t2end1(js^S,1:nwflux))
   do iw=1,nwflux
      where (tbvT2(js^S,iw) < tbvU(js^S,iw))
         uend2(js^S,iw)=t2end2(js^S,iw)
         uend1(js^S,iw)=t2end1(js^S,iw)
      end where
   end do
end if
```

**第 5 步：写回界面值 + 物理合法性回退**

```fortran
do iw=1,nwflux
   wLCtmp(iL^S,iw)=uend2(iL^S,iw)   ! 单元 j 2 端（右端）→ 面 j+1/2 左侧
   wRCtmp(iL^S,iw)=uend1(iLp^S,iw)  ! 单元 j+1 1 端（左端）→ 面 j+1/2 右侧
end do
call fix_limiter(ixI^L,iL^L,wLCtmp,wRCtmp,wLC,wRC)
```

`fix_limiter`（来自 `mod_weno`）是安全网：对重构出的界面状态逐面检查物理合法性（密度/压力为正等），非物理处保留调用前的一阶迎风值。

---

## 5. 集成点：所有修改的文件

| 文件 | 改动 | 作用 |
|------|------|------|
| `src/limiter/mod_bvd.t` | **新建** | BVD 核心算法 |
| `src/limiter/mod_limiter.t` | `limiter_bvd=26` 常量、`'bvd'` 解析、`bvd_beta/bvd_beta2/bvd_twostage/bvd_highorder/bvd_selection` 参数、dwlimiter2 报错分支 | 注册限制器 |
| `src/mod_finite_volume.t` | `reconstruct_LR` 加 `use mod_bvd` 和 `case (limiter_bvd)` | 分发入口 |
| `src/io/mod_input_output.t` | namelist `bvd_beta,bvd_beta2,bvd_twostage,bvd_highorder,bvd_selection`、默认值、`nghostcells≥3` | 参数解析与 ghost cell |
| `src/limiter/makefile` | `FOBJECTS += mod_bvd.t` | 编译注册 |
| `src/makefile` | `mod_bvd.o` 依赖、`mod_finite_volume.o: mod_bvd.mod` | 依赖排序 |
| `tests/schemes/3step_hll_bvd.par` | **新建** | 测试用 par 文件 |
| `tests/rho/auto_1d/test.make` | SCHEMES 加 `3step_hll_bvd` | 测试注册 |

### 5.1 限制器分发机制（reconstruct_LR）

`src/mod_finite_volume.t` 的 `reconstruct_LR` 是所有特殊限制器的统一入口，按 `type_limiter(block%level)` 分发：

```fortran
select case (type_limiter(block%level))
case (limiter_mp5)
   call MP5limiter(ixI^L,ixL^L,idims,w,wLp,wRp)
case (limiter_weno5)
   call WENO5limiter(ixI^L,ixL^L,idims,dxdim,w,wLp,wRp,1)
...
case (limiter_bvd)
   call BVDlimiter(ixI^L,ixL^L,idims,dxdim,w,wLp,wRp)
...
end select
```

BVD 与 MP5/WENO 系列一样，是"一次性重构全部变量"的特殊限制器——这就是为什么 BVD 的单元剖面逻辑能完整放进 limiter 内部。

### 5.2 ghost cell 要求

`mod_input_output.t` 中：

```fortran
if(any(limiter(1:nlevelshi)=='bvd')) then
  nghostcells=max(nghostcells,3)
end if
```

BVD 的模板（二次多项式 3 点 + 邻域 TBV 各扩一层）合计需要单元 i-2..i+3，与 WENO5 相同，3 层 ghost cell 即可。

---

## 6. 使用方法

在 `amrvac.par` 的 `&methodlist` 中：

```fortran
 &methodlist
    time_stepper = 'threestep'
    flux_scheme  = 20*'hll'
    limiter      = 20*'bvd'
    bvd_beta     = 1.05d0     ! 一级 THINC 陡度（可选，默认 1.05）
    bvd_beta2    = 1.2d0      ! 二级 THINC 陡度（可选，默认 1.2）
    bvd_twostage = .true.     ! 开启两阶段 BVD（可选，默认 .false.）
    bvd_highorder = 2         ! 高阶候选：1=二次多项式（默认），2=WENO3
    bvd_selection = .false.   ! 关闭选择机制，只用高阶候选（诊断用，默认 .true.）
 /
```

参数说明：

| 参数 | 类型 | 默认 | 含义 |
|------|------|------|------|
| `bvd_beta` | real | 1.05 | 第一级 THINC 陡度，越大越锐利（建议 0.5~3.5） |
| `bvd_beta2` | real | 1.2 | 第二级 THINC 陡度（仅 twostage 时用） |
| `bvd_twostage` | logical | .false. | 是否启用两阶段 BVD |
| `bvd_highorder` | integer | 1 | 高阶候选选择：**1=二次多项式**（与原程序 `second_degree_polynomial` 逐位一致）；**2=WENO3**（调用框架 `WENO3limiter`，与 `limiter='weno3'` 同款重构） |
| `bvd_selection` | logical | .true. | 是否启用 BVD 选择机制（THINC 候选 + TBV 择优 + 邻域传播）。设为 `.false.` 时跳过全部选择逻辑，直接采用高阶候选两端值——**诊断/对照专用** |

> **诊断开关**：`bvd_selection=.false.` 用于隔离问题——例如判断"太阳风 NaN 是否来自选择机制"。此时 THINC 与 TBV 完全不参与，仅高阶候选（二次多项式或 WENO3）经写回与 `fix_limiter` 后输出。

**兼容性**：

- 支持 1D/2D/3D
- 兼容所有走 `reconstruct_LR` 路径的通量格式（hll, hllc, hlld, tvdlf, tvdmu, cd, cd4 等）
- 用 `tvd` 通量格式时会得到明确的报错提示

---

## 7. 编译与测试

### 7.1 编译

```bash
export AMRVAC_DIR=/path/to/amrvac-BVD
cd tests/rho/auto_1d          # 或任意测试目录
$AMRVAC_DIR/setup.pl -d=1 -arch=default
make -f test.make ball_1d_3step_hll_bvd.log
```

### 7.2 测试结构

- `tests/schemes/3step_hll_bvd.par`：scheme 片段（时间步进器+通量+限制器），被 `tests/rho/auto_1d/ball_1d.par` 拼接
- 测试 = 1D 间断球（ball）在周期域上的平流，校验 `mean(w)`（守恒）与 `mean(w**2)`（耗散）序列
- 回归参考在 `tests/rho/auto_1d/correct_output/`

### 7.3 调试技巧

- 看宏展开结果：`lib/1d1_default/mod_bvd.f`
- 单独预处理器检查：`src/vacpp.pl -d=1 -v=1 src/limiter/mod_bvd.t > /tmp/test.f`
- 改算法后需**重新生成参考**：跑一遍测试 → 拷贝 log 到 `correct_output/` → 再跑确认 PASSED

---

## 8. 验证结果

### 8.1 1D ball 平流测试（`tests/rho/auto_1d`）

3 步 SSPRK + HLL，courant=0.9：

| 指标 | 结果 |
|------|------|
| vacpp 预处理（1D/2D/3D） | 全部通过 |
| 编译链接 | 通过 |
| 正向/反向回归测试 | 均 PASSED |
| 质量守恒 | mean(w) 保持 1.10000001（浮点舍入级） |
| 单调性 | mean(w²) 单调衰减，无 NaN/Inf |

### 8.2 HD Sod 激波管测试（`tests/hd/Riemann_1D`）

修正型 Sod 问题：左态 (ρ=1.5, p=1.0)，右态 (ρ=0.5, p=1/3)，γ=5/3，
3 步 SSPRK + HLLC，固定均匀网格 64 单元，t=0.2（经典波系时刻）。
对照 PPM 与 minmod：

| 指标 | BVD | PPM | minmod |
|------|-----|-----|--------|
| 稳定性 | 无 NaN，运行稳定 | 同 | 同 |
| rho 范围 | **[0.50000, 1.50000]** | [0.50023, 1.49990] | [0.50000, 1.50000] |
| 接触间断过渡 | **≈2 单元** | ≈6 单元 | ≈4 单元 |
| 激波过渡 | **≈2 单元** | ≈8 单元 | ≈3 单元 |
| 质量守恒 mean(ρ) | 0.99961 | 0.999998 | 1.000000 |
| 接触区小幅抖动 | ±0.006 | 平滑 | 平滑 |

结论：

1. **无振荡**：BVD 的密度严格落在初值界 [0.5, 1.5] 内，无 Gibbs 型过冲/下冲
2. **间断最锐利**：接触间断和激波均 ~2 单元过渡，明显优于 PPM（6/8 单元）
   和 minmod（4/3 单元）——THINC 候选在非线性问题中发挥作用
3. **代价**：接触区有 ±0.6% 的小幅抖动（THINC 锐利的典型代价），
   PPM/minmod 更平滑但更弥散
4. 质量守恒到 0.04%（minmod 为精确），总能量守恒（内能→动能转换正确）

### 8.3 光滑平流收敛阶测试（`tests/rho/analytic_test_1d`）

sin 剖面平流（iprob=1，有解析解），3 步 SSPRK + HLL，t=1.0，
`domain_nx1 = 16/32/64/128/256` 五档分辨率，L1 误差与收敛阶：

| 限制器 | L1 误差（N=256） | 各段收敛阶 | 平均 |
|--------|-----------------|-----------|------|
| **BVD** | **5.23e-6** | 2.96 / 2.99 / 3.00 / 3.00 | **≈3.0** |
| PPM | 1.29e-4 | 2.08 / 2.11 / 2.33 / 2.49 | ≈2.3 |
| minmod | 1.56e-3 | 1.35 / 1.60 / 1.78 / 1.87 | ≈1.7 |

结论：

1. **BVD 收敛阶 = 3.0**，与二次多项式候选的理论阶一致——证明光滑区
   BV 判据正确选择了多项式（若 THINC 被误选，阶数会退化）
2. BVD 在 N=256 时比 PPM 精确 **25 倍**、比 minmod 精确 **300 倍**
3. PPM/minmod 阶数低于理论值的原因：二者在极值处限幅过强；BVD 的
   二次多项式在极值处因 BV 判据胜出而**不限制**，保持全 3 阶
4. 复现命令：`cd tests/rho/analytic_test_1d && ./conv_test.sh`

### 8.4 MHD 强背景场测试：NaN 根因与修复（已解决）

**现象**：2.5D 太阳风算例（`tests/mhd/solar_wind_2.5D/xia2018.par`，
Powell divB 修正 + B0 背景场，内边界附近 β~1e-5）中，BVD 在前 ~10 步
（t<0.01）内即出现 NaN（8 个守恒量体积平均全部 NaN）。
`fix_small_values=.true.` 无改善。NaN 进入态后，CFL 中 gfortran 的
min/max 会忽略 NaN 操作数，dt 失去约束、反增至 dtmax，程序携带
NaN 空转（日志持续刷 NaN）。

**排查过程**（逐步缩小范围）：

| 实验 | 结果 | 排除项 |
|------|------|--------|
| koren / weno3（仅换限制器） | 0 NaN | 代码环境、par 配置 |
| BVD 换 WENO3 候选（`bvd_highorder=2`） | 同一步 NaN（t≈0.01） | 高阶候选的选择 |
| BVD 关选择机制（`bvd_selection=.false.`） | 同一步 NaN | THINC 候选 + TBV 择优 + 邻域传播 |
| THINC 压平（`bvd_beta=1e-10`） | 同一步 NaN | THINC 陡度 |
| **weno5（框架自带、与 BVD 零共用代码）** | **同一步 NaN** | **BVD 代码本身** |
| 径向边界改 `'cont'`（框架自动填满 ghost） | weno5 与 BVD 均 0 NaN | **证实：ghost 层填充不足** |

**根因**：测试算例的 `specialbound_usr`（case 1，内边界）按"2 层
ghost"硬编码——只写了 `w(ixOmax1,...)` 与 `w(ixOmin1,...)` 两个赋值。
koren 模板宽 1 层、weno3 宽 2 层，恰好够用；而 weno5 与 BVD 的模板
宽 3 层，会把 `nghostcells` 抬到 3——于是**中间那层 ghost 从未被
填充**（保持初值 0，即 ρ=0），第一步重构（含 THINC/多项式读取）
即产生 NaN。

**修复**：把 `specialbound_usr` case 1 的 ghost 填充推广到任意层数
——对第 k 层（k=1 为贴内边界层）用径向线性外推

```
q(k) = 2k·q_bnd − (2k−1)·q(I1)
```

k=1,2 与原公式逐位一致（`nghostcells=2` 时行为完全不变），k=3
及以后自然扩展。实现见 `tests/mhd/solar_wind_2.5D/mod_usr.t`
（原版备份 `mod_usr.t.bak_ghostfix`）。

**修复后验证**（三组同 par、t=0.1 短跑，全部穿过原崩溃时刻）：

| 配置 | 修复前 | 修复后 |
|------|--------|--------|
| weno5 | t≈0.01 NaN | **0 NaN** ✅ |
| BVD（纯二次多项式，`bvd_selection=.false.`） | t≈0.01 NaN | **0 NaN** ✅ |
| BVD（完整：THINC + 选择 + 二次多项式） | t≈0.01 NaN | **0 NaN** ✅ |

**长时程验证（t=10.0 全时长）**：单级 BVD 完整跑完，**0 NaN**。同 par
仅换限制器的对照（各时刻体积平均量的最大相对差，最大出现在磁分量 b2）：

- BVD vs weno3：t=1~10 各时刻约 1.6e-3 ~ 4.5e-3
- BVD vs koren（其数据到 t≈0.98）：约 2.6e-3 ~ 1.4e-2

**尚未验证**：

- 与原程序的整体逐位数值对照（第一级选择逻辑的语义等价性已单独验证，见 §8.5）

### 8.5 第一级选择逻辑的语义等价性验证

复核指出第一级选择的正确语义是"**先默认、后覆盖**"：

1. 先给**所有**单元（含择优范围两侧各一个外延单元）选高阶候选
   （默认二次多项式，对应原程序 `us = poly`）；
2. 单元 i 命中（`TBV_thinc < TBV_poly`）时，把 **i-1、i、i+1** 三个单元
   的备选函数都换成 THINC。

`mod_bvd.t` 的第一级已按此实现为三步（§4 第 3 步）：① 全范围默认
（`uend2` 在 iw3、`uend1` 在 iw3p 预置高阶候选）→ ② 命中单元自身换
THINC → ③ 命中单元的左右邻居也换 THINC。

**差异与修复**：修复前的写法为"逐单元 `where/elsewhere` 二选一 + 条件
传播"——默认值只覆盖择优范围，两侧外延单元（`uend2(iLmin-1)`、`uend1(iLmax+2)`）仅在邻居命中时才被写入；
邻居未命中时它们保持**未初始化值**。外延单元只在两级模式（第二级 TBV）
中被读取，因此该缺陷对单级模式不可见，对两级模式会引入未定义行为。

**穷举对照**：独立程序 `bvd_select_equiv_test.f90`（仓库根目录）对
32 种命中模式（5 个内部单元命中情况的全组合）逐一比对：

| 写法 | 择优范围内单元 | 两侧外延单元（二级读取） |
|------|----------------|--------------------------|
| 修复前（逐单元二选一 + 条件传播） | 0 处差异 | **32 处差异**（未初始化） |
| 现行（默认 + 命中 + 传播） | 0 处差异 | **0 处差异**（全部被读取单元 0 差异） |

复现：`gfortran -O0 -o /tmp/bvd_equiv bvd_select_equiv_test.f90 && /tmp/bvd_equiv`

**回归**（重构后）：1D 平流正/反测试 PASSED；2D 平流正/反测试 PASSED；
1D 与太阳风两级冒烟 0 NaN；单级模式与修复前逐位一致（太阳风工况
t=0.01~0.1 对比相对差 0.000e+00），故此前 t=10 长时程结论不受影响。

---

## 9. 二次开发指南

### 9.1 替换高阶候选（二次多项式 / WENO3）

高阶候选默认是二次多项式；运行时置 `bvd_highorder=2` 可切换为 WENO3
（见 §6）。若要换成其它公式，修改 `BVDlimiter` 中 `pend2/pend1` 的
计算即可。例如换回 5 阶多项式：

```fortran
pend2(ic^S,1:nwflux) = (2.0d0*w(icmm^S)-13.0d0*w(icm^S)+47.0d0*w(ic^S) &
     +27.0d0*w(icp^S)-3.0d0*w(icpp^S))/60.0d0
```

注意：5 点模板需要 `nghostcells≥4`（`mod_input_output.t` 同步修改）。

### 9.2 调 THINC 锐利度

运行时 namelist 调整 `bvd_beta`，无需重编译。

### 9.3 换 TBV 判据

修改 `tbvP/tbvT` 的定义。当前是跨面跳变；若改为相对单元均值：

```fortran
tbvP(js^S,1:nwflux) = abs(pend2(js^S)-w(js^S)) + abs(pend1(js^S)-w(js^S))
```

### 9.4 加特征投影（提升 MHD/HD 鲁棒性）

当前择优在原始变量上逐分量进行。对强波系统，可在候选计算前对 `w` 做特征分解、在特征空间择优、再投影回原始变量——参考 `mod_weno.t` 中 WENO 系列对特征分解的处理方式。

### 9.5 换新限制器名

若要并存多个 BVD 变体：

1. `mod_limiter.t` 加新常量（如 `limiter_bvd2=27`）和 `limiter_type` 解析分支
2. `reconstruct_LR` 加新 `case`
3. `mod_input_output.t` 加 ghost cell 分支（如需要）

---

## 参考文献

- Sun, Teng, Xiao, *J. Comput. Phys.* 320 (2016) 414-435 —— BVD 框架
- Deng, Xie, Loubère, Shimizu, Xiao, *J. Comput. Phys.* (2019) —— 两阶段 BVD
- Xiao, Ii, Chen, *Int. J. Numer. Meth. Fluids* 47 (2005) —— THINC 原始文献
