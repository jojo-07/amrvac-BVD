!> BVD（Boundary Variation Diminishing，边界变差递减）限制器
!> —— 单元剖面择优版，判据与原程序 select_function_1st/2st 一致
!>
!> 命名约定（单元端用 1/2，避免与界面左/右混淆）：
!>   1 端 = 单元左端 = 面 i-1/2 的右侧值（对应原程序 poly_q(i,1)）
!>   2 端 = 单元右端 = 面 i+1/2 的左侧值（对应原程序 poly_q(i,2)）
!>   注意：与 AMRVAC 界面侧命名（wLp/wRp）视角翻转——
!>   单元的"1 端"落在界面的"右"侧、单元的"2 端"落在界面的"左"侧。
!>
!> 对每个单元 i 构造完整剖面（1 端/2 端），按单元总边界变差 TBV 选剖面：
!>
!>     TBV(q, i) = |q_2(i) - q_1(i+1)| + |q_2(i-1) - q_1(i)|
!>
!> 即单元 i 两个界面处的跨面跳变之和（经典 BVD 定义）。
!> 严格比较：TBV_thinc < TBV_hi 时取 THINC，平局取高阶候选。
!>
!> 候选剖面：
!>   (1) 高阶候选（由 bvd_highorder 选择）：
!>       =1（默认）：二次多项式，3 点二次插值
!>         （与原程序 second_degree_polynomial 逐位一致）
!>           2 端 = (-1 U_{i-1} + 5 U_i + 2 U_{i+1})/6
!>           1 端 = ( 2 U_{i-1} + 5 U_i - 1 U_{i+1})/6
!>       =2：WENO3 重构（调用框架 WENO3limiter，与 limiter='weno3' 一致）
!>   (2) THINC：双曲正切界面捕捉重构（与原程序 thinc 逐位一致）
!>       剖面 = umin + 0.5*umax*(1 + r*tanh(βξ + atanh(A)))，
!>       其中 umin=min(U_{i-1},U_{i+1})、umax=|U_{i+1}-U_{i-1}|、
!>       r=sign(U_{i+1}-U_{i-1})、cc=(U_i-umin)/umax，
!>       A 由单元平均守恒条件解出：
!>         B = exp(r*β*(2*cc-1))，A = (B/cosh(β)-1)/tanh(β)
!>       1 端 = umin + 0.5*umax*(1 + r*A)
!>       2 端 = umin + 0.5*umax*(1 + r*(tanh(β)+A)/(1+A*tanh(β)))
!>       非单调模板（局部极值）时两端退回一阶值 U_i；
!>       剖面两端值恒落在 [umin, umin+umax] 内（保单调）。
!>   (3) THINC2：第二级 THINC，陡度 bvd_beta2（默认 1.2）
!>
!> 两级择优（bvd_twostage=.true.）：
!>   第一级（先默认、后覆盖，对应原程序 select_function_1st）：
!>           先给所有会被读取的单元默认取高阶候选（us = poly 语义），
!>           再把命中单元（TBV_thinc < TBV_hi，严格小于）及其左右邻居
!>           （i-1、i+1）换成 THINC（传播一层，不级联）
!>   第二级：中间剖面与 THINC(beta2) 按 TBV 择优（无邻域传播）
!>
!> 诊断开关（bvd_selection=.false.）：跳过 THINC 与两级 TBV 择优，
!> 直接采用高阶候选剖面（用于对照试验与问题定位）。
!>
!> 模板宽度：候选与邻域 TBV 合计需要单元 i-2..i+3 的原始值，
!> 与 WENO5 相同，故 3 层 ghost cell 即可（mod_input_output.t 已设置）。
!> 接口与 WENO5limiter 完全一致，在 reconstruct_LR 中按同模式调用。
!> 择优在原始变量上逐分量进行（含传播），对强波系统可加特征投影。
!>
!> 参考文献：
!>   Sun, Teng, Xiao, J. Comput. Phys. 320 (2016) 414-435 （BVD 框架）
!>   Deng, Xie, Loubere, Shimizu, Xiao, J. Comput. Phys. (2019)
!>   Xiao, Ii, Chen, Int. J. Numer. Meth. Fluids 47 (2005)（THINC 原始文献）
module mod_bvd
  ! 全局参数：nw, nwflux, kr, zero, half, one 等
  use mod_global_parameters
  ! BVD 控制参数（/methodlist/ namelist 可调）
  use mod_limiter, only: bvd_beta, bvd_beta2, bvd_twostage, bvd_highorder, bvd_selection
  ! 非物理解时回退一阶的安全网；WENO3 高阶候选
  use mod_weno, only: fix_limiter, WENO3limiter

  implicit none
  private

  ! THINC 归一化小量（与原程序 e=1.d-20 一致，防止 0/0）
  double precision, parameter :: bvd_eps = 1.0d-20

  ! 本模块唯一对外接口
  public :: BVDlimiter

contains

  !> 完整 THINC 剖面两端值（与原程序 thinc(q,thinc_q,bt) 逐位一致）
  !>
  !> 输入：ixI 全域范围；ic 候选单元范围（idims 方向按单元编号）；
  !>      idims 当前维度；beta 陡度参数（必须 >0）；w 单元中心原始变量
  !> 输出：tend1/tend2 单元剖面 1 端（左端）/2 端（右端）值
  subroutine bvd_thinc(ixI^L,ic^L,idims,beta,w,tend1,tend2)
    integer, intent(in)             :: ixI^L, ic^L, idims
    double precision, intent(in)    :: beta
    double precision, intent(in)    :: w(ixI^S,1:nw)
    double precision, intent(inout) :: tend1(ixI^S,1:nw), tend2(ixI^S,1:nw)

    !> 邻居最小值 umin 与极差 umax
    double precision :: uminC(ixI^S), urangeC(ixI^S)
    !> 归一化单元平均 cc、方向符号 r、中间量 B 与 A
    double precision :: ccC(ixI^S), rsignC(ixI^S), bigB(ixI^S), bigA(ixI^S)
    !> tanh(β) 与 cosh(β)
    double precision :: tbeta, cbeta
    !> 单元 ic 的左右邻居索引
    integer          :: icm^L, icp^L
    !> 变量循环指标
    integer          :: iw

    icm^L=ic^L-kr(idims,^D);
    icp^L=ic^L+kr(idims,^D);
    tbeta = tanh(beta)
    cbeta = cosh(beta)

    do iw=1,nwflux
       ! 以左右邻居定义剖面值的允许区间
       uminC(ic^S) = min(w(icm^S,iw), w(icp^S,iw))
       urangeC(ic^S) = max(w(icm^S,iw), w(icp^S,iw)) - uminC(ic^S)
       ! 方向符号：r = sign(U_{i+1} - U_{i-1})
       rsignC(ic^S) = sign(one, w(icp^S,iw)-w(icm^S,iw))
       ! 归一化单元平均（加小量防 0/0）
       ccC(ic^S) = (w(ic^S,iw)-uminC(ic^S)+bvd_eps) &
            /(urangeC(ic^S)+bvd_eps)
       ! 由单元平均守恒条件解出移位参数 A
       bigB(ic^S) = exp(rsignC(ic^S)*beta*(2.0d0*ccC(ic^S)-1.0d0))
       bigA(ic^S) = (bigB(ic^S)/cbeta-1.0d0)/tbeta
       ! 剖面两端值（1 端 = 左端，2 端 = 右端）
       tend1(ic^S,iw) = uminC(ic^S)+half*urangeC(ic^S) &
            *(1.0d0+rsignC(ic^S)*bigA(ic^S))
       tend2(ic^S,iw) = uminC(ic^S)+half*urangeC(ic^S) &
            *(1.0d0+rsignC(ic^S)*(tbeta+bigA(ic^S)) &
            /(1.0d0+bigA(ic^S)*tbeta))
       ! 非单调模板（局部极值）：两端退回一阶值
       where ((w(ic^S,iw)-w(icm^S,iw))*(w(icp^S,iw)-w(ic^S,iw)) <= 0.0d0)
          tend1(ic^S,iw) = w(ic^S,iw)
          tend2(ic^S,iw) = w(ic^S,iw)
       end where
    end do

  end subroutine bvd_thinc

  !> BVD 限制器主入口（单元剖面择优版）
  !>
  !> 输入：
  !>   ixI^L  全域（含 ghost cell）范围
  !>   iL^L   需要重构的界面范围（idims 方向按单元编号）
  !>   idims  当前维度方向
    !>   dxdim  该方向网格间距（bvd_highorder=2 时供 WENO3 候选使用）
  !>   w      单元中心的原始变量
  !> 输出（界面侧命名，注意与单元端 1/2 视角翻转）：
  !>   wLC    界面左侧重构值；wLC(j) = 单元 j 的 2 端（右端）
  !>   wRC    界面右侧重构值；wRC(j) = 单元 j+1 的 1 端（左端）
  subroutine BVDlimiter(ixI^L,iL^L,idims,dxdim,w,wLC,wRC)
    integer, intent(in)             :: ixI^L, iL^L, idims
    double precision, intent(in)    :: dxdim
    double precision, intent(in)    :: w(ixI^S,1:nw)
    double precision, intent(inout) :: wRC(ixI^S,1:nw), wLC(ixI^S,1:nw)

    !> 高阶候选剖面两端值（bvd_highorder=1 二次多项式 / =2 WENO3；1 端/2 端）
    double precision :: pend1(ixI^S,1:nw), pend2(ixI^S,1:nw)
    !> 第一级 THINC 剖面两端值（1 端/2 端）
    double precision :: tend1(ixI^S,1:nw), tend2(ixI^S,1:nw)
    !> 第二级 THINC2 剖面两端值（1 端/2 端，仅两级择优时使用）
    double precision :: t2end1(ixI^S,1:nw), t2end2(ixI^S,1:nw)
    !> 被选剖面两端值（1 端/2 端；第一级结果含邻域传播，两级模式再经第二级）
    double precision :: uend1(ixI^S,1:nw), uend2(ixI^S,1:nw)
    !> 各候选剖面的单元总边界变差 TBV
    double precision :: tbvP(ixI^S,1:nw), tbvT(ixI^S,1:nw)
    double precision :: tbvU(ixI^S,1:nw), tbvT2(ixI^S,1:nw)
    !> 第一级选中 THINC 的标记（用于邻域传播）
    logical          :: flagj(ixI^S)
    !> 最终界面值（物理合法性回退之前的暂存）
    double precision :: wLCtmp(ixI^S,1:nw), wRCtmp(ixI^S,1:nw)
    !> WENO3 高阶候选的界面左右值（仅 bvd_highorder=2 时使用）
    double precision :: wLC3(ixI^S,1:nw), wRC3(ixI^S,1:nw)
    !> THINC 陡度参数（小量下限防除零，β→0 极限即一阶迎风）
    double precision :: thinc_beta1, thinc_beta2
    !> 候选剖面单元范围 ic = [iLmin-1, iLmax+2] 及其偏移邻居
    integer          :: ic^L, icm^L, icp^L
    !> 择优单元范围 js = [iLmin, iLmax+1] 及其左右邻居
    integer          :: js^L, jsm^L, jsp^L
    !> 输出右值所需单元范围 iLp = [iLmin+1, iLmax+1]
    integer          :: iLp^L
    !> 覆盖范围 iw3 = [iLmin-1, iLmax+1]（WENO3 调用 / uend2 默认初始化）；
    !> iw3p = iw3+1 = [iLmin, iLmax+2]（uend1 默认初始化）
    integer          :: iw3^L, iw3p^L
    !> 变量循环指标
    integer          :: iw

    ! 候选剖面范围：邻域 TBV 需要左右各多一个单元
    icmin^D=iLmin^D-kr(idims,^D); icmax^D=iLmax^D+2*kr(idims,^D);
    ! 择优范围：输出所需单元 iLmin..iLmax+1
    jsmin^D=iLmin^D; jsmax^D=iLmax^D+kr(idims,^D);
    ! 各范围的偏移邻居
    icm^L=ic^L-kr(idims,^D);
    icp^L=ic^L+kr(idims,^D);
    jsm^L=js^L-kr(idims,^D); jsp^L=js^L+kr(idims,^D);
    iLp^L=iL^L+kr(idims,^D);
    ! 覆盖范围（WENO3 调用、第一级默认初始化；见下方第一级择优）
    iw3min^D=iLmin^D-kr(idims,^D); iw3max^D=iLmax^D+kr(idims,^D);
    iw3p^L=iw3^L+kr(idims,^D);

    ! THINC 陡度参数：加小量下限防止 β=0 时的除零
    thinc_beta1 = max(bvd_beta, 1.0d-10)
    thinc_beta2 = max(bvd_beta2, 1.0d-10)

    !====================================================================
    ! 高阶候选剖面两端值（由 bvd_highorder 选择，统一覆盖择优与传播需求）
    !====================================================================
    if (bvd_highorder==2) then
       ! -------- WENO3 高阶候选（bvd_highorder=2）--------
       ! 调用框架 WENO3limiter（与 limiter='weno3' 同款重构）：
       !   wLC3(j) = 面 j+1/2 的左侧值（单元 j 的 2 端）
       !   wRC3(j) = 面 j+1/2 的右侧值（单元 j+1 的 1 端）
       ! 预置一阶回退值（WENO3 内部 fix_limiter 对非物理解位置保留此值）
       wLC3(iw3^S,1:nwflux) = w(iw3^S,1:nwflux)
       wRC3(iw3^S,1:nwflux) = w(iw3p^S,1:nwflux)
       call WENO3limiter(ixI^L,iw3^L,idims,dxdim,w,wLC3,wRC3,1)
       ! 映射回单元端（视角翻转）：pend2(j)=wLC3(j)；pend1(j+1)=wRC3(j)
       pend2(iw3^S,1:nwflux) = wLC3(iw3^S,1:nwflux)
       pend1(iw3p^S,1:nwflux) = wRC3(iw3^S,1:nwflux)
    else
       ! -------- 二次多项式高阶候选（bvd_highorder=1，默认）--------
       ! 与原程序 second_degree_polynomial 逐位一致
       pend2(ic^S,1:nwflux) = (-1.0d0*w(icm^S,1:nwflux)+5.0d0*w(ic^S,1:nwflux) &
            +2.0d0*w(icp^S,1:nwflux))/6.0d0
       pend1(ic^S,1:nwflux) = (2.0d0*w(icm^S,1:nwflux)+5.0d0*w(ic^S,1:nwflux) &
            -1.0d0*w(icp^S,1:nwflux))/6.0d0
    end if
    !====================================================================
    ! BVD 选择机制开关（bvd_selection，namelist 可调）
    !   .true. （默认）：THINC 与高阶候选按 TBV 择优（完整 BVD 行为）
    !   .false.：跳过 THINC 与 TBV 择优，直接采用高阶候选（诊断/对照用）
    !====================================================================
    if (bvd_selection) then

       ! 第一级 THINC 剖面两端（完整 THINC，与原程序 thinc 逐位一致）
       call bvd_thinc(ixI^L,ic^L,idims,thinc_beta1,w,tend1,tend2)

       !====================================================================
       ! 第一级择优：全范围默认 + 命中覆盖（对应原程序 select_function_1st）
       !   默认：全部单元（含两侧各一个外延单元，即 iw3/iw3p——覆盖
       !         第二级 TBV 的全部读取范围）先取高阶候选剖面
       !   命中：TBV_thinc < TBV_poly（严格小于）的单元换成 THINC 剖面
       !   TBV(q,j) = |q_2(j)-q_1(j+1)| + |q_2(j-1)-q_1(j)|
       !====================================================================
       uend2(iw3^S,1:nwflux)=pend2(iw3^S,1:nwflux)
       uend1(iw3p^S,1:nwflux)=pend1(iw3p^S,1:nwflux)
       tbvP(js^S,1:nwflux) = abs(pend2(js^S,1:nwflux)-pend1(jsp^S,1:nwflux)) &
            +abs(pend2(jsm^S,1:nwflux)-pend1(js^S,1:nwflux))
       tbvT(js^S,1:nwflux) = abs(tend2(js^S,1:nwflux)-tend1(jsp^S,1:nwflux)) &
            +abs(tend2(jsm^S,1:nwflux)-tend1(js^S,1:nwflux))

       do iw=1,nwflux
          ! 命中单元换成 THINC（严格小于；未命中保持默认高阶候选）
          flagj(js^S)=.false.
          where (tbvT(js^S,iw) < tbvP(js^S,iw))
             uend2(js^S,iw)=tend2(js^S,iw)
             uend1(js^S,iw)=tend1(js^S,iw)
             flagj(js^S)=.true.
          end where
          ! 邻域传播：命中单元的左右邻居（i-1、i+1）也换成 THINC
          where (flagj(js^S))
             uend2(jsm^S,iw)=tend2(jsm^S,iw)
             uend1(jsm^S,iw)=tend1(jsm^S,iw)
             uend2(jsp^S,iw)=tend2(jsp^S,iw)
             uend1(jsp^S,iw)=tend1(jsp^S,iw)
          end where
       end do

       !====================================================================
       ! 第二级择优（可选，bvd_twostage=.true.）：
       ! 中间剖面与 THINC2(beta2) 再比一次 TBV，无邻域传播
       !====================================================================
       if (bvd_twostage) then
          ! THINC2 剖面两端
          call bvd_thinc(ixI^L,ic^L,idims,thinc_beta2,w,t2end1,t2end2)

          ! 中间剖面与 THINC2 各自的跨面跳变 TBV
          tbvU(js^S,1:nwflux) = abs(uend2(js^S,1:nwflux)-uend1(jsp^S,1:nwflux)) &
               +abs(uend2(jsm^S,1:nwflux)-uend1(js^S,1:nwflux))
          tbvT2(js^S,1:nwflux) = abs(t2end2(js^S,1:nwflux)-t2end1(jsp^S,1:nwflux)) &
               +abs(t2end2(jsm^S,1:nwflux)-t2end1(js^S,1:nwflux))

          ! THINC2 严格更小时覆盖中间剖面
          do iw=1,nwflux
             where (tbvT2(js^S,iw) < tbvU(js^S,iw))
                uend2(js^S,iw)=t2end2(js^S,iw)
                uend1(js^S,iw)=t2end1(js^S,iw)
             end where
          end do
       end if

    else

       ! 选择机制已关闭：跳过 THINC 与 TBV 择优，直接采用高阶候选剖面
       do iw=1,nwflux
          uend2(js^S,iw)=pend2(js^S,iw)
          uend1(js^S,iw)=pend1(js^S,iw)
       end do

    end if

    !====================================================================
    ! 写回界面值（注意单元端 1/2 与界面左/右的视角翻转）：
    !   wLC(j) = uend2(j)   ：单元 j 的 2 端（右端）= 面 j+1/2 的左侧值
    !   wRC(j) = uend1(j+1) ：单元 j+1 的 1 端（左端）= 面 j+1/2 的右侧值
    !====================================================================
    do iw=1,nwflux
       wLCtmp(iL^S,iw)=uend2(iL^S,iw)
       wRCtmp(iL^S,iw)=uend1(iLp^S,iw)
    end do

    !> 物理合法性回退：若最终重构态在某个界面位置非物理
    !> （负密度/负压力等），fix_limiter 保留调用前的一阶迎风值
    call fix_limiter(ixI^L,iL^L,wLCtmp,wRCtmp,wLC,wRC)

  end subroutine BVDlimiter

end module mod_bvd
