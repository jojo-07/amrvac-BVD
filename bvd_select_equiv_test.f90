program wheretest
  ! Compare formulations of the BVD first-stage selection logic:
  !   A: original program style -- default ALL cells to poly, then for each
  !      flagged cell i switch cells i-1, i, i+1 to THINC.
  !   B: previous mod_bvd.t style -- per-cell where (flag -> tend else pend),
  !      then masked shifted-where propagation to jsm/jsp.
  !   C: B plus full-span default init (uend2 on iw3, uend1 on iw3p = pend)
  !   D: new mod_bvd.t style -- full-span default on iw3/iw3p (no elsewhere),
  !      then flag -> THINC where, then masked shifted-where propagation.
  ! Ranges mimic iLmin=1, iLmax=4:
  !   js = 1..5, jsm = 0..4, jsp = 2..6, iw3 = 0..5, iw3p = 1..6
  ! Read cells (must match A): u2(0:5), u1(1:6).  u2(6)/u1(0): never read.
  implicit none
  integer :: ip, ibit, i, m1B, m2B, m1C, m2C, mD
  real :: t2(-1:9), p2(-1:9), t1(-1:9), p1(-1:9)
  real :: u2(-1:9), u1(-1:9), us2(-1:9), us1(-1:9), u2c(-1:9), u1c(-1:9)
  real :: u2d(-1:9), u1d(-1:9)
  logical :: flag(-1:9)
  real, parameter :: SENT = -999.0

  do i = -1, 9
     t2(i) = 100.0 + i
     p2(i) = 200.0 + i
     t1(i) = 300.0 + i
     p1(i) = 400.0 + i
  end do

  m1B = 0
  m2B = 0
  m1C = 0
  m2C = 0
  mD  = 0
  do ip = 0, 31
     flag = .false.
     do ibit = 0, 4
        if (btest(ip, ibit)) flag(1 + ibit) = .true.
     end do

     ! ---- A: original style (us = poly over full span, then switch) ----
     us2 = p2
     us1 = p1
     do i = 1, 5
        if (flag(i)) then
           us2(i) = t2(i)
           us1(i) = t1(i)
        end if
     end do
     do i = 1, 5
        if (flag(i)) then
           us2(i-1) = t2(i-1)
           us1(i-1) = t1(i-1)
           us2(i+1) = t2(i+1)
           us1(i+1) = t1(i+1)
        end if
     end do

     ! ---- B: previous mod_bvd.t (per-cell where then shifted propagation) ----
     u2 = SENT
     u1 = SENT
     where (flag(1:5))
        u2(1:5) = t2(1:5)
        u1(1:5) = t1(1:5)
     elsewhere
        u2(1:5) = p2(1:5)
        u1(1:5) = p1(1:5)
     end where
     where (flag(1:5))
        u2(0:4) = t2(0:4)
        u1(0:4) = t1(0:4)
        u2(2:6) = t2(2:6)
        u1(2:6) = t1(2:6)
     end where

     ! ---- C: B + full-span default (iw3 / iw3p set to pend first) ----
     u2c = SENT
     u1c = SENT
     u2c(0:5) = p2(0:5)
     u1c(1:6) = p1(1:6)
     where (flag(1:5))
        u2c(1:5) = t2(1:5)
        u1c(1:5) = t1(1:5)
     elsewhere
        u2c(1:5) = p2(1:5)
        u1c(1:5) = p1(1:5)
     end where
     where (flag(1:5))
        u2c(0:4) = t2(0:4)
        u1c(0:4) = t1(0:4)
        u2c(2:6) = t2(2:6)
        u1c(2:6) = t1(2:6)
     end where

     ! ---- D: new mod_bvd.t (default on iw3/iw3p, hit-where, propagation) ----
     u2d = SENT
     u1d = SENT
     u2d(0:5) = p2(0:5)
     u1d(1:6) = p1(1:6)
     where (flag(1:5))
        u2d(1:5) = t2(1:5)
        u1d(1:5) = t1(1:5)
     end where
     where (flag(1:5))
        u2d(0:4) = t2(0:4)
        u1d(0:4) = t1(0:4)
        u2d(2:6) = t2(2:6)
        u1d(2:6) = t1(2:6)
     end where

     ! compare 1: output cells 1..5
     do i = 1, 5
        if (abs(u2(i) - us2(i)) > 1.0e-5) m1B = m1B + 1
        if (abs(u1(i) - us1(i)) > 1.0e-5) m1B = m1B + 1
        if (abs(u2c(i) - us2(i)) > 1.0e-5) m1C = m1C + 1
        if (abs(u1c(i) - us1(i)) > 1.0e-5) m1C = m1C + 1
     end do
     ! compare 2: second-stage edge cells u2(0), u1(6)
     if (abs(u2(0) - us2(0)) > 1.0e-5) m2B = m2B + 1
     if (abs(u1(6) - us1(6)) > 1.0e-5) m2B = m2B + 1
     if (abs(u2c(0) - us2(0)) > 1.0e-5) m2C = m2C + 1
     if (abs(u1c(6) - us1(6)) > 1.0e-5) m2C = m2C + 1
     ! compare 3: D vs A on ALL read cells: u2(0:5) and u1(1:6)
     do i = 0, 5
        if (abs(u2d(i) - us2(i)) > 1.0e-5) mD = mD + 1
     end do
     do i = 1, 6
        if (abs(u1d(i) - us1(i)) > 1.0e-5) mD = mD + 1
     end do

     if (ip == 4) then   ! flag(3)=T only
        print *, 'pattern flag(3)=T only:'
        print *, '  B  u2(1:5) =', u2(1:5)
        print *, '  A  us2(1:5)=', us2(1:5)
     end if
     if (ip == 0) then   ! no flags
        print *, 'pattern no flags:'
        print *, '  B  u2(0)=', u2(0), '   A us2(0)=', us2(0)
        print *, '  C  u2c(0)=', u2c(0), '  C u1c(6)=', u1c(6)
        print *, '  D  u2d(0)=', u2d(0), '  D u1d(6)=', u1d(6)
     end if
  end do

  print *, '=== compare1 output cells (expect 0 / 0):  B =', m1B, '  C =', m1C
  print *, '=== compare2 edge cells (expect 32 / 0):   B =', m2B, '  C =', m2C
  print *, '=== compare3 read cells D vs A (expect 0): D =', mD
end program wheretest
