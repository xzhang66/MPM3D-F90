
!-------------------------------------------------------------------
!- MPM3D - A Three Dimensional Explicit Material Point Method Code -
!-                                                                 -
!- Developed by                                                    -
!-    Computational Dynamics Laboratory                            -
!-    School of Aerospace, Tsinghua University                     -
!-    Beijing 100084, China.                                       -
!-                                                                 -
!-    Professor Xiong Zhang                                        -
!-    Email: xzhang@tsinghua.edu.cn                                -
!-    Web  : http://www.comdyn.cn                                  -
!-                                                                 -
!-    Copyright (C) 2004 ~ 2012                                    -
!-------------------------------------------------------------------

module MaterialModel

  use MaterialData

  private  
  public Constitution

  integer:: mid       ! material set id
  integer:: etype_    ! Type of EOS
  integer:: mtype_    ! Type of material model

  real(8):: young_, poisson_
  real(8):: yield0_, tangmod_

  real(8):: den0_     ! Initial density
  real(8):: den_      ! Current density

  real(8):: vold      ! volume of step n
  real(8):: vol_      ! Current volume
  real(8):: vol0_     ! initial volume
  real(8):: dvol      ! 0.5 * volume increment 

  real(8):: dinc(6)   ! strain increment
  real(8):: sm        ! Mean stress
  real(8):: sd(6)     ! deviatoric stress      
  real(8):: sig(6)    ! stress components        
  real(8):: sold(6)   ! deviatoric stress of step n
  real(8):: dsm
  real(8):: bqf       ! bulk viscosity force

  real(8):: seqv      ! Equivalent stress
  real(8):: epeff_    ! Effective plastic strain
  real(8):: sig_y_    ! Current yield stress
  real(8):: depeff    ! increment of equivalent plastic strain
  real(8):: ratio     ! for hardening caculation

  real(8):: G2, K3, PlaMod ! 2*G, 3*K, plastic hardening modulus

  real(8):: specheat_
  real(8):: tmprt     ! temperature

  real(8):: iener     ! internal energy
  real(8):: specen    ! internal energy per initial volume
  real(8):: ieinc     ! internal energy increment

  real(8):: mu        ! mu = den_/den0_ - 1
  real(8):: rv        ! relative volume rv = vol_/vol0_
  real(8):: bfac      ! burn fraction
  real(8):: cp        ! sound speed

contains

  subroutine Constitution(de, vort, b, p)
!------------------------------------------------------------------
!-  Purpose                                                       -
!-      Update stresses by using a constitution model             -
!-  Inputs:                                                       -
!-      de - strain increment                                     -
!-             (D11, D22, D33, 2D23, 2D13, 2D12)*DT               -
!-      vort  - vorticity * DT (W32, W13, W21)*DT                 -
!-      b     - body index                                        -
!-      p     - particle index                                    -
!-    Outputs:                                                    -
!-        stress component                                        -
!-    Note:                                                       -
!-        sd(i) and de(i) comply the Voigt rule                   -
!------------------------------------------------------------------
    use ParticleData

    implicit none

    integer, intent(in):: b, p
    real(8), intent(in):: de(6), vort(3)

    real(8):: mp_
    real(8):: ltim
    logical:: failure

    ! Pick parameters    
    dinc = de

    ! Volume at time step t
    vold = particle_list(p)%VOL
    ! Current volume
    particle_list(p)%VOL = vold*(1+de(1)+de(2)+de(3))    
    vol_ = particle_list(p)%VOL
    dvol = 0.5 * (vol_ - vold)

    if (vol_ .lt. 0) then
        write(*,*) '=== warning: negative volume at particle', p
    end if
    mid = body_list(b)%mat

    epeff_ = particle_list(p)%epeff
    sig_y_ = particle_list(p)%sig_y    ! SIGY
    seqv = particle_list(p)%seqv       ! 2005-8-13

    ltim = particle_list(p)%LT         ! lighting time
    tmprt = particle_list(p)%celsius_t ! temperature

    iener = particle_list(p)%ie           ! internal energy
    mp_ = particle_list(p)%mass           ! particle mass

    failure = particle_list(p)%failure

    mtype_ = mat_list(mid)%MatType
    etype_ = mat_list(mid)%EosType

    young_ = mat_list(mid)%Young       ! Young's Modulus
    poisson_ = mat_list(mid)%Poisson   ! Poisson Ratio
    yield0_ = mat_list(mid)%Yield0     ! Yield limit
    tangmod_ = mat_list(mid)%TangMod   ! Tangential modulus

    den0_ = mat_list(mid)%Density      ! Initial density
    vol0_ = mp_/den0_                  ! Initial volume
    den_ = mp_/vol_                    ! Current density
    mu = den_/den0_ - 1
    rv = vol_/vol0_                    ! Relative volume
    specen = iener / vol0_


    young_ = young_*(1.0-particle_list(p)%DMG+epsilon(1.0d0))
    G2 = young_ / (1 + poisson_)
    K3 = young_ / (1 - 2*poisson_)
    PlaMod = young_ * tangmod_ / (young_ - tangmod_)

    sm = particle_list(p)%SM           ! Mean stress

    sd(1) = particle_list(p)%SDxx      ! deviatoric stress
    sd(2) = particle_list(p)%SDyy
    sd(3) = particle_list(p)%SDzz
    sd(4) = particle_list(p)%SDyz
    sd(5) = particle_list(p)%SDxz
    sd(6) = particle_list(p)%SDxy

    sig(1) = sd(1) + sm        ! cauchy stress
    sig(2) = sd(2) + sm
    sig(3) = sd(3) + sm
    sig(4) = sd(4)
    sig(5) = sd(5)
    sig(6) = sd(6)

    sold = sd    ! cauchy stress at time step t

    ! Select material model

    select case(mtype_)

    case(1) 
       ! elas: elastic model
       call sigrot(vort, sig, sm, sd)    ! Rotate stress
       call M3DM1()
       call lieupd()

    case(2) 
       ! pla1: elastic-perfectly plastic
       call sigrot(vort, sig, sm, sd)
       call M3DM2()
       call lieupd()

    case(3) 
       ! pla2: isotropic hardening
       call sigrot(vort, sig, sm, sd)
       call M3DM3()
       call lieupd()

    case(4) 
       ! john: Johnson-Cook plasticity model
       call sigrot(vort, sig, sm, sd)
       call M3DM4(mat_list(mid), DT, tmprt)
       call lieupd()
       ! call seleos()

       specheat_ = mat_list(mid)%SpecHeat
       particle_list(p)%celsius_t = particle_list(p)%celsius_t + &
                                    seqv*depeff/den_/specheat_

    case(5) 
       ! sjc: Simplified Johnson-Cook plasticity model
       call sigrot(vort, sig, sm, sd)
       call M3DM5(mat_list(mid), DT)
       call seleos(failure)

    case(6) 
       ! sjcf: Simplified Johnson-Cook plasticity model with 
       !       failure
       call sigrot(vort, sig, sm, sd)
       if (.not.particle_list(p)%failure) then
          call M3DM5(mat_list(mid), DT)    
          call seleos(failure)
       else
          call M3DM7()
          call seleos(failure)
       end if

       if (epeff_ .gt. mat_list(mid)%epf) then
          epeff_ = mat_list(mid)%epf + 0.0000001
          particle_list(p)%failure = .true.
       end if

    case(7) 
       ! john: Johnson-Cook plasticity model with failure
       call sigrot(vort, sig, sm, sd)
       if (.not.particle_list(p)%failure) then
          call M3DM8(mat_list(mid), DT, tmprt)
          call seleos(failure)
          specheat_ = mat_list(mid)%SpecHeat
          particle_list(p)%celsius_t = particle_list(p)%celsius_t &
                                       + seqv*depeff/den_/specheat_
       else
          call M3DM7()
          call seleos(failure)
       end if

       if (epeff_ .gt. mat_list(mid)%epf) then
          epeff_ = mat_list(mid)%epf + 0.0000001
          failure = .true.
       end if
       
       if (failure .and. sm < 0) then
          sm=0.0
          particle_list(p)%VOL = vold
       end if
       particle_list(p)%failure = failure

    case(8) 
       ! hiex: High Explosive burn
       call M3DM6(ltim, CurrentTime, mat_list(mid)%D)
       call seleos(failure)

    case(9)
       ! null: used to model air
       call M3DM7()
       call seleos(failure)

    case(10) 
       ! Drucker-Prager elastic-perfectly plastic
       call sigrot(vort, sig, sm, sd)
       call M3DM9(mat_list(mid), DT)
       call lieupd()

    case default 
       write(*, 10) mtype_
10     format(1x,'*** Stop *** material type ', i2, &
              ' has not been implemented !')
       stop

    end select

    ! Write stress result

    particle_list(p)%SM = sm
    particle_list(p)%Seqv = seqv

    particle_list(p)%SDxx = sd(1)
    particle_list(p)%SDyy = sd(2)
    particle_list(p)%SDzz = sd(3)
    particle_list(p)%SDyz = sd(4)
    particle_list(p)%SDxz = sd(5)
    particle_list(p)%SDxy = sd(6)

    particle_list(p)%sig_y = sig_y_
    particle_list(p)%epeff = epeff_

    particle_list(p)%ie = iener
    particle_list(p)%cp = cp

  end subroutine Constitution


  subroutine sigrot(vort, sig, sm, sd)
!------------------------------------------------------------------
!-  Purpose                                                       -
!-      Rotate stresses, and then update the mean stress sm       -
!-      and deviate stresses sd                                   -
!-  Input                                                         -
!-      vort - Vorticity increments (W32, W13, W21)*DT            -
!-      sig  - Caucy stresses at time step t                      -
!-             (S11, S22, S33, S23, S13, S12)                     -
!-  Output                                                        -
!-      sig  - Rotated caucy stresses                             -
!-      sd   - Rotated deviate stress                             -
!-      sm   - Rotated mean stress                                -
!-  References                                                    -
!-      Section 5.1                                               -
!------------------------------------------------------------------
    implicit none
    real(8), intent(in) :: vort(3)
    real(8), intent(inout) :: sig(6)
    real(8), intent(out) :: sm, sd(6)

    real(8) :: rot(6), q(3)

    q(1) = 2d0*sig(6)*vort(3)
    q(2) = 2d0*sig(5)*vort(2)
    q(3) = 2d0*sig(4)*vort(1)

    rot(1) = - q(1) + q(2) ! (Eq: 5.4)
    rot(2) = + q(1) - q(3)
    rot(3) = - q(2) + q(3)
    rot(4) = vort(1)*(sig(2)-sig(3)) + vort(3)*sig(5) - &
             vort(2)*sig(6) 
    rot(5) = vort(2)*(sig(3)-sig(1)) + vort(1)*sig(6) - &
             vort(3)*sig(4)
    rot(6) = vort(3)*(sig(1)-sig(2)) + vort(2)*sig(4) - &
             vort(1)*sig(5)

    sig = sig + rot ! First two terms in RHS of (Eq. 5.3)

    sm = (sig(1)+sig(2)+sig(3))/3d0

    sd(1) = sig(1) - sm
    sd(2) = sig(2) - sm
    sd(3) = sig(3) - sm
    sd(4) = sig(4)
    sd(5) = sig(5)
    sd(6) = sig(6)

  end subroutine sigrot


  subroutine elastic_devi()
!------------------------------------------------------------------
!-  Purpose                                                       -
!-      Update devitoric stress by elastic relation               -
!-  Inputs                                                        -
!-      dinc - strain increment                                   -
!-             (D11, D22, D33, 2D23, 2D13, 2D12)*DT               -
!-  Outputs                                                       -
!-    sd     - devitoric stress component                         -
!-  References                                                    -
!-      Section 5.2.1                                             -
!------------------------------------------------------------------
    implicit none
    real(8):: dem, G

    dem = (dinc(1) + dinc(2) + dinc(3)) / 3.0

    G = 0.5*G2    ! Shear modulus

    sd(1) = sd(1) + G2*(dinc(1)-dem)    ! (Eq: 5.29)
    sd(2) = sd(2) + G2*(dinc(2)-dem)
    sd(3) = sd(3) + G2*(dinc(3)-dem)
    sd(4) = sd(4) + G*dinc(4)
    sd(5) = sd(5) + G*dinc(5)
    sd(6) = sd(6) + G*dinc(6)

  end subroutine elastic_devi


  subroutine elastic_p()
!------------------------------------------------------------------
!-  Purpose                                                       -
!-      Update pressure by elastic relation                       -
!-  Inputs                                                        -
!       dinc - strain increment                                   -
!   Outputs                                                       -
!       sm   - mean stress (pressure)                             -
!------------------------------------------------------------------
    implicit none
    real(8):: dem

    dem = (dinc(1) + dinc(2) + dinc(3)) / 3.0
    dsm = K3*dem

    sm = sm + dsm        ! (Eq: 5.30)

  end subroutine elastic_p


  subroutine lieupd()
!------------------------------------------------------------------
!-  Purpose                                                       -
!-      Update energy for material models without EOS             -
!------------------------------------------------------------------
    use ParticleData, only: EngInternal
    implicit none

    real:: vavg

    vavg = vol_ + vold
    iener = iener + 0.25d0*ieinc*vavg ! (Eq. 5.10)

  end subroutine lieupd


  subroutine hieupd()
!------------------------------------------------------------------
!-  Purpose                                                       -
!-      Update energy for material models calling a EOS           -
!-  References                                                    -
!-      Section 5.1                                               -
!------------------------------------------------------------------
    implicit none
    real:: vavg

    vavg = vol_ + vold

    ieinc = dinc(1)*(sold(1)+sd(1)) + dinc(2)*(sold(2)+sd(2)) +  &
         dinc(3)*(sold(3)+sd(3)) + dinc(4)*(sold(4)+sd(4)) +     &
         dinc(5)*(sold(5)+sd(5)) + dinc(6)*(sold(6)+sd(6))

    iener = iener + 0.25*ieinc*vavg + dvol*sm    ! (Eq: 5.15)

    specen = iener / vol0_

  end subroutine hieupd


  real function EquivalentStress()
!------------------------------------------------------------------
!-  Purpose                                                       -
!-      Calculate the equivalent stress sqrt(3*J2)                -
!-  Input                                                         -
!-      sd    - the deviatoric stress components                  -
!-  Return                                                        -
!-      EquivalentStress - the equivalent stress                  -
!------------------------------------------------------------------
    use ParticleData
    implicit none

    real(8) :: J2  ! the second stress invariant

    J2 = 0.5d0*(sd(1)**2 + sd(2)**2 + sd(3)**2) + sd(4)**2 + &
         sd(5)**2 + sd(6)**2
    J2 = J2*3.0
    EquivalentStress = sqrt(J2)     

    return
  end function EquivalentStress


  subroutine M3DM1()
!------------------------------------------------------------------
!-  Purpose                                                       -
!-      Elastic material model                                    -
!-  References                                                    -
!-      Section 5.2.1                                             -
!------------------------------------------------------------------
    implicit none

    ieinc = sig(1)*dinc(1) + sig(2)*dinc(2) + sig(3)*dinc(3) +   &
         sig(4)*dinc(4) + sig(5)*dinc(5) + sig(6)*dinc(6)    

    call elastic_devi()
    call elastic_p()

    seqv = EquivalentStress()

    ieinc = ieinc + (sd(1)+sm)*dinc(1) + (sd(2)+sm)*dinc(2) +   &
            (sd(3)+sm)*dinc(3) +sd(4)*dinc(4) + sd(5)*dinc(5) + &
            sd(6)*dinc(6)     ! (part of Eq: 5.10)

  end subroutine M3DM1


  subroutine M3DM2()
!------------------------------------------------------------------
!-  Purpose                                                       -
!-      Elastic-perfectly plastic material model                  -
!-  References                                                    -
!-      Section 5.2.4                                             -
!------------------------------------------------------------------
    implicit none

    ieinc = sig(1)*dinc(1) + sig(2)*dinc(2) + sig(3)*dinc(3) +   &
         sig(4)*dinc(4) + sig(5)*dinc(5) + sig(6)*dinc(6)

    call elastic_devi()

    seqv = EquivalentStress()

    if (seqv .GT. yield0_) then
       depeff = (seqv - sig_y_) / (1.5e0*G2)    ! (Eq: 5.93)
       epeff_ = epeff_ + depeff                 ! (Eq: 5.98)

       ratio = yield0_/seqv    ! (Eq: 5.100)

       sd(1) = sd(1)*ratio     ! (Eq: 5.101)
       sd(2) = sd(2)*ratio
       sd(3) = sd(3)*ratio
       sd(4) = sd(4)*ratio
       sd(5) = sd(5)*ratio
       sd(6) = sd(6)*ratio

       seqv = seqv*ratio       ! (Eq: 5.102)
    end if

    call elastic_p()

    ieinc = ieinc + (sd(1)+sm)*dinc(1) + (sd(2)+sm)*dinc(2) + &
            (sd(3)+sm)*dinc(3) + sd(4)*dinc(4) +  &
            sd(5)*dinc(5) + sd(6)*dinc(6)   ! (part of Eq: 5.10)

  end subroutine M3DM2


  subroutine M3DM3()
!------------------------------------------------------------------
!-  Purpose                                                       -
!-      Isotropic hardening plastic material model                -
!-      LS-DYNA theorectical manual : Material Model 10           -
!-  References                                                    -
!-      Section 5.2.4                                             -
!------------------------------------------------------------------
    implicit none

    ieinc = sig(1)*dinc(1) + sig(2)*dinc(2) + sig(3)*dinc(3) +   &
         sig(4)*dinc(4) + sig(5)*dinc(5) + sig(6)*dinc(6)

    call elastic_devi()

    seqv = EquivalentStress()

    if (seqv .GT. sig_y_) then
       depeff = (seqv - sig_y_) / (1.5e0*G2 + PlaMod)    ! (Eq: 5.93)
       epeff_ = epeff_ + depeff                          ! (Eq: 5.98)

       sig_y_ = sig_y_ + PlaMod*depeff    ! (Eq: 5.99)
       ratio = sig_y_/seqv                ! (Eq: 5.100)

       sd(1) = sd(1)*ratio                ! (Eq: 5.101)
       sd(2) = sd(2)*ratio
       sd(3) = sd(3)*ratio
       sd(4) = sd(4)*ratio
       sd(5) = sd(5)*ratio
       sd(6) = sd(6)*ratio

       seqv = seqv*ratio                
    end if

    call elastic_p()

    ieinc = ieinc + (sd(1)+sm)*dinc(1) + (sd(2)+sm)*dinc(2) +  &
            (sd(3)+sm)*dinc(3) + sd(4)*dinc(4) + &
            sd(5)*dinc(5) + sd(6)*dinc(6)

  end subroutine M3DM3


  subroutine M3DM4(mat, DT, tmprt)
!------------------------------------------------------------------
!-  Purpose                                                       -
!-      Johnson-Cook plastic material model (john)                -
!-              (G.R. Johnson 1988)                               -
!-      LS-DYNA theorectical manual : Material Model 15           -
!-  References                                                    -
!-      Section 5.2.5                                             -
!------------------------------------------------------------------
    implicit none

    real(8), intent(in) :: DT, tmprt
    type(material), intent(in) :: mat
    real(8):: Bjc, njc, Cjc, mjc, epso, tstar

    ieinc = sig(1)*dinc(1) + sig(2)*dinc(2) + sig(3)*dinc(3) +  &
         sig(4)*dinc(4) + sig(5)*dinc(5) + sig(6)*dinc(6)

    Bjc = mat%B_jc
    njc = mat%n_jc
    Cjc = mat%C_jc
    mjc = mat%m_jc
    epso = mat%epso

    epeff_ = epeff_ + 0.0001
    ! Note: 
    !    PlaMod=Bjc*njc*(epeff_**(njc-1))*
    !           (1 + Cjc*log(depeff/epso/DT))*(1-1 - tstar**mjc)
    ! simplied as follow:
    PlaMod = Bjc*njc*(epeff_**(njc-1)) 
    epeff_ = epeff_ - 0.0001

    call elastic_devi()
    seqv = EquivalentStress()

    if (seqv .GT. sig_y_) then
       depeff = (seqv - sig_y_) / (1.5e0*G2 + PlaMod)
       epeff_ = epeff_ + depeff

       tstar = (tmprt-mat%roomt)/(mat%melt-mat%roomt)

       ! (Eq: 5.107)
	   sig_y_ = (yield0_ + Bjc*(epeff_**njc)) * &
                (1 + Cjc*log(depeff/epso/DT)) * (1 - tstar**mjc)

       ratio = sig_y_/seqv    ! (Eq: 5.100)

       sd(1) = sd(1)*ratio    ! (Eq: 5.101)
       sd(2) = sd(2)*ratio
       sd(3) = sd(3)*ratio
       sd(4) = sd(4)*ratio
       sd(5) = sd(5)*ratio
       sd(6) = sd(6)*ratio

       seqv = seqv*ratio      ! (Eq: 5.102)
    end if

    call elastic_p()

    ieinc = ieinc + (sd(1)+sm)*dinc(1) + (sd(2)+sm)*dinc(2) + &
            (sd(3)+sm)*dinc(3) +  sd(4)*dinc(4) + &
            sd(5)*dinc(5) + sd(6)*dinc(6)  ! (part of Eq: 5.10)

  end subroutine M3DM4


  subroutine M3DM5(mat, DT)
!------------------------------------------------------------------
!-  Purpose                                                       -
!-      Simplified Johnson-Cook plastic material model (sjc)      -
!-                 Ignoring the temperature effects in            -
!-                 Johnson-Cook model                             -
!-      LS-DYNA theorectical manual : Material Model 98           -
!-      should be used with EOS                                   -
!-  References                                                    -
!-      Section 5.2.4                                             -
!------------------------------------------------------------------
    implicit none

    real(8), intent(in) :: DT
    type(material), intent(in) :: mat
    real(8):: Bjc, njc, Cjc, epso

    Bjc = mat%B_jc
    njc = mat%n_jc
    Cjc = mat%C_jc
    epso = mat%epso

    epeff_ = epeff_ + 0.0001
    PlaMod = Bjc*njc*(epeff_**(njc-1))
    epeff_ = epeff_ - 0.0001

    call elastic_devi()

    seqv = EquivalentStress()

    if (seqv .GT. sig_y_) then
       depeff = (seqv - sig_y_) / (1.5e0*G2 + PlaMod)
       epeff_ = epeff_ + depeff

       sig_y_ = (yield0_ + Bjc*(epeff_**njc)) * &
                (1 + Cjc*log(depeff/epso/DT))

       ratio = sig_y_/seqv

       sd(1) = sd(1)*ratio
       sd(2) = sd(2)*ratio
       sd(3) = sd(3)*ratio
       sd(4) = sd(4)*ratio
       sd(5) = sd(5)*ratio
       sd(6) = sd(6)*ratio

       seqv = seqv*ratio
    end if

  end subroutine M3DM5


  subroutine M3DM6(ltim,ctim,dvelo)
!------------------------------------------------------------------
!-  Purpose                                                       -
!-      High Explosive burn material                              -
!-      an Equation of State must be defined                      -
!-      LS-DYNA theorectical manual : Material Model 8            -
!-  References                                                    -
!-      Section 5.2.13                                            -
!------------------------------------------------------------------
    implicit none
    real(8), intent(in) :: ltim, ctim, dvelo

    if (ctim.gt.ltim) then
       bfac = 1.0
    else
       bfac = 0.0
    end if

  end subroutine M3DM6


  subroutine M3DM7()
!------------------------------------------------------------------
!-  Purpose                                                       -
!-      Null material model                                       -
!------------------------------------------------------------------
    implicit none

    sd = 0.0

  end subroutine M3DM7

  subroutine M3DM8(mat, DT, tmprt)
!------------------------------------------------------------------
!-  Purpose                                                       -
!-      Johnson-Cook plastic material model (john)                -
!-      (G.R. Johnson 1988)                                       -
!-      LS-DYNA theorectical manual : Material Model 15           -
!-      should be used with EOS                                   -
!-  References                                                    -
!-      Section 5.2.5                                             -
!------------------------------------------------------------------
    implicit none

    real(8), intent(in) :: DT, tmprt
    type(material), intent(in) :: mat
    real(8):: Bjc, njc, Cjc, mjc, epso, tstar

    Bjc = mat%B_jc
    njc = mat%n_jc
    Cjc = mat%C_jc
    mjc = mat%m_jc
    epso = mat%epso

    epeff_ = epeff_ + 0.0001
    PlaMod = Bjc*njc*(epeff_**(njc-1))
    epeff_ = epeff_ - 0.0001


    call elastic_devi()
    seqv = EquivalentStress()

    if (seqv .GT. sig_y_) then
       depeff = (seqv - sig_y_) / (1.5e0*G2 + PlaMod)
       epeff_ = epeff_ + depeff

       tstar = (tmprt-mat%roomt)/(mat%melt-mat%roomt)

       sig_y_ = (yield0_ + Bjc*(epeff_**njc)) * &
                (1 + Cjc*log(depeff/epso/DT)) * (1 - tstar**mjc)

       ratio = sig_y_/seqv

       sd(1) = sd(1)*ratio
       sd(2) = sd(2)*ratio
       sd(3) = sd(3)*ratio
       sd(4) = sd(4)*ratio
       sd(5) = sd(5)*ratio
       sd(6) = sd(6)*ratio

       seqv = seqv*ratio
    end if

  end subroutine M3DM8

  subroutine M3DM9(mat, DT)
!------------------------------------------------------------------
!-  Purpose                                                       -
!-      perfect elastic-plastic Drucker-Prager model for soil     -
!-      LS-DYNA theorectical manual : Material Model 193          -
!-      should be used without EOS                                -
!-  References                                                    -
!-      Section 5.2.6                                             -
!------------------------------------------------------------------
    implicit none

    real(8), intent(in) :: DT
    type(material), intent(in) :: mat
    real(8):: qfai, kfai, qpsi, tenf, Gmod, Kmod,tenf_max 
    real(8):: J2, Tau ,dpFi, dpsig ,dlamd ,newTau
    real(8):: dp_hfai, Taup,alphap
    integer :: iplas

    Gmod = young_ / (2* (1 + poisson_))
    Kmod =  young_ / (3*(1 - 2*poisson_))
    qfai = mat%q_fai
    kfai = mat%k_fai
    qpsi = mat%q_psi

    ! internal energy incremental
    ieinc = sig(1)*dinc(1) + sig(2)*dinc(2) + sig(3)*dinc(3) +  &
         sig(4)*dinc(4) + sig(5)*dinc(5) + sig(6)*dinc(6)

    ! Give the tension stress value
    !* --- set fTension to cone apex if larger than apex --- *
    if (qfai == 0.0 )then
       tenf = 0.0
    else
       tenf_max =  kfai / qfai
       tenf = min(mat%ten_f,tenf_max)
    end if

    ! --- trial elastic stresses --- 
    call elastic_devi()
    call elastic_p()

    iplas = 0         ! elastic calculation
    J2 = 0.5d0*(sd(1)**2 + sd(2)**2 + sd(3)**2) + sd(4)**2 + &
         sd(5)**2 + sd(6)**2
    Tau = sqrt(J2)
    seqv = Tau*sqrt(3.0)
    dpFi = Tau + qfai*sm - kfai   ! D-P yield surface
    dpsig = sm - tenf             ! the spherical stress difference

    if (dpsig < 0.0) then         
       if (dpFi >0.0) then
          iplas = 1 ! shear plastic flow
          ! plastic flow coefficient
          dlamd = dpFi/(Gmod + Kmod*qfai*qpsi) 
          ! correct spherical stress
          sm = sm - Kmod*qpsi*dlamd            
          ! correct shear stress 
          newTau = kfai - qfai*sm              
          ratio = newTau/Tau

          ! correct deviatoric stress
          sd(1) = sd(1)*ratio
          sd(2) = sd(2)*ratio
          sd(3) = sd(3)*ratio
          sd(4) = sd(4)*ratio
          sd(5) = sd(5)*ratio
          sd(6) = sd(6)*ratio

          seqv =  seqv*ratio     !correct the mises stress         

          ! calculate the effective plastic strain
          depeff = dlamd*sqrt(1.0/3.0 + (2.0/9.0)*(qpsi**2))
          epeff_ = epeff_ + depeff

       end if
    else   !(dpsig >= 0.0)then
       alphap = sqrt(1 + qfai**2) - qfai
       Taup = kfai - qfai*tenf
       dp_hfai = Tau - Taup - alphap * dpsig

       if(dp_hfai > 0.0) then
          iplas = 1 ! shear plastic flow
          ! plastic flow coefficient
          dlamd = dpFi/(Gmod + Kmod*qfai*qpsi)    
          ! correct spherical stress 
          sm = sm - Kmod*qpsi*dlamd               
          ! correct shear stress
          newTau = kfai - qfai*sm                 
          ratio = newTau/Tau

          ! correct deviatoric stress
          sd(1) = sd(1)*ratio
          sd(2) = sd(2)*ratio
          sd(3) = sd(3)*ratio
          sd(4) = sd(4)*ratio
          sd(5) = sd(5)*ratio
          sd(6) = sd(6)*ratio

          seqv =  seqv*ratio     !correct the mises stress         
          ! calculate the effective plastic strain

          depeff = dlamd*sqrt(1.0/3.0 + (2.0/9.0)*(qpsi**2))
          epeff_ = epeff_ + depeff
       else
          iplas = 2 ! tension plastic flow
          dlamd = (sm - tenf)/Kmod
          sm = tenf
          depeff = dlamd*(1.0/3.0)*sqrt(2.0)
          epeff_ = epeff_ + depeff
       end if
    end if

    ieinc = ieinc + (sd(1)+sm)*dinc(1) + (sd(2)+sm)*dinc(2) + &
            (sd(3)+sm)*dinc(3) + sd(4)*dinc(4) + &
            sd(5)*dinc(5) + sd(6)*dinc(6)

  end subroutine M3DM9


  subroutine seleos(failure)
!------------------------------------------------------------------
!-  Purpose                                                       -
!-      Update pressure and energy using appropriate              -
!-      equation of state                                         -
!-  Input                                                         -
!-      etype - Type of equation of state                         -
!-      failure - failure state of particle                       -
!------------------------------------------------------------------
    implicit none
    logical:: failure

    select case(etype_)

    case(1)
       call eos1(failure)
    case(2)
       call eos2(failure)
    case(3)
       call eos3(failure)
    case default
       stop "error eos id"
    end select

  end subroutine seleos


  subroutine eos1(failure)
!------------------------------------------------------------------
!-  Purpose                                                       -
!-      Linear Equation Of State                                  -
!-  Inputs                                                        -
!-      mu      - ( = rho/rho0 - 1 )                              -
!-      iener   - internal energy                                 -
!-      den_    - density                                         -
!-  Outputs                                                       -
!-      sm      - mean stress (pressure)                          -
!-  References                                                    -
!-      Section 5.2.3                                             -
!------------------------------------------------------------------
    use ParticleData, only: EngInternal
    implicit none
    real:: A, B, pnew, c0, c1, c2, c3, c4, c5, c6
    logical:: failure

    c0 = mat_list(mid)%cEos(1)
    c1 = mat_list(mid)%cEos(2)
    c2 = mat_list(mid)%cEos(3)
    c3 = mat_list(mid)%cEos(4)
    c4 = mat_list(mid)%cEos(5)
    c5 = mat_list(mid)%cEos(6)
    c6 = mat_list(mid)%cEos(7)

    A = c0 + mu * (c1 + mu * (c2 + mu * c3))
    B = c4 + mu * (c5 + mu * c6)

    call hieupd()

    pnew = (A + B * specen) / (1 + B * dvol / vol0_) !(Eq: 5.18)
    
    if(failure .and. pnew < 0) pnew=0.0 

    iener = iener - dvol * pnew        ! (Eq: 5.14)

    sm = -pnew

  end subroutine eos1


  subroutine eos2(failure)
!------------------------------------------------------------------
!-  Purpose                                                       -
!-      Mie-Gruneisen Equation Of State                           -
!-  Inputs                                                        -
!-      mu      - ( = rho/rho0 - 1 )                              -
!-      iener   - internal energy                                 -
!-      den_    - density                                         -
!-  Outputs                                                       -
!-      sm      - mean stress (pressure)                          -
!-  References                                                    -
!-      Section 5.3.5                                             -
!------------------------------------------------------------------
    use ParticleData, only: EngInternal
    implicit none
    real:: A, B, pnew, c1, c2, c3, c4, c5
    logical:: failure

    c1 = mat_list(mid)%cEos(1)
    c2 = mat_list(mid)%cEos(2)
    c3 = mat_list(mid)%cEos(3)
    c4 = mat_list(mid)%cEos(4)
    c5 = mat_list(mid)%cEos(5)

    if (mu .gt. 0) then      ! compressed
       A = (c1*mu+c2*mu**2 + c3*mu**3) * (1-c4*mu*0.5d0/den_)
       B = c5
    else          ! expanded
       !mu = max(mu,-0.5)
       A = c1*mu
       B = 0
    end if

    call bulkq()
    call hieupd()

    pnew = (A + B * specen) / (1 + B * dvol / vol0_) !(Eq: 5.18)

    if(failure .and. pnew < 0) pnew=0.0 

    iener = iener - dvol * pnew        ! (Eq: 5.14)

    sm = -pnew

  end subroutine eos2


  subroutine eos3(failure)
!------------------------------------------------------------------
!-  Purpose                                                       -
!-      JWL Equation Of State                                     -
!-  Inputs                                                        -
!-      rv    - Relative volume rv = vol_/vol0_                   -
!-      Eng - specific internal energy                            -
!-             (internal energy per unit mass)                    -
!-  Outputs                                                       -
!-      sm  - mean stress (pressure)                              -
!-  References                                                    -
!-      Section 5.2.4                                             -
!------------------------------------------------------------------
    use ParticleData, only: EngInternal
    implicit none

    real(8):: r1v, r2v, wr1v, wr2v, wdr1v, wdr2v, er1v, er2v
    real(8):: A, B, pnew, c1, c2, r1, r2, w1
    logical:: failure

    ! initialize parameters
    c1 = mat_list(mid)%cEos(1)
    c2 = mat_list(mid)%cEos(2)
    r1 = mat_list(mid)%cEos(3)
    r2 = mat_list(mid)%cEos(4)
    w1 = mat_list(mid)%cEos(5)

    r1v = r1*rv
    r2v = r2*rv
    wr1v = c1*w1/r1v
    wr2v = c2*w1/r2v
    wdr1v = c1 - wr1v
    wdr2v = c2 - wr2v
    er1v = exp(-r1v)
    er2v = exp(-r2v)

    ! calculate bulk viscosity

    A = wdr1v*er1v + wdr2v*er2v + bqf
    B = w1/rv

    call hieupd()

    A = A * bfac
    B = B * bfac

    pnew = (A + B * specen) / (1 + B * dvol / vol0_) ! (Eq: 5.18)
    
    if(failure .and. pnew < 0) pnew=0.0 
    iener = iener - dvol * pnew        ! (Eq: 5.14)

    sm = -pnew

  end subroutine eos3


  subroutine bulkq()
!------------------------------------------------------------------
!-  Purpose                                                       -
!-      bulk viscosity                                            -
!-  Inputs                                                        -
!-      bq1,bq2                                                   - 
!-  Outputs                                                       -
!-      bqf                                                       -
!------------------------------------------------------------------
    use GridData, only: DCell
    use ParticleData, only: DT
    implicit none
    real(8):: dd    ! bulk strain rate

    dd = (dinc(1) + dinc(2) + dinc(3)) / DT
    dd = min(dd,0.0)

    bqf = den_ * DCell * DCell * bq1 * dd * dd - &
          bq2 * den_ * DCell * cp * dd      !(Eq. 2.82)

  end subroutine bulkq

end module MaterialModel
