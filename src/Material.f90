
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
! ------------------------------------------------------------------
! -                                                                -
! -  Material parameters                                           -
! -                                                                -
! ------------------------------------------------------------------
module MaterialData

  type Material
     integer:: MatType    ! Material type
     integer:: EosType    ! EOS type
     real(8):: Density    ! Initial density
     real(8):: Young      ! Young's modulus
     real(8):: Poisson    ! Poisson's ratio
     real(8):: Mp         ! Partical mass
     real(8):: Yield0     ! Initial yield stress
     real(8):: TangMod    ! Tangential modulus
     real(8):: roomt      ! Room temperature
     real(8):: Melt       ! Melting temperature
     real(8):: SpecHeat   ! Specified heat capacity
     ! Parameters for Johnson-Cook material
     real(8):: B_jc, n_jc, C_jc, m_jc    
     ! strain rate normalization factor used in J-C model
     real(8):: epso    
     real(8):: prd      ! tensile pressure to begin damage
     ! tensile pressure and equivalent plastic strain at failure
     real(8):: epf    
     real(8):: D        ! Detonation velocity
     real(8):: cEos(10) = 0    ! Constants in Equation of State
     real(8):: Wavespd         ! wavespeed of Null material
     ! Parameters for Drucker-Prager soil material
     real(8):: q_fai, k_fai, q_psi, ten_f  
  end type Material
  ! Linear polynomial - C0, C1, C2, C3, C4, C5, C6
  ! Gruneisen - C1, C2, C3, C4
  ! JWL - A, B, R1, R2, W, E0, V0

  integer:: nb_mat = 0    ! number of material sets
  type(Material), allocatable:: mat_list(:)

  logical:: Jaum = .true.    ! use Jaumann rate ?


  integer, parameter:: maxDeto = 10
  integer:: nDeto = 0
  real(8):: DetoX(maxDeto) = 0.0,  & ! detonation point
            DetoY(maxDeto) = 0.0,  &
            DetoZ(maxDeto) = 0.0    

  ! artificial bulk viscosity coefficients
  real(8):: bq1=0.0, bq2=0.0    

end module MaterialData
