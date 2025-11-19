
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

!------------------------------------------------------------------
!-   Data input procedures                                        -
!------------------------------------------------------------------

module DataIn

  use ParticleData
  use GridData
  use DataOut
  use MaterialData

  integer:: bodyCounter = 0 ! body counter
  integer:: comCounter = 0  ! component counter
  integer:: parCounter = 0  ! particle counter

contains

  subroutine InputPara()
!------------------------------------------------------------------
!-    purpose: input data both from file and by program           -
!------------------------------------------------------------------

    use ParticleData
    use GridData   
    use FFI

    implicit none

    integer :: idx

    call GETARG(1,jobname)

    idx = index(jobname,".")
    if (idx > 0) then
        fName = jobname(1:index(jobname,".")-1)
    else
        fName = jobname
    endif

    FileInp = trim(fName) // ".mpm"    ! Input data
    FileOut = trim(fName) // ".out"    ! Message file
    FileCurv = trim(fName) // "_curv.dat" ! Curve data (TecPlot)

    open(iow2, file = FileCurv, status = 'unknown')

    call FFIOpen()
    call GetData()
    close(iord)        ! Close input file

    ! Write results in TecPlot format
    if (WriteTecPlot) then
        FileAnim = trim(fName) // "_anim.dat" ! Animation data (TecPlot)
        open(iow1, file = FileAnim, status = 'unknown')
    endif

    ! Write results in ParaView format
    if (WriteParaView) then
        FileAnim = trim(fName) // "_anim.vtu"
        FileAnimColl = trim(fName) // "_anim.pvd"
        open(iow12, file = FileAnimColl, status = 'unknown')

        ! Open XML Tags in Collection file
		write(iow12,300)
        write(iow12,301)
        write(iow12,302)
300     format('<?xml version="1.0"?>')
301     format('<VTKFile type="Collection" version="0.1" byte_order="BigEndian" compressor="vtkZLibDataCompressor">')
302     format('    <Collection>')
    end if

    open(iow03, file = 'EnergyPlot.dat', status = 'unknown')
    open(iow04, file = 'MomentumPlot.dat', status = 'unknown')
    open(iow05, file = 'ContforcPlot.dat', status = 'unknown')

    call Initial()

    call SetGridData()
    if (WriteParaView) call OutGrid()

    ! Create mutiple grid for contact
    if(contact) call SetContact_GridNodeData() 

    call SetDT()

    call statinfo()

  end subroutine InputPara

  subroutine GetData()
!------------------------------------------------------------------
!-  purpose: Get input data using FFI module                      -
!-           data are stored in ParticleData module               -
!------------------------------------------------------------------
    use FFI
    implicit none
    integer key, i, tempI
    logical:: updateMethod = .false.

    integer,parameter:: nbkw = 43
    character(4),parameter:: kw(nbkw) = (/ &
         'endi','mpm3','nbmp','endt','grid','spx ', &
         'spy ','spz ','dcel','dtsc','outt','rptt', &
         'fixe','nmat','mate','part','musl','jaum', &
         'load','velo','outr','curv','seos','pt2d', &
         'curx','deto','tecp','bulk','gimp','cont', &
         'usf ','usl ','nbco','nbbo','para','drda', &
         'bimp','sgmp','smoo','TLmp','UVFp','QuLo',&
         'Damp'/)

    do while(.true.)
       key = keyword(kw,nbkw)
       select case(key)

       case(1)    ! end input
          exit    ! Terminates execution of Do loop

       case(2)    ! mpm3 (title)
          call GetString(Title)

       case(3)    ! nbmp - number of material particles
          nb_particle = GetInt()
          write(*,"(a,i12)") 'Number of particles = ',nb_particle
          write(iomsg,"(a,i12)") 'Number of particles = ',nb_particle
          allocate(particle_list(nb_particle))
          call InitParticle()

       case(4)    ! endtime - Simulation end time
          EndTime = GetReal()
          OutTime = EndTime

       case(5)    ! grid - Define computational grid
          SpanX(1) = GetReal()
          SpanX(2) = GetReal()

          SpanY(1) = GetReal()
          SpanY(2) = GetReal()

          SpanZ(1) = GetReal()
          SpanZ(2) = GetReal()

       case(6)    ! spx - Define computational grid in x direction
          SpanX(1) = GetReal()
          SpanX(2) = GetReal()

       case(7)    ! spy - Define computational grid in y direction
          SpanY(1) = GetReal()
          SpanY(2) = GetReal()

       case(8)    ! spz - Define computational grid in z direction
          SpanZ(1) = GetReal()
          SpanZ(2) = GetReal()

       case(9)    ! dcel - Define computational grid space
          DCell = GetReal()

          if (abs(DCell) .le. 1.0e-15)  &
               write(*,"(a50)") '***Warning*** Computational grid &
                                 is too small'

       case(10)   ! DTScale - Define scale for time step size
          DTScale = GetReal()
          write(*,"(a,f10.3)") 'DTScale = ',DTScale
          write(iomsg,"(a,f10.3)") 'DTScale = ',DTScale
          if (DTScale .gt. 2.5) then 
             call ErrorMsg()
             stop '*** Error *** DTScale too large'
          end if

       case(11)   ! outTime - Define time interval 
                  !           for writting plot file
          OutTime = GetReal()

       case(12)   ! ReportTime - Define time interval 
                  ! for reporting status
          ReportTime = GetReal()

       case(13)   ! fixe - Define fixed rigid plane
          do i = 1, 6
             FixS(i) = GetInt()
          end do

       case(14)   ! nmat - Number of material sets
          nb_mat = GetInt()
          allocate(mat_list(nb_mat))
          mat_list%EosType = 0 ! no EOS as default
          write(*,"(a,i12)") 'Number of material sets  = ',nb_mat
          write(iomsg,"(a,i12)") 'Number of material sets  = ',nb_mat

       case(15)   ! material - Read in material properties
          call SetMaterial()

       case(16)   ! particle - Read in particles
          call SetParticle()

       case(17)   ! MUSL
          MUSL = SetOnOff()

          if(MUSL) then
             write(*,"(a)") ' MUSL method is used '
             write(iomsg,"(a)") 'MUSL method is used '
             updateMethod = .true.
          endif

       case(18)   ! Jaum - Use Jaumann rate ?
          Jaum = SetOnOff()

       case(19)   ! load - Define external load
          call SetLoad()

       case(20)   ! velo - Read in initial velocity
          call SetVelocity()

       case(21)   ! outres - Define variables which 
                  !          will be saved to plot file
          nAnimate = nAnimate + 1
          AnimOption(nAnimate) = SetResOption()

       case(22)   ! curve - Define variables which 
                  !         will be saved to plot file
          nCurves = nCurves + 1
          CurveOption(nCurves) = SetResOption()
          if(nb_word.gt.0) then
             tempI = GetInt()
             if (tempI .le. 0 .or. tempI .gt. nb_particle) then
                call ErrorMsg()
                stop '*** Error *** Invalid curved particle ID!'
             else
                CurvePoint(nCurves) = tempI
             end if
          else    ! No curve point specified
             CurvePoint(nCurves) = 1    
          end if

       case(23)   ! seos - Define equation of state
          call SetEos()

       case(24)   ! pt2d
          plot2dTrue = .true.
          plot2d(1) = GetReal()
          plot2d(2) = GetReal()
          plot2d(3) = GetReal()
          plot2d(4) = GetReal()
          plot2d(5) = GetReal()
          plot2d(6) = GetReal()

       case(25)   ! curx - curve x y z
          call SetCurX()

       case(26)   ! deto - Read in detonation point
          nDeto = nDeto + 1
          DetoX(nDeto) = GetReal()
          DetoY(nDeto) = GetReal()
          DetoZ(nDeto) = GetReal()

       case(27)   ! tecp - Write results to TecPlot
          WriteTecPlot = .true.

       case(28)   ! bulk - Set artificial bulk viscosity
          bq1 = GetReal()
          bq2 = GetReal()

       case(29)   ! GIMP 
          GIMP = .true.
          write(*,"(a)") 'GIMP is used'
          write(iomsg,"(a)") 'GIMP is used'

       case(30)   ! contact ( = 0 by default )
          contact = .true.   !using contact method
          if(contact)then
             call Setcontact
          end if

       case(31)
          USF = SetOnOff()   !using USF method    
          if(USF) then
             write(*,"(a)") ' USF method is used '
             write(iomsg,"(a)") 'USF method is used '
             updateMethod = .true.
       end if
   
       case(32)
          USL = SetOnOff()   !using USL method    
          if(USL) then
             write(*,"(a)") ' USL method is used '
             write(iomsg,"(a)") 'USL method is used '
             updateMethod = .true.
          end if

       case(33)    !nbco - number of components
          nb_component = GetInt() 
          write(*,"(a,i12)") 'Number of components = ',nb_component
          write(iomsg,"(a,i12)") 'Number of components = ',nb_component

       case(34)    ! nbbody - number of bodies
          nb_body = GetInt()
          write(*,"(a,i12)") 'Number of bodies = ',nb_body
          write(iomsg,"(a,i12)") 'Number of bodies = ',nb_body
          allocate(body_list(nb_body))
		  call InitBody()
           
       case(35)   ! para - Write results to ParaView
          WriteParaView = .true.

       case(36)   ! drda - Activate damping for static solution using dynamic relaxation
          DR_DAMPING = .true.
          
       case(37)   ! BSMPM
          Bspline = .true.
          write(*,"(a)") 'BSMPM is used'
          write(iomsg,"(a)") 'BSMPM is used'
          
       case(38)   ! SGMP
          SGMP = .true.
          write(*,"(a)") 'SGMP is used'
          write(iomsg,"(a)") 'SGMP is used'
          
       case(39)   ! SmoothStressbyGrid
           SmoothStress = .true.
           write(*,"(a)") 'SmoothStressbyGrid is used'
           write(iomsg,"(a)") 'SmoothStressbyGrid is used'
           
       case(40)   ! TLMPM
           TLMPM = .true.
           write(*,"(a)") 'TLMPM is used'
           write(iomsg,"(a)") 'TLMPM is used'
       case(41)   ! UpdateVbyFp
            UpdateVbyFp = .true.
            write(*,"(a)") 'UpdateVolumebyFp is used'
            write(iomsg,"(a)") 'UpdateVolumebyFp is used'
        case(42)   !  'QuasiLoad'
            QuasiLoad = .true.
            QuasiLoadTime = GetReal()
            write(*,"(a)") 'QuasiLoad is on'
            write(iomsg,"(a)") 'QuasiLoad is on'
        case(43)   ! ,'QuasiDamp'
            QuasiDamp = .true.
            DampFactor = GetReal()
            write(*,"(a)") 'QuasiDamp is on'
            write(iomsg,"(a)") 'QuasiDamp is on'
        
       case default ! error
          stop 'STOP - Error encountered in reading data'

       end select
    end do

    ! Set MUSL as the default option
    if(.NOT. updateMethod ) MUSL = .true. 

  end subroutine GetData

  subroutine SetMaterial()
!------------------------------------------------------------------
!-    purpose: Get material data using FFI module                 -
!-             data are stored in mat_list                        -
!------------------------------------------------------------------
    use FFI
    implicit none

    integer i, t, mtype, matcount

    integer,parameter:: nbkw = 11
    character(4),parameter:: kw(nbkw) = (/&
         'elas','pla1','pla2','john','sjc ',&
         'sjcf','jcf ','hiex','null','dpm ',&
         'neoh'/)

    if(nb_mat.eq.0 .or. nb_particle.eq.0) then
       stop '*** Error *** nb_mat/nb_particle must be defined in advance!'
    end if

    write(iomsg,"(a)") 'Material:'
    i = 0
    do while(i.lt.nb_mat)
       i = GetInt()
       write(*,*) trim(sss)
       write(iomsg,"(a)") trim(sss)

       if(i.gt.nb_mat) then
          call ErrorMsg()
          print *, '*** Error *** Too many material sets'
          print *, 'required : ',nb_mat
          print *, 'inputed  : ',i
          stop 
       end if

       t = KeyWord(kw,nbkw)
       select case(t)

       case(1)    ! elastic
          mat_list(i)%MatType = 1
          mat_list(i)%Density = GetReal()
          mat_list(i)%Young = GetReal()
          mat_list(i)%Poisson = GetReal()

       case(2) ! pla1: elastic-perfectly plastic
          mat_list(i)%MatType = 2
          mat_list(i)%Density = GetReal()
          mat_list(i)%Young = GetReal()
          mat_list(i)%Poisson = GetReal()
          mat_list(i)%Yield0 = GetReal()

       case(3) ! pla2: isotropic hardening
          mat_list(i)%MatType = 3
          mat_list(i)%Density = GetReal()
          mat_list(i)%Young = GetReal()
          mat_list(i)%Poisson = GetReal()
          mat_list(i)%Yield0 = GetReal()
          mat_list(i)%TangMod = GetReal()

       case(4) ! johnson-cook
          mat_list(i)%MatType = 4
          mat_list(i)%Density = GetReal()
          mat_list(i)%Young = GetReal()
          mat_list(i)%Poisson = GetReal()
          mat_list(i)%Yield0 = GetReal()
          mat_list(i)%B_jc = GetReal()
          mat_list(i)%n_jc = GetReal()
          mat_list(i)%C_jc = GetReal()
          mat_list(i)%m_jc = GetReal()
          mat_list(i)%roomt = GetReal()
          mat_list(i)%melt = GetReal()
          mat_list(i)%SpecHeat = GetReal()
          mat_list(i)%epso = GetReal()

       case(5) ! Simplified johnson-cook
          mat_list(i)%MatType = 5
          mat_list(i)%Density = GetReal()
          mat_list(i)%Young = GetReal()
          mat_list(i)%Poisson = GetReal()
          mat_list(i)%Yield0 = GetReal()
          mat_list(i)%B_jc = GetReal()
          mat_list(i)%n_jc = GetReal()
          mat_list(i)%C_jc = GetReal()
          mat_list(i)%epso = GetReal()

       case(6) ! Simplified johnson-cook with failure
          mat_list(i)%MatType = 6
          mat_list(i)%Density = GetReal()
          mat_list(i)%Young = GetReal()
          mat_list(i)%Poisson = GetReal()
          mat_list(i)%Yield0 = GetReal()
          mat_list(i)%B_jc = GetReal()
          mat_list(i)%n_jc = GetReal()
          mat_list(i)%C_jc = GetReal()
          mat_list(i)%epso = GetReal()
          mat_list(i)%epf = GetReal()

       case(7) ! johnson-cook with fracture
          mat_list(i)%MatType = 7
          mat_list(i)%Density = GetReal()
          mat_list(i)%Young = GetReal()
          mat_list(i)%Poisson = GetReal()
          mat_list(i)%Yield0 = GetReal()
          mat_list(i)%B_jc = GetReal()
          mat_list(i)%n_jc = GetReal()
          mat_list(i)%C_jc = GetReal()
          mat_list(i)%m_jc = GetReal()
          mat_list(i)%roomt = GetReal()
          mat_list(i)%melt = GetReal()
          mat_list(i)%SpecHeat = GetReal()
          mat_list(i)%epso = GetReal()    
          mat_list(i)%epf = GetReal()

       case(8) ! High explosive burn
          mat_list(i)%MatType = 8
          mat_list(i)%Density = GetReal()
          mat_list(i)%D = GetReal()

       case(9) ! null(used to model water)
          mat_list(i)%MatType = 9
          mat_list(i)%Density = GetReal()
          mat_list(i)%BulkModulus = GetReal()
          mat_list(i)%miu     = GetReal()
          mat_list(i)%Wavespd = SQRT(mat_list(i)%BulkModulus/mat_list(i)%Density)

       case(10) ! drucker-prager model
          mat_list(i)%MatType = 10
          mat_list(i)%Density = GetReal()
          mat_list(i)%Young =   GetReal()
          mat_list(i)%Poisson = GetReal()
          mat_list(i)%q_fai = GetReal()
          mat_list(i)%k_fai = GetReal()
          mat_list(i)%q_psi = GetReal()
          mat_list(i)%ten_f = GetReal()    
       case(11)    ! NeoHookean
          mat_list(i)%MatType = 11
          mat_list(i)%Density = GetReal()
          mat_list(i)%Young = GetReal()
          mat_list(i)%Poisson = GetReal()
          UpdateVbyFp = .true.

       case default
          call ErrorMsg()
          stop '*** Error *** Invalid Material Type!'

       end select

    end do

    ! set CutOff value
    CutOff = 0.0
    matcount =0
    do i = 1, nb_mat
       CutOff = mat_list(i)%density + CutOff
       matcount = matcount + 1
    end do
    CutOff = CutOff * DCell**3 * 1e-5 / matcount

    ! print *, 'Grid mass CutOff value: ', CutOff
    write(*,"(a,i12)") 'the effective number of material',matcount
    write(iomsg,"(a,i12)") 'the effective number of material',matcount
    write(*,"(a,e12.4)") 'Grid mass CutOff value: ', CutOff
    write(iomsg,"(a,e12.4)") 'Grid mass CutOff value: ', CutOff

  end subroutine SetMaterial

  subroutine SetParticle()
!------------------------------------------------------------------
!-  purpose: Get particle data using FFI module                   -
!-           data are stored in particle_list                     -
!------------------------------------------------------------------
    use FFI
    implicit none

    integer p, newi, PtOption, pnum
    integer matID
    integer ix, iy, iz
    integer nx, ny, nz,nr
    real(8) ox, oy, oz ! coordinates of origin point
    real(8) oox,ooy,ooz
    real(8) tx, ty, tz
    real(8) dp, radius2, pmass

    integer,parameter:: nbkw = 5
    character(4),parameter:: kw(nbkw) = (/  &
         'poin','bloc','sphe','blod','cicle'  /)

    if(nb_particle.eq.0) then
       stop '*** Error *** nbmp must be defined in advance!'
    end if

    PtOption = keyword(kw,nbkw)

    select case(PtOption)

    case(1)    ! point (comID, pnum, pID, matID, pmass, x, y, z)
       bodyCounter = bodyCounter + 1
       body_list(bodyCounter)%comID = GetInt()

       pnum = GetInt()
       body_list(bodyCounter)%par_begin = parCounter + 1
       body_list(bodyCounter)%par_end = parCounter + pnum
       do p = 1, pnum
          newi = GetInt()
          parCounter = parCounter + 1    ! count the particle number
          if (newi .lt. parCounter) then
             call ErrorMsg()
             stop 'particle ID conflict'
          end if

          if (parCounter.gt.nb_particle) then
             call ErrorMsg()
             print *, '*** Error *** Too many particles'
             print *, 'required : ',nb_particle
             print *, 'inputed  : ',parCounter
             stop 
          end if

          matID = GetInt()

          if (matID .gt. nb_mat .or. matID .le. 0) then
             call ErrorMsg()
             stop '*** Error *** Improper material ID'
          end if

          particle_list(parCounter)%mass = GetReal()

          if(matID.gt.nb_mat) then
             call ErrorMsg()
             stop '*** Error *** material ID greater than nb_mat'
          end if

          particle_list(parCounter)%Xp(1) = GetReal()
          particle_list(parCounter)%Xp(2) = GetReal()
          particle_list(parCounter)%Xp(3) = GetReal()

          if(plot2dTrue) then
             call setskip(parCounter)
          end if
       end do
       body_list(bodyCounter)%mat = matID

    case(2)    ! block (comID, matID, density, dp, ox, oy, oz, nx, ny, nz)
       bodyCounter = bodyCounter + 1
       body_list(bodyCounter)%comID = GetInt()

       matID = GetInt()
       body_list(bodyCounter)%mat = matID
       body_list(bodyCounter)%par_begin = parCounter +1

       if(matID.gt.nb_mat .or. matID.le.0) then
          call ErrorMsg()
          stop '*** Error *** improper material ID'
       end if

       pmass = GetReal()
       dp = GetReal()    ! distance between particles
	   pmass = pmass*dp**3
       ox = GetReal()    ! origin of the block
       oy = GetReal() 
       oz = GetReal()
       !repair
       ox = ox + dp*0.5  ! offset the origin to the first particle
       oy = oy + dp*0.5
       oz = oz + dp*0.5
       nx = GetInt()     ! number of particle in each direction
       ny = GetInt()
       nz = GetInt()

       do iz = 1, nz
          do iy = 1, ny
             do ix = 1, nx
                ! count the particle number			    
                parCounter = parCounter + 1 
                if(parCounter.gt.nb_particle) then
                   call ErrorMsg()
                   call ErrorPt()
                end if

                particle_list(parCounter)%mass = pmass

                if(matID.gt.nb_mat) then
                   call ErrorMsg()
                   stop '*** Error *** material ID greater than nb_mat'
                end if

                particle_list(parCounter)%Xp(1) = (ix - 1)*dp + ox
                particle_list(parCounter)%Xp(2) = (iy - 1)*dp + oy
                particle_list(parCounter)%Xp(3) = (iz - 1)*dp + oz

                if(plot2dTrue) then
                   call setskip(parCounter)
                end if

             end do
          end do
       end do

       body_list(bodyCounter)%par_end = parCounter

    case(3)  ! sphere (comID, matID, density, dp, ox, oy, oz, nRadius)
       bodyCounter = bodyCounter + 1
       body_list(bodyCounter)%comID = GetInt()

       matID = GetInt()
       body_list(bodyCounter)%mat = matID
       body_list(bodyCounter)%par_begin = parCounter +1
       if(matID.gt.nb_mat .or. matID.le.0) then
          call ErrorMsg()
          stop '*** Error *** improper material ID'
       end if

       pmass = GetReal()
       dp = GetReal()    ! distance between particles
	   pmass = pmass*dp**3
       ox = GetReal()    ! origin of the sphere
       oy = GetReal() 
       oz = GetReal()
       nx = GetReal()
       radius2 = (nx * dp ) ** 2
       do ix = 1-nx, nx
          do iy = 1-nx, nx
             do iz = 1-nx, nx
                tx = (ix - 1)*dp + ox + 0.5*dp
                ty = (iy - 1)*dp + oy + 0.5*dp
                tz = (iz - 1)*dp + oz + 0.5*dp
                if (((tx - ox)**2 + (ty - oy)**2 + &
                     (tz - oz)**2) .lt. radius2) then
                   ! count the particle number
                   parCounter = parCounter + 1    
                   if(parCounter.gt.nb_particle) then
                      call ErrorMsg()
                      call ErrorPt()
                   end if

                   particle_list(parCounter)%mass = pmass

                   if(matID.gt.nb_mat) then
                      call ErrorMsg()
                      stop '*** Error *** material ID greater than nb_mat'
                   end if

                   particle_list(parCounter)%Xp(1) = tx
                   particle_list(parCounter)%Xp(2) = ty
                   particle_list(parCounter)%Xp(3) = tz

                   if(plot2dTrue) then
                      call setskip(parCounter)
                   end if

                end if
             end do
          end do
       end do

       body_list(bodyCounter)%par_end = parCounter

    case(4)    ! blod (comID, matID, density, dp, ox, oy, oz, nx, ny, nz,oox,ooy,ooz,nr)
       bodyCounter = bodyCounter + 1
       body_list(bodyCounter)%comID = GetInt()

       matID = GetInt()
       body_list(bodyCounter)%mat = matID
       body_list(bodyCounter)%par_begin = parCounter +1

       if(matID.gt.nb_mat .or. matID.le.0) then
          call ErrorMsg()
          stop '*** Error *** improper material ID'
       end if

       pmass = GetReal()
       dp = GetReal()    ! distance between particles
	   pmass = pmass*dp**3
       ox = GetReal()    ! origin of the block
       oy = GetReal() 
       oz = GetReal()
       !repair
       ox = ox + dp*0.5  ! offset the origin to the first particle
       oy = oy + dp*0.5
       oz = oz + dp*0.5
       nx = GetInt()     ! number of particle in each direction
       ny = GetInt()
       nz = GetInt()
       oox=Getreal()
       ooy=Getreal()
       ooz=Getreal()
       nr=Getint()
       radius2 = (nr * dp ) ** 2

       do iz = 1, nz
          do iy = 1, ny
             do ix = 1, nx
                 tx=(ix - 1)*dp + ox
                 ty=(iy - 1)*dp + oy
                 tz=(iz - 1)*dp + oz
                 if (((tx - oox)**2 + (ty - ooy)**2).GE. radius2) then
                   ! count the particle number
                ! count the particle number			    
                parCounter = parCounter + 1 
                if(parCounter.gt.nb_particle) then
                   call ErrorMsg()
                   call ErrorPt()
                end if

                particle_list(parCounter)%mass = pmass

                if(matID.gt.nb_mat) then
                   call ErrorMsg()
                   stop '*** Error *** material ID greater than nb_mat'
                end if

                particle_list(parCounter)%Xp(1) = tx
                particle_list(parCounter)%Xp(2) = ty
                particle_list(parCounter)%Xp(3) = tz

                if(plot2dTrue) then
                   call setskip(parCounter)
                end if
                end if

             end do
          end do
       end do

       body_list(bodyCounter)%par_end = parCounter
       
    case(5)  ! plate (comID, matID, density, dp, ox, oy, oz, nRadius)
       bodyCounter = bodyCounter + 1
       body_list(bodyCounter)%comID = GetInt()

       matID = GetInt()
       body_list(bodyCounter)%mat = matID
       body_list(bodyCounter)%par_begin = parCounter +1
       if(matID.gt.nb_mat .or. matID.le.0) then
          call ErrorMsg()
          stop '*** Error *** improper material ID'
       end if

       pmass = GetReal()
       dp = GetReal()    ! distance between particles
	   pmass = pmass*dp**3
       ox = GetReal()    ! origin of the sphere
       oy = GetReal() 
       oz = GetReal()
       nx = GetInt()
       radius2 = (nx * dp ) ** 2
       do ix = 1-nx, nx
          do iy = 1-nx, nx
                tx = (ix - 1)*dp + ox + 0.5*dp
                ty = (iy - 1)*dp + oy + 0.5*dp
                tz =  oz + dp*0.5
                if (((tx - ox)**2 + (ty - oy)**2 ) .lt. radius2) then
                   ! count the particle number
                   parCounter = parCounter + 1    
                   if(parCounter.gt.nb_particle) then
                      call ErrorMsg()
                      call ErrorPt()
                   end if

                   particle_list(parCounter)%mass = pmass

                   if(matID.gt.nb_mat) then
                      call ErrorMsg()
                      stop '*** Error *** material ID greater than nb_mat'
                   end if

                   particle_list(parCounter)%Xp(1) = tx
                   particle_list(parCounter)%Xp(2) = ty
                   particle_list(parCounter)%Xp(3) = tz

                   if(plot2dTrue) then
                      call setskip(parCounter)
                   end if

                end if
             end do
          end do

       body_list(bodyCounter)%par_end = parCounter

    case default    ! error
       call ErrorMsg()
       stop 'particle option should be point/block/sphere'

    end select

    if(bodyCounter == nb_body) then
        write(*,"(a,i12)") 'particles read:', parCounter
    end if 

    ! set particle initial position
    do p = 1, nb_particle
      particle_list(p)%Xo(1) = particle_list(p)%Xp(1)
      particle_list(p)%Xo(2) = particle_list(p)%Xp(2)
      particle_list(p)%Xo(3) = particle_list(p)%Xp(3)
    end do

  end subroutine SetParticle

  subroutine setskip(p)
    implicit none
    integer:: p

    if (particle_list(p)%Xp(1)<plot2d(1) .or. &
        particle_list(p)%Xp(1)>plot2d(2)) then 
       particle_list(p)%SkipThis = .true.
    else if (particle_list(p)%Xp(2)<plot2d(3) .or. &
             particle_list(p)%Xp(2)>plot2d(4)) then 
       particle_list(p)%SkipThis = .true.
    else if (particle_list(p)%Xp(3)<plot2d(5) .or. &
             particle_list(p)%Xp(3)>plot2d(6)) then 
       particle_list(p)%SkipThis = .true.
    end if

  end subroutine setskip

  subroutine SetLoad()
!------------------------------------------------------------------
!-   purpose: Get Load by component or by node                    -
!------------------------------------------------------------------
    use FFI
    implicit none

    integer k, ipart, i, ibody
    real(8):: fxp, fyp, fzp
    integer,parameter:: nbkw = 5
    character(4),parameter:: kw(nbkw) = &
                             (/'endl','part','body','grav','node' /)

    if(nb_body.eq.0) then
       stop '*** Error *** nbby must be defined !'
    end if

    do while(.true.)
       
       k = keyword(kw,nbkw)
       select case(k)
       
       case(1)    ! end load
          exit

       case(2)    ! by particle 
          ipart = GetInt()
          particle_list(ipart)%FXp(1) = GetReal()
          particle_list(ipart)%FXp(2) = GetReal()
          particle_list(ipart)%FXp(3) = GetReal()
          
       case(3)    ! by body
          ibody = GetInt()    ! body number
          fxp = GetReal()
          fyp = GetReal()
          fzp = GetReal()

          do i = body_list(ibody)%par_begin, body_list(ibody)%par_end
             particle_list(i)%FXp(1) = fxp
             particle_list(i)%FXp(2) = fyp
             particle_list(i)%FXp(3) = fzp
          end do

       case(4) ! by body
	      Gravity = .true.
          ibody = GetInt()    ! body number
          body_list(ibody)%Gravp(1) = GetReal()
          body_list(ibody)%Gravp(2) = GetReal()
          body_list(ibody)%Gravp(3) = GetReal()

       case (5) ! Legacy keyword
         write(*,*) '*** Error *** : Please use the keyword part instead of node to apply a load to a particle'
       
       case default    ! error
          call ErrorMsg()
          stop 'error GetLoad'

       end select
    end do

  end subroutine SetLoad

  subroutine SetVelocity()
!------------------------------------------------------------------
!-  purpose: Get Initial Velocity by component or by node         -
!------------------------------------------------------------------
    use FFI
    implicit none

    integer k, inode, ibody, cpl, i,j,m,n,parBegin,parEnd,b,p
    real(8):: vxp, vyp, vzp,v0,dc,y,A,W,q,c1,c2,x
    integer,parameter:: nbkw = 6
    character(4),parameter:: kw(nbkw) = (/'endv','node','body','ksin','kbonuli','vsin'/)
    type(Particle),POINTER::pt

    if(nb_body.eq.0) then
       stop '*** Error *** nbby must be defined !'
    end if

    do while(.true.)
       k = keyword(kw,nbkw)
       select case(k)
       case(1)    ! endload
          exit

       case(2)    ! by node
          inode = GetInt()
          particle_list(inode)%VXp(1) = GetReal()
          particle_list(inode)%VXp(2) = GetReal()
          particle_list(inode)%VXp(3) = GetReal()

       case(3) ! by body
          ibody = GetInt()    ! body number
          vxp = GetReal()
          vyp = GetReal()
          vzp = GetReal()

          do i = body_list(ibody)%par_begin, body_list(ibody)%par_end
             particle_list(i)%VXp(1) = vxp
             particle_list(i)%VXp(2) = vyp
             particle_list(i)%VXp(3) = vzp
          end do
       case(4)!bysin
          v0=GetReal()
          ibody=GetInt()
          n=320
          dc=0.003125
          if(ibody==1)then
        do j=1,4
          do i=n*(j-1)+1,n*j
            m=mod(i,n)
            if(m==0)m=n
            particle_list(i)%VXp(1) = v0*sin(3.1415926535*((m-1)*dc+dc/2.0))
            particle_list(i)%VXp(2) = 0
            particle_list(i)%VXp(3) = 0
          end do
        end do
          end if
          case(5)!EUler beam
          y=GetReal()
          n=384
          dc=0.00015625
          q=78.83
          c1=cosh(q*0.06)-cos(q*0.06)
          c2=sinh(q*0.06)-sin(q*0.06)
          A=y/(2*c2)
          W=q*q*sqrt(318e9*0.01*0.01/(12*1845*(1-0.054*0.054)))
        do j=1,128
          do i=n*(j-1)+1,n*j
            m=mod(i,n)
            if(m==0)m=n
            x=((m-1)*dc+dc/2.0)
            particle_list(i)%VXp(1) = 0
            particle_list(i)%VXp(2) = A*W*(c1*(sinh(q*x)+sin(q*x))-c2*(cosh(q*x)+cos(q*x)))
            particle_list(i)%VXp(3) = 0
          end do
        end do
        
        case(6)!used for Free softbeam
        v0 = GetReal()
            do b = 1,nb_body
                parBegin = body_list(b)%par_begin
                parEnd   = body_list(b)%par_End
                do p = parBegin,parEnd
                    pt => particle_list(p)
                    pt%VXp(1)=v0*sin(3.1415926535*pt%xo(2)/6)
                    pt%VXp(2)=0
                    pt%VXp(3)=0
                end do
            end do
        
       case default    ! error
          call ErrorMsg()
          stop 'error GetVelocity'

       end select
    end do

  end subroutine SetVelocity

  logical function SetOnOff()
!------------------------------------------------------------------
!-   purpose: Get On/Off switch                                   -
!------------------------------------------------------------------
    use FFI
    implicit none

    integer:: k
    integer,parameter:: nbkw = 2
    character(4),parameter:: kw(nbkw) = (/'on ','off'/)

    k = keyword(kw,nbkw)

    select case(k)
    case(1)    ! on
       SetOnOff = .true.

    case(2)    ! off
       SetOnOff = .false.

    case default    ! error
       call ErrorMsg()
       stop 'switch should be ON/OFF'

    end select

  end function SetOnOff

  integer function SetResOption()
!------------------------------------------------------------------
!-    purpose: Get result option                                  -
!------------------------------------------------------------------
    use FFI
    implicit none
    integer:: k
    integer,parameter:: nbkw = 18
    character(4),parameter:: kw(nbkw) = (/ &
         'seqv','epef','mat ','pres','volu',&
         'engk','engi','velx','vely','velz',&
         'cels','fail','sspd','damg','stressx','disx',&
         'disy','disz'/)

    SetResOption = keyword(kw,nbkw)

    if(SetResOption.le.0 .or. SetResOption.gt.nbkw) then
       call ErrorMsg()
       stop 'Result option error'
    end if

  end function SetResOption

  subroutine SetEos()
!------------------------------------------------------------------
!-        purpose: set EOS parameter                              -
!------------------------------------------------------------------
    use FFI
    implicit none
    integer mset
    real(8):: c0, lambda, gamma0

    mset = GetInt()
    ! tie a EosType with given Material set
    mat_list(mset)%EosType = GetInt() 

    select case(mat_list(mset)%EosType)

    case(1) ! Linear polynomial EOS

       mat_list(mset)%cEos(1) = GetReal()    ! c0
       mat_list(mset)%cEos(2) = GetReal()    ! c1
       mat_list(mset)%cEos(3) = GetReal()    ! c2
       mat_list(mset)%cEos(4) = GetReal()    ! c3
       mat_list(mset)%cEos(5) = GetReal()    ! c4
       mat_list(mset)%cEos(6) = GetReal()    ! c5
       mat_list(mset)%cEos(7) = GetReal()    ! c6
       mat_list(mset)%cEos(10) = GetReal()   ! E0

    case(2) ! Gruneisen EOS

       c0 = GetReal()
       lambda = GetReal()
       gamma0 = GetReal()
       ! C1
       mat_list(mset)%cEos(1) = mat_list(mset)%density * (c0**2)
       ! C2
       mat_list(mset)%cEos(2) = mat_list(mset)%cEos(1) * (2*lambda - 1)
       ! C3
       mat_list(mset)%cEos(3) = mat_list(mset)%cEos(1) * &
                                (lambda - 1) * (3*lambda - 1)
       ! C4
       mat_list(mset)%cEos(4) = mat_list(mset)%density * gamma0                
       mat_list(mset)%cEos(5) = gamma0            ! C5

       mat_list(mset)%cEos(10) = GetReal()        ! E0

    case(3) ! JWL EOS
       mat_list(mset)%cEos(1) = GetReal()        ! A
       mat_list(mset)%cEos(2) = GetReal()        ! B
       mat_list(mset)%cEos(3) = GetReal()        ! R1
       mat_list(mset)%cEos(4) = GetReal()        ! R2
       mat_list(mset)%cEos(5) = GetReal()        ! w
       mat_list(mset)%cEos(10) = GetReal()       ! E0
    
    case(4) !Weakly Compressible EOS
      mat_list(mset)%cEos(1) = GetReal()         !C Numerical Speed

    case default
       call ErrorMsg()
       stop 'Error set EOS parameter (seos)'

    end select

  end subroutine SetEos

  subroutine SetCurX()
!-------------------------------------------------------------------
!-  purpose:                                                       -
!-     find particle ID according to coordinates                   -
!-------------------------------------------------------------------
    use FFI
    implicit none

    integer:: pid, bid
    real(8):: tx, ty, tz, maxx

    if (parCounter .ne. nb_particle) then
       call ErrorMsg()
       stop 'curx should be used after particle input'
    end if

    nCurves = nCurves + 1
    CurveOption(nCurves) = SetResOption()

    tx = GetReal()
    ty = GetReal()
    tz = GetReal()

    do bid = 1, nb_body
       do pid = body_list(bid)%par_begin, body_list(bid)%par_end
          maxx = MAX(ABS(particle_list(pid)%Xp(1)-tx), &
                     ABS(particle_list(pid)%Xp(2)-ty),  &
                     ABS(particle_list(pid)%Xp(3)-tz))

          if (maxx .le. (0.25*DCell)) then
             CurvePoint(nCurves) = pid
             exit
          end if
       end do
    end do
    
  end subroutine SetCurX


  subroutine ErrorPt()
!-------------------------------------------------------------------
!-                                                                 -
!-  purpose:                                                       -
!-     output error message and line number                        -
!-                                                                 -
!-------------------------------------------------------------------
    implicit none
    print *, '*** Error *** Too many particles'
    print *, 'required : ',nb_particle
    print *, 'inputed  : ',parCounter
    stop 
  end subroutine ErrorPt

  subroutine Initial()
!------------------------------------------------------------------
!-  purpose: Initialize specific variable used after reading all  -
!-           particle data                                        -
!-    including: VOL, sig_y, LT, ie                               -
!------------------------------------------------------------------
    use ParticleData
    use MaterialData,only:nb_mat
    use FFI
    implicit none
    integer i, j, m, c, mtype, b
    integer parBegin, parEnd
    real temp
    real(8)  nu, E, ro
    type(particle), pointer:: pt
    type(body), pointer:: bd
    

    do b = 1, nb_body
       bd => body_list(b)
       parBegin = bd%par_begin
       parEnd = bd%par_end
       do i = parBegin, parEnd
          pt => particle_list(i)
          m = bd%mat
          mtype = mat_list(m)%MatType
          !The deformable body
          if (mtype .ne.12)then         
             pt%VOL = pt%mass / mat_list(m)%Density
             pt%VOL0 =  pt%VOL
             pt%sig_y = mat_list(m)%Yield0
             pt%ie = mat_list(m)%cEos(10) * pt%VOL
             if (mtype==1.OR.mtype==11 )then
                 E = mat_list(m)%Young
                 nu = mat_list(m)%Poisson
                 ro = mat_list(m)%Density
                 pt%cp = sqrt(E*(1-nu)/(1+nu)/(1-2*nu)/ro) 
             end if
             if (mtype==9 ) then 
                pt%cp = mat_list(m)%Wavespd 
                mat_list(m)%cEos(1)=mat_list(m)%Wavespd
             end if
             pt%LT = 10e6
             do j = 1, nDeto
                temp = sqrt((pt%Xp(1)-DetoX(j))**2 + &
                            (pt%Xp(2)-DetoY(j))**2 + &
                            (pt%Xp(3)-DetoZ(j))**2) /   &
                            mat_list(m)%D
                pt%LT = min(temp,pt%LT)
             end do
          end if
          
          
       end do
    end do

  end subroutine Initial

  subroutine SetDT()
!------------------------------------------------------------------
!-  purpose:  Computational time step size                        -
!------------------------------------------------------------------
    use ParticleData
    implicit none
    integer m, mtype
    real(8) cp, nu, E, ro

    DT = 1.0e6

    do m = 1, nb_mat

       mtype = mat_list(m)%MatType
       if (mtype .ne.12)then         !The deformable body

          E = mat_list(m)%Young
          nu = mat_list(m)%Poisson
          ro = mat_list(m)%Density
          cp = sqrt(E*(1-nu)/(1+nu)/(1-2*nu)/ro)    ! sound speed
          cp = max(cp,mat_list(m)%D)
          cp = max(cp,mat_list(m)%Wavespd)
          DT = min(DT, DCell/cp)
       end if
    end do

    DT = DT*DTscale

  end subroutine SetDT


  subroutine statinfo()
!------------------------------------------------------------------
!-   purpose: Output statistical data to message file             -
!------------------------------------------------------------------
    use ParticleData
    use FFI

    implicit none
    integer:: i, cpl, inode, ibody, parBegin, parEnd
    real:: cmass, tmass=0, cke, cie, tke=0, tie=0
    type(Particle), pointer :: pt

    ! tmass = 0
    ! tke = 0
    write(iomsg,*)
    do ibody = 1, nb_body
       write(iomsg,"(a13,i3,a4)") '--- body', ibody, ' ---'

       ! write mass and energy information
       cmass = 0
       cke = 0
       cie = 0

       parBegin = body_list(ibody)%par_begin
       parEnd = body_list(ibody)%par_end 

       do i = parBegin, parEnd
          pt => particle_list(i)
          cmass = cmass + pt%mass
          cke = cke + dot_product(pt%VXp, pt%VXp)*pt%mass*0.5d0
          cie = cie + pt%ie
       end do
       tmass = tmass + cmass
       tke = tke + cke
       tie = tie + cie
       write(iomsg,"(a25,e12.4)") ' mass:', cmass
       write(iomsg,"(a25,e12.4)") '  kinetic energy:', cke
       write(iomsg,"(a25,e12.4)") ' internal energy:', cie
    end do

    ! write total infomation
    write(iomsg,"(a20)") '--- total ---'
    write(iomsg,"(a25,e12.4)") 'total mass:', tmass
    write(iomsg,"(a25,e12.4)") ' total kinetic energy:', tke
    write(iomsg,"(a25,e12.4)") 'total internal energy:', tie

  end subroutine statinfo

  subroutine Setcontact()
!------------------------------------------------------------------
!-    purpose: Get contact data using FFI module                  -
!------------------------------------------------------------------
    use FFI
    implicit none

    integer k, t
    integer,parameter:: nbkw = 2
    character(4),parameter:: kw(nbkw) = (/'lagr  ','others'/)

    k = keyword(kw,nbkw)
    select case(k)
    case(1)
       contact_type= 1     ! using Bardenhagens method
       fricfa = GetReal()  ! the frictional coefficient
       ! the computational method of contact normal
       normbody = GetInt()    
       write(*,"(a)") 'Notice: Bardenhagens contact method &
                       is used with USF '
       write(iomsg,"(a)") 'Notice: Bardenhagens contact method &
                           is used with USF'
       write(*,"(a,e10.3)") 'frictional factor = ',fricfa
       write(iomsg,"(a,e10.3)") 'frictional factor = ',fricfa

    case(2)            
       !using another contact method such as penalty method 
       ! implementing by users 
       contact_type= 2    
       fricfa = GetReal() ! the frictional coefficient
       ! the computational method of contact normal
       normbody = GetInt()    
       write(*,"(a)") 'Notice: another contact method &
                       is used with USF '
       write(iomsg,"(a)") 'Notice: another contact method &
                           is used with USF'
       write(*,"(a,e10.3)") 'frictional factor = ',fricfa
       write(iomsg,"(a,e10.3)") 'frictional factor = ',fricfa

    case default
       call ErrorMsg()
       stop 'Error set contact parameter (lagr/pena)'

    end select

  end subroutine Setcontact

end module DataIn
