
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
! -  Main procedures                                               -
! -                                                                -
! ------------------------------------------------------------------
program MPM3D
  use ParticleData
  use FFI, only: iomsg, iow1, iow2, iow12, FileAnim
  use DataIn
  use DataOut

  implicit none

  real:: t_begin, t_end, t_cpu=0.0, TPZC
  real:: t_bg, t_ed, t_elapsed, t0, t1
  real:: t_sec, t_s
  integer:: t_min, t_m
  integer:: iplotstep = 1

  real:: plt = 0.0
  real:: prt = 0.0

  call cpu_time( t_bg )

  call InputPara()    ! Input data

  call calcEnergy()   ! Calculate kinetic energy

  write(iomsg,*)
  write(iomsg,"(a,e10.3)") 'DT = ', DT

  plt = plt + OutTime

  ! Write results in TecPlot format: initial step
  if (WriteTecPlot) call OutAnim()

  ! Write results in ParaView format: initial step
  if (WriteParaView) call OutAnimPV(iplotstep)
  iplotstep = iplotstep + 1

  call cpu_time( t_ed )
  print *, '** Time for preprocessing is', t_ed - t_bg, ' seconds'

  write(*,"(a,e10.3)") 'DT = ', DT
  write(*,*) 'solving...'

  call cpu_time( t_bg )
  t0 = secnds(0.0)

  ! Solving
  do while(CurrentTime .le. EndTime)
     call cpu_time( t_begin )

     istep = istep+1
     CurrentTime = CurrentTime + DT
     EngInternal = 0.0

     ! Step 1: Initialize background grid nodal mass and Momentum
     call GridMomentumInitial()  ! Eq.(3.46) and Eq.(3.47)

     if(USF) then
	 ! Step 2: Apply boundary conditions 
	    call ApplyBoundaryConditions()
	 ! Step 3: Update particles stress 
        call ParticleStressUpdate() ! Eq.(3.49-3.50 ...)
     end if

     ! Step 4: Calculate the grid nodal force
     call GridMomentumUpdate() ! Eq.(3.52-3.54)

     ! Step 5: Integrate momentum equations on background grids
     call IntegrateMomentum()  ! Eq.(3.55)

     ! Step 6: Detect contact grid node, calculate contact force and
     !          adjust nodal momentum
     if(Contact_type == 1) then
        call Lagr_NodContact()
     end if

     ! Step 7: Update particles position and velocity
     call ParticlePositionUpdate() ! Eq.(3.56) and Eq.(3.57)

     ! Step 8: Recalculate the grid node momentum for MUSL
     if(MUSL) then
        call GridMomentumMUSL()    ! Eq.(3.58)
        call ApplyBoundaryConditions()
     end if

     ! Step 9: Update particles stress for both USL and MUSL
     if(.NOT. USF) then
        call ParticleStressUpdate() ! Eq.(3.60-3.62 ...)
     end if

     call calcEnergy()              ! Calculate kinetic energy

     call cpu_time( t_end )
     t_cpu = t_cpu + t_end - t_begin

     call OutCurve()    ! out put curve and animation data

     if (CurrentTime.ge.plt) then
        plt = plt + OutTime
        write(*,*) 'Write output data'
        write(iomsg,*) 'Write output data'

        ! Write results in TecPlot format: current step
	if (WriteTecPlot) call OutAnim()
        
	! Write results in ParaView format: current step
	if (WriteParaView) call OutAnimPV(iplotstep)
	iplotstep = iplotstep + 1
     end if

     ! report current computational progress
     if (CurrentTime.ge.prt) then
        prt = prt + ReportTime
        write(*,100) istep, CurrentTime, EngKinetic, &
                     EngKinetic+EngInternal
        write(iomsg,100) istep, CurrentTime, EngKinetic, &
                         EngKinetic+EngInternal
100     format(1x, "Step=", i6, 1p, "  T=", e10.3, "  &
                K.E.=", e10.3, "  T.E.=", e10.3)

        write(*,300) Momentum(1),Momentum(2),Momentum(3)
        write(iomsg,300) Momentum(1),Momentum(2),Momentum(3)
300     format(14x,"Mx=", e10.3, "   My=", e10.3, "   Mz=", e10.3)
     end if

  end do

  call cpu_time( t_ed )
  t1 = secnds(t0)

  t_elapsed = t_ed - t_bg

  t_min = floor(t_cpu/60)
  t_sec = t_cpu - t_min*60
  t_m = floor(t_elapsed/60)
  t_s = t_elapsed - t_m*60

  TPZC = t_cpu / nb_particle / istep

  write(*,200) t_cpu, t_min, t_sec
  write(*,201) t_elapsed, t_m, t_s

  write(iomsg,*)
  write(iomsg,200) t_cpu, t_min, t_sec

  write(iomsg,201) t_elapsed, t_m, t_s

  write(*,202) TPZC *1.0e9

  write(iomsg,202) TPZC*1.0e9

  write(*,203) t1

  write(iomsg,203) t1

  ! Close OutPut Files
  if (WriteTecPlot) close(iow1)  
  
  ! Close XML HEADER TAGS
  if (WriteParaview) then
	write(iow12,400)
	write(iow12,401)
	close(iow12)
  end if

  close(iow2)
  close(iomsg)

200 format("** Total CPU Time is ", &
           f12.4, " seconds (", i5, " minutes", f8.4, " seconds)")
201 format("** Elapsed time is   ", &
           f12.4, " seconds (", i5, " minutes", f8.4, " seconds)")
202 format("** Time per particle cycle: ", f12.4, " nanoseconds")
203 format("** Total elasped time : ", f12.4, " seconds")
400 format('    </Collection>')
401 format('</VTKFile>')

end program MPM3D
