
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
! -  Calculation procedures                                        -
! -                                                                -
! ------------------------------------------------------------------

subroutine GridMomentumInitial()
!-------------------------------------------------------------------
!-  Purpose                                                        -
!-      1. The variables of particle is mapped to the grid node    -
!-------------------------------------------------------------------

  use ParticleData
  use GridData
  use MaterialData
  implicit none

  integer:: b, p, n, c, parBegin, parEnd ! loop counter
  integer:: icell, inode, ix, iy, iz, mat_, comID = 1
  real(8):: sxx, syy, szz, sxy, syz, sxz
  real(8):: fx(3), f_int(3), f_ext(3), mp_, vol_
  real(8):: shm, SHPn, DNDXn, DNDYn, DNDZn

  type(Particle), POINTER :: pt
  type(GridNodeProperty), POINTER :: gd

  ! Calculate the grid nodal masses, moemntum only 
  ! Reset Grid data
  grid_list%Mg = 0.0d0         ! Grid nodal mass

  grid_list%PXg(1) = 0.0d0;    ! Nodal momentum
  grid_list%PXg(2) = 0.0d0;    
  grid_list%PXg(3) = 0.0d0;

  do b = 1, nb_body     ! Loop over all bodies
     parBegin = body_list(b)%par_begin
     parEnd = body_list(b)%par_End

     if(contact) comID = body_list(b)%comID ! Get comID from body

     do p = parBegin, parEnd    ! Loop over all particles (1)
        pt => particle_list(p)

        pt%icell = InWhichCell(pt%Xp)
        icell = pt%icell
        ! Particle p is out of the computational region
        if (icell < 0) cycle    

        vol_ = pt%VOL
        mp_ = pt%Mass

        ! Calculate the shape functions and their derivatives 
        InflNode(1:8)=CellsNode(icell,:)
        if (GIMP) then
           call FindInflNode(p,icell)
           call NShape_GIMP(p)
        else
           call NShape(InflNode(1),p,0)
        end if

        ! Loop over the grid nodes of the hexhedron 
        !   in which the particle is located
        do n = 1, nb_InflNode 
           ! out of the computational grid
           if (InflNode(n) .gt. nb_gridnode .or. &
               InflNode(n) .le. 0) cycle  

           gd => grid_list(comID, InflNode(n))

           SHPn = SHP(n)
           shm = SHPn*mp_

           gd%Mg = gd%Mg + shm            ! the nodal mass
           gd%PXg = gd%PXg + pt%VXp*shm   ! the nodal momentum
        end do !n

     end do !p
  end do    !b
end subroutine GridMomentumInitial

subroutine ParticleStressUpdate()
!------------------------------------------------------------------
!-  Purpose                                                       -
!-      1. Calculate the strain rate and spin tensor              -
!-      2. Update stresses by appropriate constitution law        -
!------------------------------------------------------------------
  use ParticleData
  use GridData
  use MaterialModel, only: Constitution
  use MaterialData
  implicit none

  integer:: b, p, n, parBegin, parEnd ! loop counter
  integer:: icell, inode, ix, iy, iz, comID = 1
  real(8):: xx(3), vx(3), ax(3), vgx(3)
  real(8):: de(6), vort(3)
  real(8):: mp_, shm, SHPn, DNDXn, DNDYn, DNDZn

  type(Particle), POINTER :: pt
  type(GridNodeProperty), POINTER :: gd

  ! Calculate the increment strain and vorticity
  ! de(i) comply the Voigt rule (d11, d22, d33, 2*d23, 2*d13, 2*d12)
  do b = 1, nb_body
     parBegin = body_list(b)%par_begin
     parEnd = body_list(b)%par_End

     if(contact) comID = body_list(b)%comID  ! Get comID from body

     do p = parBegin, parEnd    ! Loop over all particles (4)
        pt => particle_list(p)

        de   = 0d0    ! Incremental strain
        vort = 0d0    ! Incremental vorticity

        icell = pt%icell    ! use old position
        ! Particle p is out of the computational region
        if (icell < 0) cycle    

        ! Calculate the shape functions and their derivatives 
        InflNode(1:8)=CellsNode(icell,:)
        if (GIMP) then
           call FindInflNode(p,icell)
           call NShape_GIMP(p)
        else
           call NShape(InflNode(1),p,1)
        end if

        ! Loop over all grid nodes of the hexhedron 
        !  in which particle p is located
        do n = 1, nb_InflNode
           if (InflNode(n) .gt. nb_gridnode .or. InflNode(n) .le. 0) &
              cycle  ! out of the computational grid
           gd => grid_list(comID, InflNode(n))
           ! If the nodal mass is not too small
           if (gd%Mg > CutOff) then    
              vgx = gd%PXg / gd%Mg    ! Grid nodal velocity

              DNDXn = DNDX(n);  DNDYn = DNDY(n);  DNDZn = DNDZ(n)
              de(1) = de(1) + DNDXn*vgx(1)            ! D11
              de(2) = de(2) + DNDYn*vgx(2)            ! D22
              de(3) = de(3) + DNDZn*vgx(3)            ! D33
              ! 2*D23
              de(4) = de(4) + (DNDYn*vgx(3) + DNDZn*vgx(2))    
              ! 2*D13
              de(5) = de(5) + (DNDZn*vgx(1) + DNDXn*vgx(3))    
              ! 2*D12
              de(6) = de(6) + (DNDXn*vgx(2) + DNDYn*vgx(1))    

              ! W32
              vort(1) = vort(1) + (DNDYn*vgx(3) - DNDZn*vgx(2))    
              ! W13
              vort(2) = vort(2) + (DNDZn*vgx(1) - DNDXn*vgx(3))    
              ! W21
              vort(3) = vort(3) + (DNDXn*vgx(2) - DNDYn*vgx(1))    
           end if
        end do ! n

        de = de * DT    
        vort = vort * DT / 2d0

        ! Update stress by constitution law
        call Constitution(de, vort, b, p)

        if(.NOT.USF) pt%Xp = pt%XX     ! the next particle position

     end do !p
  end do    !b

end subroutine ParticleStressUpdate

subroutine GridMomentumUpdate()
!------------------------------------------------------------------
!-  Purpose                                                       -
!-      1. calculate the background grid nodal force              -
!------------------------------------------------------------------
  use ParticleData
  use GridData
  use MaterialData
  implicit none

  integer:: b, p, n, c, parBegin, parEnd ! loop counter
  integer:: icell, inode, ix, iy, iz, mat_, comID = 1
  real(8):: sxx, syy, szz, sxy, syz, sxz
  real(8):: fx(3), f_int(3), f_ext(3), mp_, vol_
  real(8):: shm, SHPn, DNDXn, DNDYn, DNDZn

  type(Particle), POINTER :: pt
  type(GridNodeProperty), POINTER :: gd
  type(ContactGridNodeProperty), POINTER :: CP

  ! Calculate the grid nodal forces only

  ! Reset nodal forces
  grid_list%FXg(1) = 0.0d0;    ! Nodal forces
  grid_list%FXg(2) = 0.0d0;
  grid_list%FXg(3) = 0.0d0;

  if(contact) then
     CP_list%ndir(1) = 0.0d0    
     CP_list%ndir(2) = 0.0d0
     CP_list%ndir(3) = 0.0d0
  end if

  do b = 1, nb_body             ! Loop over all bodies
     parBegin = body_list(b)%par_begin
     parEnd = body_list(b)%par_End

     if(contact) comID = body_list(b)%comID  ! Get comID from body
     do p = parBegin, parEnd    ! Loop over all particles
        pt => particle_list(p)

        icell = pt%icell        ! using old position

        ! Particle p is out of the computational region
        if (icell < 0) cycle    

        sxx = pt%SM + pt%SDxx   ! Stresses
        syy = pt%SM + pt%SDyy
        szz = pt%SM + pt%SDzz
        sxy = pt%SDxy
        syz = pt%SDyz
        sxz = pt%SDxz
        ! External forces
        fx = pt%FXp
		if (Gravity) then
			fx = fx + pt%Mass * (body_list(b)%Gravp)    
		end if

        vol_ = pt%VOL
        mp_ = pt%Mass

        ! Calculate the shape functions and their derivatives 
        InflNode(1:8)=CellsNode(icell,:)
        if (GIMP) then
           call FindInflNode(p,icell)
           call NShape_GIMP(p)
        else
           call NShape(InflNode(1),p,2)
        end if

        ! Loop over the grid nodes of the hexhedron 
        !  in which the particle is located
        do n = 1, nb_InflNode 
           if (InflNode(n) .gt. nb_gridnode .or. InflNode(n) .le. 0) &
               cycle  ! out of the computational grid

           gd => grid_list(comID, InflNode(n))

           SHPn = SHP(n)
           DNDXn = DNDX(n)
           DNDYn = DNDY(n)
           DNDZn = DNDZ(n)

           f_int(1) = - (sxx*DNDXn + sxy*DNDYn + sxz*DNDZn)*vol_
           f_int(2) = - (sxy*DNDXn + syy*DNDYn + syz*DNDZn)*vol_
           f_int(3) = - (sxz*DNDXn + syz*DNDYn + szz*DNDZn)*vol_

           f_ext = fx*SHPn

           gd%FXg = gd%FXg + f_int + f_ext  !nodal force

           if(contact) then
              CP => CP_list(comID, InflNode(n)) 
              CP%ndir(1) = CP%ndir(1) + DNDXn*mp_
              CP%ndir(2) = CP%ndir(2) + DNDYn*mp_
              CP%ndir(3) = CP%ndir(3) + DNDZn*mp_
           end if

        end do !n

     end do !p
  end do    !b
end subroutine GridMomentumUpdate

subroutine IntegrateMomentum()
!------------------------------------------------------------------
!-  Purpose                                                       -
!-    1. Integrating  the momentum equations on the               -
!=        computational grid                                      -
!-    2. Apply  boundary conditions                               -
!------------------------------------------------------------------
  use ParticleData
  use GridData
  implicit none

  integer:: c, n ! loop counter
  type(GridNodeProperty), POINTER :: gd
  type(GridNode), POINTER :: node

  do c = 1, nb_component
     do n = 1, nb_gridnode
        gd => grid_list(c, n)
        node => node_list(n)

        ! Integrate momentum equation
        if(istep == 1) then
           gd%PXg = gd%PXg + gd%FXg * DT*0.5
        else
           gd%PXg = gd%PXg + gd%FXg * DT
        end if

        ! Apply boundary conditions on computational grid
        if (node%Fix_x) then  ! Grid node n is fixed in x direction
           gd%PXg(1) = 0.0
           gd%FXg(1) = 0.0
        end if

        if (node%Fix_y) then  ! Grid node n is fixed in y direction
           gd%PXg(2) = 0.0
           gd%FXg(2) = 0.0
        end if

        if (node%Fix_z) then  ! Grid node n is fixed in z direction
           gd%PXg(3) = 0.0
           gd%FXg(3) = 0.0
        end if

     end do !n
  end do    !c
end subroutine IntegrateMomentum

subroutine Lagr_NodContact()
!------------------------------------------------------------------
!-  Purpose                                                       -
!-      1. Establishing the nodal contact criteria,          and  -
!-      2. correct  the normal vectors                       and  -
!-      3. Apply the contact force and adjust nodal velocities    -
!------------------------------------------------------------------
  use  ParticleData
  use  GridData
  use  MaterialData

  implicit none

  integer:: p, n, c ! loop counter
  real(8):: nx, ny, nz, tt,tta,ttb,crit, crita, critb
  real(8):: nomforce, val_fslip, val_fstick, val_ffric, nodtolmg
  real(8):: fstick(3), fslip(3), cforce(3)
  integer:: abody,bbody

  type(GridNodeProperty), POINTER :: gd1
  type(GridNodeProperty), POINTER :: gd2
  type(ContactGridNodeProperty), POINTER :: CP1
  type(ContactGridNodeProperty), POINTER :: CP2

  tot_cont_for = 0.0 ! the total contact force between of bodies

  ! calculate contact force and adjust the nodal force and momentum
  do n = 1, nb_gridnode 
     CP1 =>  CP_list(1,n)
     CP2 =>  CP_list(2,n)
     gd1 => grid_list(1, n)
     gd2 => grid_list(2, n)

     ! recalculate the nodal normal direction    
     ! if normbody 0 then using average method; 
     ! if 1,using abody; if 2,using bbody
     if(normbody == 0)then        
        nx = CP1%ndir(1)  - CP2%ndir(1)
        ny = CP1%ndir(2)  - CP2%ndir(2)
        nz = CP1%ndir(3)  - CP2%ndir(3)
     end if

     if(normbody == 1)then         
        nx = CP1%ndir(1) 
        ny = CP1%ndir(2) 
        nz = CP1%ndir(3) 
     end if

     if(normbody == 2)then        
        nx = - CP2%ndir(1)
        ny = - CP2%ndir(2)
        nz = - CP2%ndir(3)    
     end if

     ! unitize normal vector        
     tt = sqrt(nx*nx + ny*ny + nz*nz)        
     if(tt > epsilon(tt)) then
        nx = nx / tt
        ny = ny / tt
        nz = nz / tt
     end if

     CP1%ndir(1) = nx; ! Nodal direction for contact
     CP1%ndir(2) = ny;    
     CP1%ndir(3) = nz;        
     CP2%ndir = -CP1%ndir        

     crit = 0.0
     ! contact criteria using the unit normal vectors
     if ( gd1%Mg > CutOff .AND. gd2%Mg > CutOff) then
        ! Eq.(3.245)
        crit = (gd1%Pxg(1)*gd2%Mg - gd2%Pxg(1)*gd1%Mg)*nx +&
             (gd1%Pxg(2)*gd2%Mg - gd2%Pxg(2)*gd1%Mg)*ny +&  
             (gd1%Pxg(3)*gd2%Mg - gd2%Pxg(3)*gd1%Mg)*nz  
     end if

     if(crit > epsilon(crit)) then

        tt = (gd1%Mg + gd2%Mg)*Dt

        ! calculate the normal contact force
        nomforce =crit/tt   ! Eq.(3.252)

        ! for friction contact
        if(fricfa > epsilon(fricfa)) then

           ! calculate the contact force   Eq.(3.250)
           cforce = (gd1%Pxg*gd2%Mg - gd2%Pxg*gd1%Mg)/tt 

           ! calculate the tangent contact force
           fstick = cforce - nomforce*CP1%ndir
           val_fstick = sqrt( fstick(1)*fstick(1) +  &
                        fstick(2)*fstick(2) + fstick(3)*fstick(3) )
           val_fslip = fricfa*abs(nomforce)
           if(val_fslip < val_fstick) then
              cforce = nomforce*CP1%ndir + val_fslip*(fstick /val_fstick)
           end if

           ! for contact without friction
        else 
           cforce = nomforce*CP1%ndir
        end if

        ! add contact force to nodal force
        gd1%Fxg = gd1%Fxg - cforce
        gd2%Fxg = gd2%Fxg + cforce

        ! adjust the nodal component by contact force
        gd1%Pxg = gd1%Pxg - cforce * Dt
        gd2%Pxg = gd2%Pxg + cforce * Dt

        tot_cont_for = tot_cont_for + cforce
     end if
  end do  !n

end subroutine Lagr_NodContact

subroutine ParticlePositionUpdate()
!------------------------------------------------------------------
!-  Purpose                                                       -
!-      1. Update particle position and velocity                  -
!------------------------------------------------------------------
  use ParticleData
  use GridData
  use MaterialModel, only: Constitution
  use MaterialData
  implicit none

  integer:: b, p, n, parBegin, parEnd ! loop counter
  integer:: icell, inode, ix, iy, iz, comID = 1
  real(8):: xx(3), vx(3), ax(3), vgx(3)
  real(8):: de(6), vort(3)
  real(8):: mp_, shm, SHPn, DNDXn, DNDYn, DNDZn

  type(Particle), POINTER :: pt
  type(GridNodeProperty), POINTER :: gd

  ! Update particle position and velocity
  do b = 1, nb_body
     parBegin = body_list(b)%par_begin
     parEnd = body_list(b)%par_End
     ! Get comID from body
     if(contact)  comID = body_list(b)%comID

     do p = parBegin, parEnd    ! Loop over all particles (2)
        pt => particle_list(p)

        icell = pt%icell
        ! Particle p is out of the computational region
        if (icell < 0) cycle

        xx = pt%Xp;  ! Particle position at time step k    

        vx = 0d0
        ax = 0d0

        ! Mapping from grid to particle

        ! Calculate the shape functions and their derivatives 
        InflNode(1:8)=CellsNode(icell,:)
        if (GIMP) then
           call FindInflNode(p,icell)
           call NShape_GIMP(p)
        else
           call NShape(InflNode(1),p,2)
        end if

        ! Loop over all grid nodes of the hexhedron 
        !  in which particle p is located
        do n = 1, nb_InflNode
           if(InflNode(n) .gt. nb_gridnode .or. InflNode(n) .le. 0) &
              cycle    ! out of the computational grid

           gd => grid_list(comID, InflNode(n))
           if (gd%Mg > CutOff) then ! The nodal mass is not too small
              SHPn = SHP(n)
              vx = vx + SHPn * (gd%PXg / gd%Mg)
              ax = ax + SHPn * gd%FXg / gd%Mg
           end if

        end do ! n

        ! Time integration
        pt%XX = xx + vx * DT       ! Update particle position
        if(istep == 1) then
           pt%VXp = pt%VXp + ax * DT * 0.5  ! Update particle velocity
        else
           pt%VXp = pt%VXp + ax * DT
        end if
        if(USF)  pt%Xp = pt%XX     ! the next particle position

     end do ! p
  end do    ! b

end subroutine  ParticlePositionUpdate

subroutine GridMomentumMUSL()
!------------------------------------------------------------------
!-  Purpose                                                       -
!-      1. recalculate the grid node momentum by mapping          -
!-         the updated particle information                       -
!-      2. apply boundary condition                               -
!------------------------------------------------------------------
  use ParticleData
  use GridData
  use MaterialModel, only: Constitution
  use MaterialData
  implicit none

  integer:: b, c, p, n, parBegin, parEnd ! loop counter
  integer:: icell, comID = 1
  real(8):: de(6), vort(3)
  real(8):: mp_, shm, SHPn

  type(Particle), POINTER :: pt
  type(GridNodeProperty), POINTER :: gd
  type(GridNode), POINTER :: node

  grid_list%PXg(1) = 0.0d0;
  grid_list%PXg(2) = 0.0d0;
  grid_list%PXg(3) = 0.0d0;

  ! Recalculate the grid node momentum
  do b = 1, nb_body        ! Loop over all bodies
     parBegin = body_list(b)%par_begin
     parEnd = body_list(b)%par_End

     if(contact) comID = body_list(b)%comid   ! Get comID from body

     do p = parBegin, parEnd ! Loop over all particles (3)
        pt => particle_list(p)

        icell = pt%icell
        ! Particle p is out of the computational region
        if (icell < 0) cycle  

        mp_ = pt%mass

        ! Calculate the shape functions
        InflNode(1:8)=CellsNode(icell,:)
        if (GIMP) then
           call FindInflNode(p,icell)
           call NShape_GIMP(p)
        else
           call NShape(InflNode(1),p,0)
        end if

        ! Loop over all grid nodes of the hexhedron 
        !  in which particle p is located
        do n = 1, nb_InflNode
           if(InflNode(n) .gt. nb_gridnode .or. InflNode(n) .le. 0) &
              cycle  ! out of the computational grid
           gd => grid_list(comID, InflNode(n))

           shm = SHP(n)*mp_
           gd%PXg = gd%PXg + pt%VXp*shm

        end do ! n
     end do    ! p
  end do       ! b

  ! Applying essential boundary conditions
  do c = 1, nb_component
     do n = 1, nb_gridnode
        gd => grid_list(c, n)
        node => node_list(n)
        if (node%Fix_x) then
           gd%PXg(1) = 0.0
           gd%FXg(1) = 0.0
        end if

        if (node%Fix_y) then
           gd%PXg(2) = 0.0
           gd%FXg(2) = 0.0
        end if

        if (node%Fix_z) then
           gd%PXg(3) = 0.0
           gd%FXg(3) = 0.0
        end if

     end do !n
  end do    ! c
end subroutine GridMomentumMUSL

subroutine ApplyBoundaryConditions()
!------------------------------------------------------------------
!-  Purpose                                                       -
!-      1. apply boundary condition                               -
!------------------------------------------------------------------
  use GridData
  use ParticleData, only: nb_component
  implicit none

  integer:: n, c ! loop counter
  type(GridNodeProperty), POINTER :: gd
  type(GridNode), POINTER :: node

  do c = 1, nb_component
     do n = 1, nb_gridnode
        gd => grid_list(c, n)
        node => node_list(n)
        ! Apply boundary conditions on computational grid
        if (node%Fix_x) then  ! Grid node n is fixed in x direction
           gd%PXg(1) = 0.0
           gd%FXg(1) = 0.0
        end if

        if (node%Fix_y) then  ! Grid node n is fixed in y direction
           gd%PXg(2) = 0.0
           gd%FXg(2) = 0.0
        end if

        if (node%Fix_z) then  ! Grid node n is fixed in z direction
           gd%PXg(3) = 0.0
           gd%FXg(3) = 0.0
        end if
     end do ! n
  end do    ! c

end subroutine ApplyBoundaryConditions
