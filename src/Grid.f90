
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
! -  Background grid procedures                                    -
! -                                                                -
! -  NOTE: only 8 node cube element is available                   -
! -        only cuboid computational region is available           -
! -  computational region is:                                      -
! -   [SpanX(1),SpanX(2)]*[SpanY(1),SpanY(2)]*[SpanZ(1),SpanZ(2)]  -
! -  Grid,Cell(element) infomation will be computed by program     -
! -  GridNode Numbering is keypoint of searching algorithm         -
! ------------------------------------------------------------------

module GridData

  type GridNode
     real(8):: Xg(3)     ! grid node coordinate
     logical:: Fix_x, Fix_y, Fix_z    ! BC
  end type GridNode

  type GridNodeProperty
     real(8):: Mg        ! mass on grid node
     real(8):: PXg(3)    ! momentum on grid node
     real(8):: FXg(3)    ! internal/external force on gride node
  end type GridNodeProperty

  type ContactGridNodeProperty
     ! the normal direction of contact grid node
     real(8):: ndir(3)  
     ! the tangential unit vetors of contact grid node
     real(8):: sdir(3)  
  end type ContactGridNodeProperty

  type(GridNodeProperty), target, allocatable:: grid_list(:,:)
  type(GridNode), target, allocatable:: node_list(:)
  type(ContactGridNodeProperty), target, allocatable:: CP_list(:,:)
  real(8):: fricfa = 0.0    ! the frictional coefficient
  integer:: normbody = 0    ! the flag of computaional normal
  ! the flag of contact type: 0-no,1-langrange,2-penalty
  integer:: contact_type = 0   
  ! the total contact force between of bodies
  real(8):: tot_cont_for(3)    

  real(8):: SpanX(2)=0.0,  &   ! computational region
            SpanY(2)=0.0,  &
            SpanZ(2)=0.0    

  real(8):: DCell  = 0.0       ! Grid node interval

  real(8):: CutOff = 0.0       ! Grid mass cutoff value

  integer:: nb_gridnode = 0    ! number of gridnodes
  integer:: FixS(6) = 0

  ! Number of cells
  integer:: NumCell=0, NumCellx=0, NumCelly=0, &
            NumCellz=0, NumCellxy=0
  integer:: NGx, NGy, NGz, NGxy
  ! NumCell * 8 - Define computational grid
  integer, allocatable:: CellsNode(:,:)    
  integer:: nb_InflNode = 8
  integer:: InflNode(27)    ! influence node list for each particle
  real(8):: rpg(27,3)       ! distance between particle and node

  real(8):: iJacobi, iJacobi4, iDCell ! 1/Jacobi
  ! shape function and its derivative
  real(8), allocatable:: SHP(:), DNDX(:), DNDY(:), DNDZ(:) 

  ! sign constant used by NShape
  integer,parameter:: SNX(8) = (/-1, 1, 1, -1, -1, 1, 1, -1/), &
                      SNY(8) = (/-1, -1, 1, 1, -1, -1, 1, 1/), &
                      SNZ(8) = (/-1, -1, -1, -1, 1, 1, 1, 1/)

contains

  integer function InWhichCell(xx)
! -----------------------------------------------------------------
! - Purpose                                                       -
! -    Determine which cell the point (xx,yy,zz) is located in    -
! -                                                               -
! - Input                                                         -
! -    xx(3) - Coordinates of a point                             -
! -                                                               -
! - Return values                                                 -
! -    >0 : Cell number in which the point is located (MP)        -
! -    <0 : point out off the cell                                -
! -----------------------------------------------------------------
    implicit none
    real(8), intent(in):: xx(3)

    integer ix, iy, iz

    if (xx(1)<SpanX(1) .or. xx(1)>SpanX(2) .or. xx(2)<SpanY(1) .or. &
        xx(2)>SpanY(2) .or.xx(3)<SpanZ(1) .or. xx(3)>SpanZ(2)) then
       InWhichCell = -1
       return
    end if

    ix = int((xx(1)-SpanX(1))/DCell) + 1
    iy = int((xx(2)-SpanY(1))/DCell) + 1
    iz = int((xx(3)-SpanZ(1))/DCell) + 1
    InWhichCell = (iz - 1)*NumCellxy + (iy - 1)*NumCellx + ix

    if(InWhichCell.gt.NumCell .or. InWhichCell.le.0) then
       InWhichCell = -1
    end if

  end function InWhichCell


  subroutine SetGridData()
! -----------------------------------------------------------------
! - Purpose                                                       -
! -    Create computational grid                                  -
! -----------------------------------------------------------------
    use ParticleData
    use FFI

    implicit none

    integer:: i, j, ix, iy, iz, icell, inode     ! loop counter
    integer:: NP2, Node1, Node5
    character(len=2) Velo
    real(8):: TempV, spx, spy, spz
    real(8):: mat_, mp_

    spx = SpanX(2)-SpanX(1)
    spy = SpanY(2)-SpanY(1)
    spz = SpanZ(2)-SpanZ(1)

    if(DCell.eq.0) then
       stop '*** Error *** DCell must be defined !'
    end if

    if(spx.le.0 .or. spy.le.0 .or. spz.le.0) then
       stop '*** Error *** SPX/SPY/SPZ must be defined !'
    end if

    NumCellx = int(spx/DCell + 0.5)
    NumCelly = int(spy/DCell + 0.5)
    NumCellz = int(spz/DCell + 0.5)

    SpanX(2) = SpanX(1) + NumCellx*DCell
    SpanY(2) = SpanY(1) + NumCelly*DCell
    SpanZ(2) = SpanZ(1) + NumCellz*DCell

    NumCellxy = NumCellx * NumCelly
    NumCell = NumCellxy * NumCellz

    NGx = NumCellx + 1
    NGy = NumCelly + 1
    NGz = NumCellz + 1

    NGxy = NGx*NGy
    nb_gridnode = NGxy*NGz

    print *, 'Number of grid nodes = ', nb_gridnode
    write(iomsg,*)
    write(iomsg,"(a14,i10)") 'Number of grid nodes = ', nb_gridnode
    write(iomsg,"(a14,3i4)") 'cells (x,y,z) ', &
                              numcellx,numcelly,numcellz

    allocate(grid_list(nb_component, nb_gridnode))
    allocate(node_list(nb_gridnode))

    node_list%Fix_x = .false.
    node_list%Fix_y = .false.
    node_list%Fix_z = .false.

    ! create grid node info
    do iz = 1, NGz
       do iy = 1, NGy
          do ix = 1, NGx
             i = (iz - 1)*NGxy + (iy - 1)*NGx + ix
             node_list(i)%Xg(1) = (ix - 1)*DCell + SpanX(1)
             node_list(i)%Xg(2) = (iy - 1)*DCell + SpanY(1)
             node_list(i)%Xg(3) = (iz - 1)*DCell + SpanZ(1)

             if ((ix==1.and.FixS(1)==1).or.(ix==NGx.and.FixS(2)==1).or. &
                 (iy==1.and.FixS(3)==1).or.(iy==NGy.and.FixS(4)==1).or. &
                 (iz==1.and.FixS(5)==1).or.(iz==NGz.and.FixS(6)==1)) then
                node_list(i)%Fix_x = .true.
                node_list(i)%Fix_y = .true.
                node_list(i)%Fix_z = .true.
             end if

             if ((ix==1.and.FixS(1)==2).or.(ix==NGx.and.FixS(2)==2)) then
                node_list(i)%Fix_x = .true.
             end if

             if ((iy==1.and.FixS(3)==2).or.(iy==NGy.and.FixS(4)==2)) then
                node_list(i)%Fix_y = .true.
             end if

             if ((iz==1.and.FixS(5)==2).or.(iz==NGz.and.FixS(6)==2)) then
                node_list(i)%Fix_z = .true.
             end if
          end do
       end do
    end do

    ! create CellsNode
    allocate(CellsNode(NumCell,8))

    ! loop over every cell to create CellsNode Matrix
    do iz = 1, NumCellz
       do iy = 1, NumCelly
          do ix = 1, NumCellx
             i = (iz - 1)*NumCellxy + (iy - 1)*NumCellx + ix
             Node1 = (iz - 1)*NGxy + (iy - 1)*NGx + ix
             Node5 = Node1 + NGxy
             CellsNode(i,1) = Node1
             CellsNode(i,2) = Node1 + 1
             CellsNode(i,3) = Node1 + 1 + NGx
             CellsNode(i,4) = Node1 + NGx
             CellsNode(i,5) = Node5
             CellsNode(i,6) = Node5 + 1
             CellsNode(i,7) = Node5 + 1 + NGx
             CellsNode(i,8) = Node5 +NGx
          end do
       end do
    end do

    ! set parameter used in NShape
    if (GIMP) then
       allocate(SHP(27), DNDX(27), DNDY(27), DNDZ(27))
    else
       allocate(SHP(8), DNDX(8), DNDY(8), DNDZ(8))
    end if

    iJacobi  = 2d0/DCell    !Jacobi = DCell/2d0, iJacobi = 1/Jacobi
    iJacobi4 = iJacobi * 0.125d0    !iJacobi4 = iJacobi/8
    iDCell = 1.0/DCell

  end subroutine SetGridData

  subroutine SetContact_GridNodeData()
! -----------------------------------------------------------------
! - Purpose                                                       -
! -    Create computational grid for simulating contact problem   -
! -----------------------------------------------------------------
    use ParticleData, only:nb_component
    use FFI

    implicit none

    allocate(CP_list(nb_component, nb_gridnode))

  end subroutine SetContact_GridNodeData

  subroutine NShape(node1, p, ider)
!------------------------------------------------------------------
!-  Purpose                                                       -
!-      Evaluate the shape functions and/or their derivatives     -
!-      at particle p associated with nodes of the cell in        -
!-      which the particle p is located                           -
!-  Inputs                                                        -
!-      node1 - number of the first node of the cell in which the -
!-              particle p is located                             -
!-      p     - particle number                                   -
!-      ider  - flag for shape function and derivative calculation-
!-              0 - shape functions only                          -
!-              1 - derivatives only                              -
!-              2 - both shape functions and their derivatives    -
!-  Outputs                                                       -
!-      SHP(8)   - value of shape functions                       -
!-      DNDX(8)  - value of derivative with respect to            -
!-                 X of shape function                            -
!-      DNDY(8)  - value of derivative with respect to            -
!-                 Y of shape function                            -
!-      DNDZ(8)  - value of derivative with respect to            -
!-                 Z of shape function                            -   
!------------------------------------------------------------------
    use ParticleData
    implicit none

    integer, intent(in):: node1, p, ider

    real(8):: x(3) ! nature coordinate ( -1 < x,y,z < 1 )
    real(8):: sx(8), sy(8), sz(8)

    type(Particle), POINTER :: pt
    type(GridNode), POINTER :: node

    node => node_list(node1)
    pt => particle_list(p)

    x = (pt%Xp - node%Xg)*iJacobi - 1d0

    sx = SNX*x(1) + 1d0    ! 1 + xi(i)  * xi
    sy = SNY*x(2) + 1d0    ! 1 + eta(i) * eta
    sz = SNZ*x(3) + 1d0    ! 1 + zeta(i)* zeta

    if (ider .NE. 1) SHP = sx * sy * sz * 0.125d0

    if (ider .NE. 0) then
       DNDX = SNX * sy * sz * iJacobi4
       DNDY = SNY * sx * sz * iJacobi4
       DNDZ = SNZ * sx * sy * iJacobi4
    end if

  end subroutine NShape

  subroutine NShape_GIMP(p)
!------------------------------------------------------------------
!-  Purpose                                                       -
!-      Evaluate the shape functions and/or their derivatives     -
!-      at particle p associated with p's influence nodes (GIMP)  -   
!-  Inputs                                                        -
!-      p     - particle number                                   -
!-  Outputs                                                       -
!-      SHP(27)   - value of shape functions                      -
!-      DNDX(27)  - value of derivative with respect to           -
!-                  X of shape function                           -
!-      DNDY(27)  - value of derivative with respect to           -
!-                  Y of shape function                           -
!-      DNDZ(27)  - value of derivative with respect to           -
!-                  X of shape function                           -
!------------------------------------------------------------------
    use ParticleData
    implicit none

    integer, intent(in):: p
    integer:: i, j

    real(8):: sh(3) = 0.0, dn(3) = 0.0
    real(8):: r(3)

    do i = 1, nb_InflNode
       ! out of the computational grid
       if (InflNode(i) .gt. nb_gridnode .or. &
           InflNode(i) .le. 0) then
           cycle    
       end if
       r = abs(rpg(i,:))
       do j = 1, 3
          if (r(j) < 0.25) then
             sh(j) = (7 - 16 * r(j) * r(j)) * 0.125
             dn(j) = -4 * rpg(i,j)
          else if (r(j) < 0.75) then
             sh(j) = 1 - r(j)
             dn(j) = -sign(real(1.0,8),rpg(i,j))
          else if (r(j) < 1.25) then
             sh(j) = ((5 - 4 * r(j)) * 0.25) ** 2
             dn(j) = 2 * rpg(i,j) - sign(real(1.0,8),rpg(i,j)) * 2.5
          else
             sh(j) = 0.0
             dn(j) = 0.0
          end if
       end do
       SHP(i) = sh(1)*sh(2)*sh(3)
       DNDX(i) = dn(1)*sh(2)*sh(3)*iDCell
       DNDY(i) = sh(1)*dn(2)*sh(3)*iDCell
       DNDZ(i) = sh(1)*sh(2)*dn(3)*iDCell
    end do

  end subroutine NShape_GIMP

  subroutine FindInflNode(p, icell)
!------------------------------------------------------------------
!-  Purpose: Find Influence Nodes for particle p                  -
!-  Inputs                                                        -
!-      p     - particle number                                   -
!-      icell - cell id containing particle p                     -
!-  Outputs                                                       -
!-      InflNode(27)                                              -
!-      rpg(27,3)                                                 -
!------------------------------------------------------------------
    use ParticleData

    integer, intent(in):: p, icell
    integer i
    real(8):: x(3) ! nature coordinate ( -1 < x,y,z < 1 )
    type(Particle), POINTER :: pt

    pt => particle_list(p)

    do i=1, 8
       InflNode(i) = CellsNode(icell,i)
       rpg(i,:) = (pt%Xp - node_list(InflNode(i))%Xg) * iDCell
    end do

    ! lower z and upper z (almost the same code for symmetry)
    if (rpg(1,3) < 0.25 .or. rpg(1,3) > 0.75) then
       if (rpg(1,3) < 0.25) then
          InflNode(9) = InflNode(1) - NGxy
          InflNode(10) = InflNode(2) - NGxy
          InflNode(11) = InflNode(3) - NGxy
          InflNode(12) = InflNode(4) - NGxy
       else
          InflNode(9) = InflNode(5) + NGxy
          InflNode(10) = InflNode(6) + NGxy
          InflNode(11) = InflNode(7) + NGxy
          InflNode(12) = InflNode(8) + NGxy
       end if

       ! lower y
       if (rpg(1,2) < 0.25) then
          InflNode(13) = InflNode(9) - NGx
          InflNode(14) = InflNode(10) - NGx
          InflNode(15) = InflNode(1) - NGx
          InflNode(16) = InflNode(2) - NGx
          InflNode(17) = InflNode(5) - NGx
          InflNode(18) = InflNode(6) - NGx

          ! lower x (1)(10)
          if (rpg(1,1) < 0.25) then
             InflNode(19) = InflNode(1) - 1
             InflNode(20) = InflNode(4) - 1
             InflNode(21) = InflNode(5) - 1
             InflNode(22) = InflNode(8) - 1
             InflNode(23) = InflNode(9) - 1
             InflNode(24) = InflNode(12) - 1
             InflNode(25) = InflNode(13) - 1
             InflNode(26) = InflNode(15) - 1
             InflNode(27) = InflNode(17) - 1
             nb_InflNode = 27

          ! upper x (2)(11)
          else if (rpg(1,1) > 0.75) then
             InflNode(19) = InflNode(2) + 1
             InflNode(20) = InflNode(3) + 1
             InflNode(21) = InflNode(6) + 1
             InflNode(22) = InflNode(7) + 1
             InflNode(23) = InflNode(10) + 1
             InflNode(24) = InflNode(11) + 1
             InflNode(25) = InflNode(14) + 1
             InflNode(26) = InflNode(16) + 1
             InflNode(27) = InflNode(18) + 1
             nb_InflNode = 27

          ! middle x (3)(12)
          else 
             nb_InflNode = 18
          end if

          ! upper y
       else if (rpg(1,2) > 0.75) then
          InflNode(13) = InflNode(11) + NGx
          InflNode(14) = InflNode(12) + NGx
          InflNode(15) = InflNode(3) + NGx
          InflNode(16) = InflNode(4) + NGx
          InflNode(17) = InflNode(7) + NGx
          InflNode(18) = InflNode(8) + NGx

          ! lower x (4)(13)
          if (rpg(1,1) < 0.25) then
             InflNode(19) = InflNode(1) - 1
             InflNode(20) = InflNode(4) - 1
             InflNode(21) = InflNode(5) - 1
             InflNode(22) = InflNode(8) - 1
             InflNode(23) = InflNode(9) - 1
             InflNode(24) = InflNode(12) - 1
             InflNode(25) = InflNode(14) - 1
             InflNode(26) = InflNode(16) - 1
             InflNode(27) = InflNode(18) - 1
             nb_InflNode = 27

          ! upper x (5)(14)
          else if (rpg(1,1) > 0.75) then
             InflNode(19) = InflNode(2) + 1
             InflNode(20) = InflNode(3) + 1
             InflNode(21) = InflNode(6) + 1
             InflNode(22) = InflNode(7) + 1
             InflNode(23) = InflNode(10) + 1
             InflNode(24) = InflNode(11) + 1
             InflNode(25) = InflNode(13) + 1
             InflNode(26) = InflNode(15) + 1
             InflNode(27) = InflNode(17) + 1
             nb_InflNode = 27

          ! middle x (6)(15)
          else 
             nb_InflNode = 18
          end if

          ! middle y
       else

          ! lower x (7)(16)
          if (rpg(1,1) < 0.25) then
             InflNode(13) = InflNode(9) - 1
             InflNode(14) = InflNode(12) - 1
             InflNode(15) = InflNode(1) - 1
             InflNode(16) = InflNode(4) - 1
             InflNode(17) = InflNode(5) - 1
             InflNode(18) = InflNode(8) - 1
             nb_InflNode = 18

          ! upper x (8)(17)
          else if (rpg(1,1) > 0.75) then
             InflNode(13) = InflNode(10) + 1
             InflNode(14) = InflNode(11) + 1
             InflNode(15) = InflNode(2) + 1
             InflNode(16) = InflNode(3) + 1
             InflNode(17) = InflNode(6) + 1
             InflNode(18) = InflNode(7) + 1
             nb_InflNode = 18

          ! middle x (9)(18)
          else 
             nb_InflNode = 12
          end if

       end if
    else
       ! lower y
       if (rpg(1,2) < 0.25) then
          InflNode(9) = InflNode(1) - NGx
          InflNode(10) = InflNode(2) - NGx
          InflNode(11) = InflNode(5) - NGx
          InflNode(12) = InflNode(6) - NGx

          ! lower x (19)
          if (rpg(1,1) < 0.25) then
             InflNode(13) = InflNode(1) - 1
             InflNode(14) = InflNode(4) - 1
             InflNode(15) = InflNode(5) - 1
             InflNode(16) = InflNode(8) - 1
             InflNode(17) = InflNode(9) - 1
             InflNode(18) = InflNode(11) - 1
             nb_InflNode = 18

          ! upper x (20)
          else if (rpg(1,1) > 0.75) then
             InflNode(13) = InflNode(2) + 1
             InflNode(14) = InflNode(3) + 1
             InflNode(15) = InflNode(6) + 1
             InflNode(16) = InflNode(7) + 1
             InflNode(17) = InflNode(10) + 1
             InflNode(18) = InflNode(12) + 1
             nb_InflNode = 18

          ! middle x (21)
          else 
             nb_InflNode = 12
          end if

          ! upper y
       else if (rpg(1,2) > 0.75) then
          InflNode(9) = InflNode(3) + NGx
          InflNode(10) = InflNode(4) + NGx
          InflNode(11) = InflNode(7) + NGx
          InflNode(12) = InflNode(8) + NGx

          ! lower x (22)
          if (rpg(1,1) < 0.25) then
             InflNode(13) = InflNode(1) - 1
             InflNode(14) = InflNode(4) - 1
             InflNode(15) = InflNode(5) - 1
             InflNode(16) = InflNode(8) - 1
             InflNode(17) = InflNode(10) - 1
             InflNode(18) = InflNode(12) - 1
             nb_InflNode = 18

          ! upper x (23)
          else if (rpg(1,1) > 0.75) then
             InflNode(13) = InflNode(2) + 1
             InflNode(14) = InflNode(3) + 1
             InflNode(15) = InflNode(6) + 1
             InflNode(16) = InflNode(7) + 1
             InflNode(17) = InflNode(9) + 1
             InflNode(18) = InflNode(11) + 1
             nb_InflNode = 18

          ! middle x (24)
          else 
             nb_InflNode = 12
          end if

          ! middle y
       else

          ! lower x (25)
          if (rpg(1,1) < 0.25) then
             InflNode(9) = InflNode(1) - 1
             InflNode(10) = InflNode(4) - 1
             InflNode(11) = InflNode(5) - 1
             InflNode(12) = InflNode(8) - 1
             nb_InflNode = 12

          ! upper x (26)
          else if (rpg(1,1) > 0.75) then
             InflNode(9) = InflNode(2) + 1
             InflNode(10) = InflNode(3) + 1
             InflNode(11) = InflNode(6) + 1
             InflNode(12) = InflNode(7) + 1
             nb_InflNode = 12

          ! middle x (27)
          else 
             nb_InflNode = 8
          end if
       end if
    end if

    do i = 9, nb_InflNode
       if (InflNode(i) .le. nb_gridnode .and. InflNode(i) .gt. 0) then
          rpg(i,:)= pt%Xp - node_list(InflNode(i))%Xg
       else    ! out of the computational grid
          rpg(i,:) = 100.0
       end if
       rpg(i,:) = rpg(i,:) * iDCell
    end do

  end subroutine FindInflNode

end module GridData

