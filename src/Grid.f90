
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
     integer:: bd_type(3)=0       ! the type of the boundary 1-left boundary 2-right boundary 0-not boundary
  end type GridNode

  type GridNodeProperty
     real(8):: Mg        ! mass on grid node
     real(8):: PXg(3)    ! momentum on grid node
     real(8):: FXg(3)    ! internal/external force on grid node
     real(8):: Gpre=0.0  ! pressure on grid node
     real(8):: GV_all=0.0!Volume on grid node mapped by particles
     logical:: mapped = .false. !Is there a mapping of material points to this node during post-processing?
  end type GridNodeProperty
  
 type CellData
     real(8)::Cxg(3)    !the node coordinate of the cell center
 end type CellData
 
 type CellDataproperty
     real(8)::Cmg       ! mass on the center node of the cell
     real(8)::Cpxg(3)   ! momentun on the center node of the cell
     real(8)::Cfxg(3)   ! total grid force on the center node of the cell
     real(8)::Cfint(6)  ! internal/external force on the center gride node
     real(8)::Cfext(3)
     real(8)::Cvx(3)    ! the velocity of the center node in the t+1/2
     real(8)::Cax(3)    ! the accleration of the center node in the t+1/2
     real(8)::Co(6)     ! the strain increment and vorticity increment at the center node
     real(8)::Cw(3)
     real(8)::CdeFp(3,3)    
     real(8)::CPKint(9) !1st PK stress for TLMPM
     real(8)::CGpre     ! pressure on auxiliary node
     real(8):: CGV_all=0.0 !Volume on auxiliary node mapped by particles
     logical:: Cmapped = .false. !Is there a mapping of material points to this node during post-processing?
  end type CellDataproperty

  type ContactGridNodeProperty
     ! the normal direction of contact grid node
     real(8):: ndir(3)  
     ! the tangential unit vetors of contact grid node
     real(8):: sdir(3)  
  end type ContactGridNodeProperty

  type(GridNodeProperty), target, allocatable:: grid_list(:,:)
  type(GridNode), target, allocatable:: node_list(:)
  type(ContactGridNodeProperty), target, allocatable:: CP_list(:,:)
  type(CellData),target,allocatable::cell_list(:)
  type(CellDataproperty),target,allocatable::cellp_list(:)
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
  integer:: nb_centernode = 0  ! number of the center node of the cell
  integer:: FixS(6) = 0

  ! Number of cells
  integer:: NumCell=0, NumCellx=0, NumCelly=0, &
            NumCellz=0, NumCellxy=0,CenterNumCell=0,&
            CenterNumCellx=0,CenterNumCelly=0,CenterNumCellz,&
            CenterNumCellxy=0
  integer:: NGx, NGy, NGz, NGxy
  ! NumCell * 8 - Define computational grid
  integer, allocatable:: CellsNode(:,:),CenterCellNode(:,:)    
  integer:: nb_InflNode = 8
  integer:: centernb_InflNode=8
  integer:: InflNode(27)    ! influence node list for each particle
  integer:: CenterInflNode(8) !influence node list for the Auxiliary grid
  real(8):: rpg(27,3)       ! distance between particle and node

  real(8):: iJacobi, iJacobi4, iDCell ! 1/Jacobi
  ! shape function and its derivative
  real(8), allocatable:: SHP(:), DNDX(:), DNDY(:), DNDZ(:) ,CSHP(:),CDNDX(:),CDNDY(:),CDNDZ(:)

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
  
  integer function CenterInWhichCell(xxx)
! - find the position which the material point in auxiliary grid
! -----------------------------------------------------------------
! - Purpose                                                       -
! -    Determine which cell the point (xxx,yyy,zzz) is located in    -
! -                                                               -
! - Input                                                         -
! -    xxx(3) - Coordinates of a point                             -
! -                                                               -
! - Return values                                                 -
! -    >0 : Cell number in which the point is located (MP)        -
! -    <0 : point out off the cell                                -
! -----------------------------------------------------------------
    implicit none
    real(8), intent(in):: xxx(3)

    integer ix, iy, iz

    !if (xxx(1)<SpanX(1)+DCell/2.0 .or. xxx(1)>SpanX(2)-DCell/2.0 .or. xxx(2)<SpanY(1)+DCell/2.0 .or. &
    !    xxx(2)>SpanY(2)-DCell/2.0 .or.xxx(3)<SpanZ(1)+DCell/2.0 .or. xxx(3)>SpanZ(2)-DCell/2.0) then
    !   CenterInWhichCell = -1
    !   return
    !end if

    ix = int((xxx(1)-SpanX(1)-DCell/2.0)/DCell) + 1
    iy = int((xxx(2)-SpanY(1)-DCell/2.0)/DCell) + 1
    iz = int((xxx(3)-SpanZ(1)-DCell/2.0)/DCell) + 1
    CenterInWhichCell = (iz - 1)*CenterNumCellxy + (iy - 1)*CenterNumCellx + ix

    if(CenterInWhichCell.gt.CenterNumCell .or. CenterInWhichCell.le.0) then
       CenterInWhichCell = -1
    end if

  end function CenterInWhichCell

   
function AuxiliaryGridSHP(XN,XXN,p_type)result(XINB)
! -----------------------------------------------------------------
! - Purpose                                                       -
! -    Calculate grid SHP in the boundary auxiliary grid          -
! -                                                               -
! - Input                                                         -
! -    XN(3) - particle natural coordinates in  background grid   -
! -    XXN(3) - nature Coordinates of a point in auxiliary grid   -
! -    p_type(3)   the type of the boundary of particle           -                                          -
! - Return values                                                 -
! -    XINB(3) : grid SHP in the boundary auxiliary grid          -  
! -----------------------------------------------------------------
    implicit none
    real(8), intent(in):: XXN(3),XN(3)
    real(8):: XINB(3)
    integer:: i,p_type(3)
    do i = 1,3
        if(p_type(i))then
            if(XN(i)<0) then
                XINB(i)=1+2*XN(i)
            else
                XINB(i)=2*XN(i)-1
            end if
        else
            XINB(i)=XXN(i)
        end if
    end do
  end function  AuxiliaryGridSHP 

  subroutine SetGridData()
! -----------------------------------------------------------------
! - Purpose                                                       -
! -    Create computational grid                                  -
! -----------------------------------------------------------------
    use ParticleData
    use FFI

    implicit none

    integer:: i, j, ix, iy, iz, icell, inode,centericell     ! loop counter
    integer:: NP2, Node1, Node5,CenterNode1,CenterNode5
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
    
    CenterNumCellx=NumCellx-1
    CenterNumCelly=NumCelly-1
    CenterNumCellz=NumCellz-1
    if (CenterNumCellx == 0) CenterNumCellx=1
    if (CenterNumCelly == 0) CenterNumCelly=1
    if (CenterNumCellz == 0) CenterNumCellz=1

    SpanX(2) = SpanX(1) + NumCellx*DCell
    SpanY(2) = SpanY(1) + NumCelly*DCell
    SpanZ(2) = SpanZ(1) + NumCellz*DCell

    NumCellxy = NumCellx * NumCelly
    CenterNumCellxy=CenterNumCellx*CenterNumCelly
    NumCell = NumCellxy * NumCellz
    CenterNumCell=CenterNumCellxy*CenterNumCellz
    nb_centernode=NumCell

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
    allocate(cell_list(nb_centernode))
    if(SGMP)then
    allocate(cellp_list(nb_centernode))
    end if

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
             if(ix==1.and.FixS(1))    node_list(i)%bd_type(1)=1
             if(ix==NGx.and.FixS(2))  node_list(i)%bd_type(1)=2
             if(iy==1.and.FixS(3))    node_list(i)%bd_type(2)=1
             if(iy==NGy.and.FixS(4))  node_list(i)%bd_type(2)=2
             if(iz==1.and.FixS(5))    node_list(i)%bd_type(3)=1
             if(iz==NGz.and.FixS(6))  node_list(i)%bd_type(3)=2
             if(sgmp) then
                  if(ix==2.and.FixS(1)/=0)      node_list(i)%bd_type(1)=1
                  if(ix==NGx-1.and.FixS(2)/=0)  node_list(i)%bd_type(1)=2
                  if(iy==2.and.FixS(3)/=0)      node_list(i)%bd_type(2)=1
                  if(iy==NGy-1.and.FixS(4)/=0)  node_list(i)%bd_type(2)=2
                  if(iz==2.and.FixS(5)/=0)      node_list(i)%bd_type(3)=1
                  if(iz==NGz-1.and.FixS(6)/=0)  node_list(i)%bd_type(3)=2
             end if
        if(Bspline.or.gimp.or.sgmp)then
                 if (((ix==1.or.ix==2).and.FixS(1)==1).or.((ix==NGx.or.ix==NGx-1).and.FixS(2)==1).or. &
                 ((iy==1.or.iy==2).and.FixS(3)==1).or.((iy==NGy.or.iy==NGy-1).and.FixS(4)==1).or. &
                 ((iz==1.or.iz==2).and.FixS(5)==1).or.((iz==NGz.or.iz==NGz-1).and.FixS(6)==1)) then
                node_list(i)%Fix_x = .true.
                node_list(i)%Fix_y = .true.
                node_list(i)%Fix_z = .true.
                end if

             if (((ix==1.or.ix==2).and.FixS(1)==2).or.((ix==NGx.or.ix==NGx-1).and.FixS(2)==2)) then
                node_list(i)%Fix_x = .true.
             end if
             
             if (((iy==1.or.iy==2).and.FixS(3)==2).or.((iy==NGy.or.iy==NGy-1).and.FixS(4)==2)) then
                node_list(i)%Fix_y = .true.
             end if
             
             if (((iz==1.or.iz==2).and.FixS(5)==2).or.((iz==NGz.or.iz==NGz-1).and.FixS(6)==2)) then
                node_list(i)%Fix_z = .true.
             end if
        end if
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
    if(SGMP)allocate(CenterCellNode(CenterNumCell,8))

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
             cell_list(i)%Cxg=node_list(CellsNode(i,1))%Xg+DCell/2.0
          end do
       end do
    end do
    
if(SGMP)then
    do iz = 1, CenterNumCellz
       do iy = 1, CenterNumCelly
          do ix = 1, CenterNumCellx
             i = (iz - 1)*CenterNumCellxy + (iy - 1)*CenterNumCellx + ix
             CenterNode1 = (iz - 1)*NumCellxy + (iy - 1)*NumCellx + ix
             CenterNode5 = CenterNode1 + NumCellxy
             CenterCellNode(i,1) = CenterNode1
             CenterCellNode(i,2) = CenterNode1 + 1
             CenterCellNode(i,3) = CenterNode1 + 1 + NumCellx
             CenterCellNode(i,4) = CenterNode1 + NumCellx
             CenterCellNode(i,5) = CenterNode5
             CenterCellNode(i,6) = CenterNode5 + 1
             CenterCellNode(i,7) = CenterNode5 + 1 + NumCellx
             CenterCellNode(i,8) = CenterNode5 +NumCellx
          end do
       end do
    end do
end if

    ! set parameter used in NShape
    if (GIMP.or.Bspline) then
       allocate(SHP(27), DNDX(27), DNDY(27), DNDZ(27))
    else if(SGMP)then
       allocate(SHP(8), DNDX(8), DNDY(8), DNDZ(8),CSHP(8),CDNDX(8),CDNDY(8),CDNDZ(8))
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
  
  subroutine SGNShape(xxn)
!------------------------------------------------------------------
!- Shape function for SGMP
!-  Purpose                                                       -
!-      Evaluate the shape functions and                          -
!-      at particle p associated with nodes of the cell in        -
!-      which the particle p is located                           -
!-  Inputs                                                        -
!-      xxn - particle natural coordinates in the auxiliary grid  -
!-  Outputs                                                       -
!-      SHP(8)   - value of shape functions                       -
!------------------------------------------------------------------
    use ParticleData
    implicit none
    real(8), intent(in):: xxn(3)! nature coordinate ( -1 < x,y,z < 1 )

    real(8):: sx(8), sy(8), sz(8)

    sx = SNX*xxn(1) + 1d0    ! 1 + xi(i)  * xi
    sy = SNY*xxn(2) + 1d0    ! 1 + eta(i) * eta
    sz = SNZ*xxn(3) + 1d0    ! 1 + zeta(i)* zeta

   SHP = sx * sy * sz * 0.125d0

  end subroutine SGNShape
  

  subroutine NShape(XN,ider) 
!------------------------------------------------------------------
!-  Purpose                                                       -
!-      Evaluate the shape functions and/or their derivatives     -
!-      at particle p associated with nodes of the cell in        -
!-      which the particle p is located                           -
!-  Inputs                                                        -
!-  XN(3) - particle natural coordinates in the background grid   -
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
    implicit none

    integer, intent(in)::  ider
    real(8),intent(in)::XN(3)! nature coordinate ( -1 < x,y,z < 1 )

    real(8):: sx(8), sy(8), sz(8)
   
    sx = SNX*XN(1) + 1d0    ! 1 + xi(i)  * xi
    sy = SNY*XN(2) + 1d0    ! 1 + eta(i) * eta
    sz = SNZ*XN(3) + 1d0    ! 1 + zeta(i)* zeta

    if (ider .NE. 1) SHP = sx * sy * sz * 0.125d0

    if (ider .NE. 0) then
       DNDX = SNX * sy * sz * iJacobi4
       DNDY = SNY * sx * sz * iJacobi4
       DNDZ = SNZ * sx * sy * iJacobi4
    end if

  end subroutine NShape
  
  subroutine CNShape(ider)
!------------------------------------------------------------------
!- Shape function for SGMP
!-  Purpose                                                       -
!-      Evaluate the shape functions and/or their derivatives     -
!-      at the centerbode associated with nodes of the cell in        -
!-      which the particle p is located                           -
!-  Inputs                                                        -                        -                               -
!-      ider  - flag for shape function and derivative calculation-
!-              0 - shape functions only                          -
!-              1 - derivatives only                              -
!-              2 - both shape functions and their derivatives    -
!-  Outputs                                                       -
!-      CSHP(8)   - value of shape functions                       -
!-      CDNDX(8)  - value of derivative with respect to            -
!-                  X of shape function                            -
!-      CDNDY(8)  - value of derivative with respect to            -
!-                  Y of shape function                            -
!-      CDNDZ(8)  - value of derivative with respect to            -
!-                  Z of shape function                            -   
!------------------------------------------------------------------
    use ParticleData
    implicit none

    integer, intent(in)::  ider

    real(8):: x(3)=(/0,0,0/) ! nature coordinate ( -1 < x,y,z < 1 )
    real(8):: sx(8), sy(8), sz(8)

    sx = SNX*x(1) + 1d0    ! 1 + xi(i)  * xi
    sy = SNY*x(2) + 1d0    ! 1 + eta(i) * eta
    sz = SNZ*x(3) + 1d0    ! 1 + zeta(i)* zeta

    !if (ider .NE. 1) CSHP = sx * sy * sz * 0.125d0

    if (ider .NE. 0) then
       CDNDX = SNX * sy * sz * iJacobi4
       CDNDY = SNY * sx * sz * iJacobi4
       CDNDZ = SNZ * sx * sy * iJacobi4
    end if

  end subroutine CNShape
  
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


    subroutine BFindInflNode(icell)
!------------------------------------------------------------------
!-  Purpose: Find Influence Nodes of B-spline for particle p                  -
!-  Inputs                                                        -                                  -
!-      icell - cell id containing particle p                     -
!-  Outputs                                                       -
!-      InflNode(27)                                              -                                               -
!------------------------------------------------------------------
    use ParticleData

    integer, intent(in):: icell
    integer i

    do i=1, 8
       InflNode(i) = CellsNode(icell,i)
    end do
        ! lower z
          InflNode(9) = InflNode(1) - NGxy
          InflNode(10) = InflNode(2) - NGxy
          InflNode(11) = InflNode(3) - NGxy
          InflNode(12) = InflNode(4) - NGxy

       ! lower y
          
          InflNode(13) = InflNode(9) - NGx
          InflNode(14) = InflNode(10) - NGx
          InflNode(15) = InflNode(1) - NGx
          InflNode(16) = InflNode(2) - NGx
          InflNode(17) = InflNode(5) - NGx
          InflNode(18) = InflNode(6) - NGx

        ! lower x 
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

    end subroutine BFindInflNode

    
      subroutine NShape_Bspline(p)
!------------------------------------------------------------------
!-  Purpose                                                       -
!-      Evaluate the shape functions and/or their derivatives     -
!-      at particle p associated with p's influence nodes (BSplinemethod)  -   
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
    type(Particle), POINTER :: pt
    
    real(8):: sh(3,3) = 0.0, dn(3,3) = 0.0
    real(8):: r(3)

    pt => particle_list(p)


    r = (pt%Xp - node_list(InflNode(1))%Xg) * iDCell*0.2+0.4
    sh(1,:)=0.5*(3-5*r)**2
    sh(2,:)=-25*r**2+25*r-5.5
    sh(3,:)=0.5*(5*r-2)**2
    !cauclate the SHP
    SHP(1)=sh(2,1)*sh(2,2)*sh(2,3)
    SHP(2)=sh(3,1)*sh(2,2)*sh(2,3)
    SHP(3)=sh(3,1)*sh(3,2)*sh(2,3)
    SHP(4)=sh(2,1)*sh(3,2)*sh(2,3)
    SHP(5)=sh(2,1)*sh(2,2)*sh(3,3)
    SHP(6)=sh(3,1)*sh(2,2)*sh(3,3)
    SHP(7)=sh(3,1)*sh(3,2)*sh(3,3)
    SHP(8)=sh(2,1)*sh(3,2)*sh(3,3)
    SHP(9)=sh(2,1)*sh(2,2)*sh(1,3)
    SHP(10)=sh(3,1)*sh(2,2)*sh(1,3)
    SHP(11)=sh(3,1)*sh(3,2)*sh(1,3)
    SHP(12)=sh(2,1)*sh(3,2)*sh(1,3)
    SHP(13)=sh(2,1)*sh(1,2)*sh(1,3)
    SHP(14)=sh(3,1)*sh(1,2)*sh(1,3)
    SHP(15)=sh(2,1)*sh(1,2)*sh(2,3)
    SHP(16)=sh(3,1)*sh(1,2)*sh(2,3)
    SHP(17)=sh(2,1)*sh(1,2)*sh(3,3)
    SHP(18)=sh(3,1)*sh(1,2)*sh(3,3)
    SHP(19)=sh(1,1)*sh(2,2)*sh(2,3)
    SHP(20)=sh(1,1)*sh(3,2)*sh(2,3)
    SHP(21)=sh(1,1)*sh(2,2)*sh(3,3)
    SHP(22)=sh(1,1)*sh(3,2)*sh(3,3)
    SHP(23)=sh(1,1)*sh(2,2)*sh(1,3)
    SHP(24)=sh(1,1)*sh(3,2)*sh(1,3)
    SHP(25)=sh(1,1)*sh(1,2)*sh(1,3)
    SHP(26)=sh(1,1)*sh(1,2)*sh(2,3)
    SHP(27)=sh(1,1)*sh(1,2)*sh(3,3)

    dn(1,:)=(5*r-3)*5*0.2*iDCell
    dn(2,:)=(-50*r+25)*0.2*iDCell
    dn(3,:)=(5*r-2)*5*0.2*iDCell
    !cauclate the DNDX
    DNDX(1)=dn(2,1)*sh(2,2)*sh(2,3)
    DNDX(2)=dn(3,1)*sh(2,2)*sh(2,3)
    DNDX(3)=dn(3,1)*sh(3,2)*sh(2,3)
    DNDX(4)=dn(2,1)*sh(3,2)*sh(2,3)
    DNDX(5)=dn(2,1)*sh(2,2)*sh(3,3)
    DNDX(6)=dn(3,1)*sh(2,2)*sh(3,3)
    DNDX(7)=dn(3,1)*sh(3,2)*sh(3,3)
    DNDX(8)=dn(2,1)*sh(3,2)*sh(3,3)
    DNDX(9)=dn(2,1)*sh(2,2)*sh(1,3)
    DNDX(10)=dn(3,1)*sh(2,2)*sh(1,3)
    DNDX(11)=dn(3,1)*sh(3,2)*sh(1,3)
    DNDX(12)=dn(2,1)*sh(3,2)*sh(1,3)
    DNDX(13)=dn(2,1)*sh(1,2)*sh(1,3)
    DNDX(14)=dn(3,1)*sh(1,2)*sh(1,3)
    DNDX(15)=dn(2,1)*sh(1,2)*sh(2,3)
    DNDX(16)=dn(3,1)*sh(1,2)*sh(2,3)
    DNDX(17)=dn(2,1)*sh(1,2)*sh(3,3)
    DNDX(18)=dn(3,1)*sh(1,2)*sh(3,3)
    DNDX(19)=dn(1,1)*sh(2,2)*sh(2,3)
    DNDX(20)=dn(1,1)*sh(3,2)*sh(2,3)
    DNDX(21)=dn(1,1)*sh(2,2)*sh(3,3)
    DNDX(22)=dn(1,1)*sh(3,2)*sh(3,3)
    DNDX(23)=dn(1,1)*sh(2,2)*sh(1,3)
    DNDX(24)=dn(1,1)*sh(3,2)*sh(1,3)
    DNDX(25)=dn(1,1)*sh(1,2)*sh(1,3)
    DNDX(26)=dn(1,1)*sh(1,2)*sh(2,3)
    DNDX(27)=dn(1,1)*sh(1,2)*sh(3,3)

    !cauclate the DNDY
    DNDY(1)=sh(2,1)*dn(2,2)*sh(2,3)
    DNDY(2)=sh(3,1)*dn(2,2)*sh(2,3)
    DNDY(3)=sh(3,1)*dn(3,2)*sh(2,3)
    DNDY(4)=sh(2,1)*dn(3,2)*sh(2,3)
    DNDY(5)=sh(2,1)*dn(2,2)*sh(3,3)
    DNDY(6)=sh(3,1)*dn(2,2)*sh(3,3)
    DNDY(7)=sh(3,1)*dn(3,2)*sh(3,3)
    DNDY(8)=sh(2,1)*dn(3,2)*sh(3,3)
    DNDY(9)=sh(2,1)*dn(2,2)*sh(1,3)
    DNDY(10)=sh(3,1)*dn(2,2)*sh(1,3)
    DNDY(11)=sh(3,1)*dn(3,2)*sh(1,3)
    DNDY(12)=sh(2,1)*dn(3,2)*sh(1,3)
    DNDY(13)=sh(2,1)*dn(1,2)*sh(1,3)
    DNDY(14)=sh(3,1)*dn(1,2)*sh(1,3)
    DNDY(15)=sh(2,1)*dn(1,2)*sh(2,3)
    DNDY(16)=sh(3,1)*dn(1,2)*sh(2,3)
    DNDY(17)=sh(2,1)*dn(1,2)*sh(3,3)
    DNDY(18)=sh(3,1)*dn(1,2)*sh(3,3)
    DNDY(19)=sh(1,1)*dn(2,2)*sh(2,3)
    DNDY(20)=sh(1,1)*dn(3,2)*sh(2,3)
    DNDY(21)=sh(1,1)*dn(2,2)*sh(3,3)
    DNDY(22)=sh(1,1)*dn(3,2)*sh(3,3)
    DNDY(23)=sh(1,1)*dn(2,2)*sh(1,3)
    DNDY(24)=sh(1,1)*dn(3,2)*sh(1,3)
    DNDY(25)=sh(1,1)*dn(1,2)*sh(1,3)
    DNDY(26)=sh(1,1)*dn(1,2)*sh(2,3)
    DNDY(27)=sh(1,1)*dn(1,2)*sh(3,3) 
    
    !caculate the DNDZ
    DNDZ(1)=sh(2,1)*sh(2,2)*dn(2,3)
    DNDZ(2)=sh(3,1)*sh(2,2)*dn(2,3)
    DNDZ(3)=sh(3,1)*sh(3,2)*dn(2,3)
    DNDZ(4)=sh(2,1)*sh(3,2)*dn(2,3)
    DNDZ(5)=sh(2,1)*sh(2,2)*dn(3,3)
    DNDZ(6)=sh(3,1)*sh(2,2)*dn(3,3)
    DNDZ(7)=sh(3,1)*sh(3,2)*dn(3,3)
    DNDZ(8)=sh(2,1)*sh(3,2)*dn(3,3)
    DNDZ(9)=sh(2,1)*sh(2,2)*dn(1,3)
    DNDZ(10)=sh(3,1)*sh(2,2)*dn(1,3)
    DNDZ(11)=sh(3,1)*sh(3,2)*dn(1,3)
    DNDZ(12)=sh(2,1)*sh(3,2)*dn(1,3)
    DNDZ(13)=sh(2,1)*sh(1,2)*dn(1,3)
    DNDZ(14)=sh(3,1)*sh(1,2)*dn(1,3)
    DNDZ(15)=sh(2,1)*sh(1,2)*dn(2,3)
    DNDZ(16)=sh(3,1)*sh(1,2)*dn(2,3)
    DNDZ(17)=sh(2,1)*sh(1,2)*dn(3,3)
    DNDZ(18)=sh(3,1)*sh(1,2)*dn(3,3)
    DNDZ(19)=sh(1,1)*sh(2,2)*dn(2,3)
    DNDZ(20)=sh(1,1)*sh(3,2)*dn(2,3)
    DNDZ(21)=sh(1,1)*sh(2,2)*dn(3,3)
    DNDZ(22)=sh(1,1)*sh(3,2)*dn(3,3)
    DNDZ(23)=sh(1,1)*sh(2,2)*dn(1,3)
    DNDZ(24)=sh(1,1)*sh(3,2)*dn(1,3)
    DNDZ(25)=sh(1,1)*sh(1,2)*dn(1,3)
    DNDZ(26)=sh(1,1)*sh(1,2)*dn(2,3)
    DNDZ(27)=sh(1,1)*sh(1,2)*dn(3,3)
      end subroutine NShape_Bspline
      
end module GridData
