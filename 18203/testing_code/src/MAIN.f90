PROGRAM MAIN

!USE TYPES
!USE PARAM
!USE SIM_PARAM
!USE IFPORT
USE mpi

IMPLICIT NONE

! NPROC IS HARD DEFINED HERE...
INTEGER, PARAMETER :: nproc = 12

INTEGER, PARAMETER :: RPREC = KIND(1.D0) 
INTEGER :: ip, np, coords(1), kkk(1)
INTEGER :: status(MPI_STATUS_SIZE)
CHARACTER (8) :: chcoord  !--holds character representation of coord
INTEGER :: ierr
INTEGER :: comm
INTEGER :: up, down
INTEGER :: MPI_RPREC, MPI_CPREC
INTEGER :: rank = 0
INTEGER :: coord = 0 
INTEGER :: rank_of_coord(0:nproc-1), coord_of_rank(0:nproc-1)

!######################################################################
!GEOMETRY AND SPATIAL DISCRETIZATION
REAL(RPREC) :: L_X,L_Y,L_Z,T1,T2,T1A,T2A
INTEGER, PARAMETER :: NX=512,NY=512,NZ_TOT=49
INTEGER, PARAMETER :: NZ=(NZ_TOT-1)/nproc + 1
INTEGER, PARAMETER :: NX2=3*NX/2,NY2=3*NY/2
INTEGER, PARAMETER :: LH=NX/2+1,LD=2*LH,LH_BIG=NX2/2+1,LD_BIG=2*LH_BIG
INTEGER :: JX,JY,JZ
CHARACTER(256) :: HEADERS,PHI_FILE,SRF_FILE,sscr1,sscr2
INTEGER :: NSTEPS,DT,JT, I
REAL(RPREC), DIMENSION(LD,NY,0:NZ) :: U, V, W
INTEGER, DIMENSION(1:12) :: REQ


CALL MPI_INIT (ierr)
CALL MPI_COMM_SIZE (MPI_COMM_WORLD, np, ierr)
CALL MPI_COMM_RANK (MPI_COMM_WORLD, rank, ierr)


!--set up a 1d cartesian topology 
!  CALL mpi_cart_create (MPI_COMM_WORLD, 1, (/ nproc /), (/ .false. /),  &
!  CALL mpi_cart_create (MPI_COMM_WORLD, 1, np, (/ .false. /),  &
!                        .true., comm, ierr)
!  CALL mpi_cart_shift (comm, 0, 1, down, up, ierr)
!  CALL mpi_cart_coords (comm, rank, 1, coords, ierr)
!  coord = coords(1)  !--use coord (NOT rank) to determine global position
!
!  write (chcoord, '(a,i0,a)') '(', coord, ')'  !--() make easier to use
!
!--rank->coord and coord->rank conversions
!  do ip = 0, np-1
!    CALL mpi_cart_rank (comm, (/ ip /), rank_of_coord(ip), ierr)
!    CALL mpi_cart_coords (comm, ip, 1, kkk, ierr)
!    coord_of_rank(ip) = kkk(1)
!  end do
!

! TO AVOID USAGE OF MPI CART FUNCTIONS
DOWN=rank-1
UP=rank+1
IF (rank .eq. np-1) THEN
UP=-1
ENDIF
write (*, *) 'np=', np, ' rank=', rank, ' up=', up, ' down=', down

! FOR NON BLOCKING SEND/RECV
DO I=1,12
REQ(I)=MPI_REQUEST_NULL  
ENDDO


MPI_RPREC = MPI_DOUBLE_PRECISION
COMM = MPI_COMM_WORLD

!---------------------------------------------------------------------------------
!INITIAL CONDITIONS
U=0.0_RPrec
V=0.0_RPrec
W=0.0_RPrec

!HORIZONTAL DERIVATIVES
U=U+1
V=V+1
W=W+1

NSTEPS=5000
T1= MPI_WTIME()

!MAIN LOOP
DO JT=1,NSTEPS


   CALL MPI_SENDRECV (U(1, 1, 1), LD*NY, MPI_RPREC, DOWN, 1,  &
                      U(1, 1, NZ), LD*NY, MPI_RPREC, UP, 1,   &
                      COMM, STATUS, IERR)
   CALL MPI_SENDRECV (V(1, 1, 1), LD*NY, MPI_RPREC, DOWN, 2,  &
                      V(1, 1, NZ), LD*NY, MPI_RPREC, UP, 2,   &
                      COMM, STATUS, IERR)
   CALL MPI_SENDRECV (W(1, 1, 1), LD*NY, MPI_RPREC, DOWN, 3,  &
                      W(1, 1, NZ), LD*NY, MPI_RPREC, UP, 3,   &
                      COMM, STATUS, IERR)                  
   CALL MPI_SENDRECV (U(1, 1, NZ-1), LD*NY, MPI_RPREC, UP, 4,  &
                      U(1, 1, 0), LD*NY, MPI_RPREC, DOWN, 4,   &
                      COMM, STATUS, IERR)
   CALL MPI_SENDRECV (V(1, 1, NZ-1), LD*NY, MPI_RPREC, UP, 5,  &
                      V(1, 1, 0), LD*NY, MPI_RPREC, DOWN, 5,   &
                      COMM, STATUS, IERR)
   CALL MPI_SENDRECV (W(1, 1, NZ-1), LD*NY, MPI_RPREC, UP, 6,  &
                      W(1, 1, 0), LD*NY, MPI_RPREC, DOWN, 6,   &
                      COMM, STATUS, IERR)


! SAME CODE BUT WITH MPI_IRECV AND MPI_ISEND
!
!   CALL MPI_IRECV (U(1, 1, NZ), LD*NY, MPI_RPREC, UP, 1, &
!                      COMM, REQ(3), IERR)
!   CALL MPI_IRECV (U(1, 1, 0), LD*NY, MPI_RPREC, DOWN, 4, &
!                      COMM, REQ(4), IERR)
!
!   CALL MPI_IRECV (V(1, 1, NZ), LD*NY, MPI_RPREC, UP, 2, &
!                      COMM, REQ(7), IERR)
!   CALL MPI_IRECV (V(1, 1, 0), LD*NY, MPI_RPREC, DOWN, 5, &
!                      COMM, REQ(8), IERR)
!
!   CALL MPI_IRECV (W(1, 1, NZ), LD*NY, MPI_RPREC, UP, 3, &
!                      COMM, REQ(11), IERR)
!   CALL MPI_IRECV (W(1, 1, 0), LD*NY, MPI_RPREC, DOWN, 6, &
!                      COMM, REQ(12), IERR)
!
!
!   CALL MPI_ISEND (U(1, 1, 1), LD*NY, MPI_RPREC, DOWN, 1, &
!                      COMM, REQ(1), IERR)
!   CALL MPI_ISEND (U(1, 1, NZ-1), LD*NY, MPI_RPREC, UP, 4, &
!                      COMM, REQ(2), IERR)
!
!   CALL MPI_ISEND (V(1, 1, 1), LD*NY, MPI_RPREC, DOWN, 2, &
!                      COMM, REQ(5), IERR)
!   CALL MPI_ISEND (V(1, 1, NZ-1), LD*NY, MPI_RPREC, UP, 5, &
!                      COMM, REQ(6), IERR)
!
!   CALL MPI_ISEND (W(1, 1, 1), LD*NY, MPI_RPREC, DOWN, 3, &
!                      COMM, REQ(9), IERR)
!   CALL MPI_ISEND (W(1, 1, NZ-1), LD*NY, MPI_RPREC, UP, 6, &
!                      COMM, REQ(10), IERR)
!
!
!   CALL MPI_WAITALL(12,REQ,STATUS,IERR)

ENDDO

T2= MPI_WTIME()
PRINT*,rank,' TOTAL LOOP TIME =', T2-T1

CALL MPI_FINALIZE(IERR)

END PROGRAM MAIN
