! COPYRIGHT (c) 2000 Council for the Central Laboratory
!               of the Research Councils
!
! Version 5.2.0
! See ChangeLog for history
!
module hsl_ma57_double
   use hsl_zd11_double
   implicit none
   integer, parameter, private :: wp = kind(0.0d0)

   type ma57_factors
     private
      integer, allocatable :: keep(:)
      integer, allocatable :: iw(:)
      real(wp), allocatable :: val(:)
      integer :: n        ! Matrix order
      integer :: nrltot   ! Size for val without compression
      integer :: nirtot   ! Size for iw without compression
      integer :: nrlnec   ! Size for val with compression
      integer :: nirnec   ! Size for iw with compression
      integer :: pivoting ! Set to pivoting option used in factorize
      integer :: scaling  ! Set to scaling option used in factorize
      integer :: static   ! Set to indicate if static pivots chosen
      integer :: rank     ! Set to indicate the rank of the factorization
      integer :: nirbdu   ! Set to number of integers in factors
      integer :: nebdu    ! Set to number of entries in factors
   end type ma57_factors

   type ma57_control
      real(wp) :: multiplier ! Factor by which arrays sizes are to be
                        ! increased if they are too small
      real(wp) :: reduce ! if previously allocated internal workspace arrays
                        !  are greater than reduce times the currently
                        !  required sizes, they are reset to current requirments
      real(wp) :: u     ! Pivot threshold
      real(wp) :: static_tolerance ! used for setting static pivot level
      real(wp) :: static_level ! used for switch to static
      real(wp) :: tolerance ! anything less than this is considered zero
      real(wp) :: convergence ! used to monitor convergence in iterative
                              ! refinement
      real(wp) :: consist ! used in test for consistency when using
                          ! fredholm alternative
      integer :: lp     ! Unit for error messages
      integer :: wp     ! Unit for warning messages
      integer :: mp     ! Unit for monitor output
      integer :: sp     ! Unit for statistical output
      integer :: ldiag  ! Controls level of diagnostic output
      integer :: nemin  ! Minimum number of eliminations in a step
      integer :: factorblocking ! Level 3 blocking in factorize
      integer :: solveblocking ! Level 2 and 3 blocking in solve
      integer :: la     ! Initial size for real array for the factors.
                        ! If less than nrlnec, default size used.
      integer :: liw    ! Initial size for integer array for the factors.
                        ! If less than nirnec, default size used.
      integer :: maxla  ! Max. size for real array for the factors.
      integer :: maxliw ! Max. size for integer array for the factors.
      integer :: pivoting  ! Controls pivoting:
                 !  1  Numerical pivoting will be performed.
                 !  2  No pivoting will be performed and an error exit will
                 !     occur immediately a pivot sign change is detected.
                 !  3  No pivoting will be performed and an error exit will
                 !     occur if a zero pivot is detected.
                 !  4  No pivoting is performed but pivots are changed to
                 !     all be positive.
      integer :: thresh ! Controls threshold for detecting full rows in analyse
                 !     Registered as percentage of N
                 ! 100 Only fully dense rows detected (default)
      integer :: ordering  ! Controls ordering:
                 !  Note that this is overridden by using optional parameter
                 !  perm in analyse call with equivalent action to 1.
                 !  0  AMD using MC47
                 !  1  User defined
                 !  2  AMD using MC50
                 !  3  Min deg as in MA57
                 !  4  Metis_nodend ordering
                 ! >4  Presently equivalent to 0 but may chnage
      integer :: scaling  ! Controls scaling:
                 !  0  No scaling
                 ! >0  Scaling using MC64 but may change for > 1
      integer :: rank_deficient  ! Controls handling rank deficiency:
                 !  0  No control
                 ! >0  Small entries removed during factorization

   end type ma57_control

   type ma57_ainfo
      real(wp) :: opsa = -1.0_wp  ! Anticipated # operations in assembly
      real(wp) :: opse = -1.0_wp  ! Anticipated # operations in elimination
      integer :: flag = 0     ! Flags success or failure case
      integer :: more = 0     ! More information on failure
      integer :: nsteps = -1  ! Number of elimination steps
      integer :: nrltot = -1  ! Size for a without compression
      integer :: nirtot = -1  ! Size for iw without compression
      integer :: nrlnec = -1  ! Size for a with compression
      integer :: nirnec = -1  ! Size for iw with compression
      integer :: nrladu = -1  ! Number of reals to hold factors
      integer :: niradu = -1  ! Number of integers to hold factors
      integer :: ncmpa = -1   ! Number of compresses
      integer :: ordering = -1! Indicates the ordering actually used
      integer :: oor = 0      ! Number of indices out-of-range
      integer :: dup = 0      ! Number of duplicates
      integer :: maxfrt = -1  ! Forecast maximum front size
      integer :: stat = 0     ! STAT value after allocate failure
   end type ma57_ainfo

   type ma57_finfo
      real(wp) :: opsa = -1.0_wp  ! Number of operations in assembly
      real(wp) :: opse = -1.0_wp  ! Number of operations in elimination
      real(wp) :: opsb = -1.0_wp  ! Additional number of operations for BLAS
      real(wp) :: maxchange = -1.0_wp
                  ! Largest pivot modification when pivoting=4
      real(wp) :: smin = -1.0_wp  ! Minimum scaling factor
      real(wp) :: smax = -1.0_wp  ! Maximum scaling factor
      integer :: flag = 0     ! Flags success or failure case
      integer :: more = 0     ! More information on failure
      integer :: maxfrt = -1  ! Largest front size
      integer :: nebdu = -1   ! Number of entries in factors
      integer :: nrlbdu = -1  ! Number of reals that hold factors
      integer :: nirbdu = -1  ! Number of integers that hold factors
      integer :: nrltot = -1  ! Size for a without compression
      integer :: nirtot = -1  ! Size for iw without compression
      integer :: nrlnec = -1  ! Size for a with compression
      integer :: nirnec = -1  ! Size for iw with compression
      integer :: ncmpbr = -1  ! Number of compresses of real data
      integer :: ncmpbi = -1  ! Number of compresses of integer data
      integer :: ntwo = -1    ! Number of 2x2 pivots
      integer :: neig = -1    ! Number of negative eigenvalues
      integer :: delay = -1   ! Number of delayed pivots (total)
      integer :: signc = -1   ! Number of pivot sign changes (pivoting=3)
      integer :: static = -1  ! Number of static pivots chosen
      integer :: modstep = -1 ! First pivot modification when pivoting=4
      integer :: rank = -1    ! Rank of original factorization
      integer :: stat = 0     ! STAT value after allocate failure
   end type ma57_finfo

   type ma57_sinfo
      real(wp) :: cond = -1.0_wp
          ! Condition number of matrix (category 1 equations)
      real(wp) :: cond2 = -1.0_wp
          ! Condition number of matrix (category 2 equations)
      real(wp) :: berr = -1.0_wp
          ! Condition number of matrix (category 1 equations)
      real(wp) :: berr2 = -1.0_wp
          ! Condition number of matrix (category 2 equations)
      real(wp) :: error = -1.0_wp
          ! Estimate of forward error using above data
      integer :: flag = 0    ! Flags success or failure case
      integer :: stat = 0    ! STAT value after allocate failure
   end type ma57_sinfo
   interface ma57_solve
! ma57_solve1 for 1 rhs  and ma57_solve2 for more than 1.
      module procedure ma57_solve1,ma57_solve2
   end interface

   interface ma57_part_solve
! ma57_part_solve1 for 1 rhs  and ma57_part_solve2 for more than 1.
      module procedure ma57_part_solve1,ma57_part_solve2
   end interface

contains

   subroutine ma57_initialize(factors,control)
      type(ma57_factors), intent(out), optional :: factors
      type(ma57_control), intent(out), optional :: control
      integer icntl(20),stat
      real(wp) cntl(5)
      if (present(factors)) then
        factors%n = 0
        deallocate(factors%keep,factors%val,factors%iw,stat=stat)
      end if
      if (present(control)) then
          call ma57id(cntl,icntl)
          control%u = cntl(1)
          control%tolerance = cntl(2)
          control%convergence = cntl(3)
          control%static_tolerance = cntl(4)
          control%static_level = cntl(5)
          control%lp = icntl(1)
          control%wp = icntl(2)
          control%mp = icntl(3)
          control%sp = icntl(4)
          control%ldiag = icntl(5)
          control%pivoting = icntl(7)
          control%ordering = icntl(6)
          control%scaling = icntl(15)
          control%factorblocking = icntl(11)
          control%nemin = icntl(12)
          control%solveblocking = icntl(13)
          control%thresh = icntl(14)
          control%rank_deficient = icntl(16)
          control%la = 0
          control%liw = 0
          control%maxla = huge(0)
          control%maxliw = huge(0)
          control%multiplier = 2.0
          control%reduce     = 2.0
          control%consist    = 1.0e-20_wp
      end if
    end subroutine ma57_initialize

   subroutine ma57_analyse(matrix,factors,control,ainfo,perm)
      type(zd11_type), intent(in) :: matrix
      type(ma57_factors), intent(inout) :: factors
      type(ma57_control), intent(in) :: control
      type(ma57_ainfo), intent(out) :: ainfo
      integer, intent(in), optional :: perm(matrix%n) ! Pivotal sequence

      integer, allocatable :: iw1(:)
      integer :: lkeep,n,ne,stat,icntl(20),info(40),rspace
      real(wp) rinfo(20)

      icntl(1) = control%lp
      icntl(2) = control%wp
      icntl(3) = control%mp
      icntl(4) = control%sp
      icntl(5) = control%ldiag
      icntl(6) = control%ordering
      icntl(12)  = control%nemin
      icntl(14) = control%thresh
! The next two are necessary so that the correct storage can be allocated
! No longer needed for internal purposes but icntl(7) is also used by
!     analyse to determine ordering if system is deemed positive definite
      icntl(7)  = control%pivoting
      icntl(15) = control%scaling
      n = matrix%n
      ne = matrix%ne
      stat = 0
      lkeep = 5*n+ne+max(n,ne)+42

      if(allocated(factors%keep)) then
         if(size(factors%keep)/=lkeep) then
            deallocate(factors%keep,stat=stat)
            if (stat/=0) go to 100
            allocate(factors%keep(lkeep),stat=stat)
            if (stat/=0) go to 100
         end if
      else
         allocate(factors%keep(lkeep),stat=stat)
         if (stat/=0) go to 100
      end if

! Override value of icntl(6) if perm present
      if (present(perm)) then
         factors%keep(1:n) = perm(1:n)
         icntl(6)=1
      end if

      allocate (iw1(5*n),stat=stat)
      if (stat/=0) go to 100

      call ma57ad(n,ne,matrix%row,matrix%col, &
               lkeep,factors%keep,iw1,icntl,info,rinfo)

! Adjust so that we don't allocate too much space in factorize
      rspace = 0
      if (control%pivoting == 4) rspace = rspace + n + 5
      if (control%scaling == 1)  rspace = rspace + n
      factors%n = n
      factors%nrltot = info(9)  - rspace
      factors%nirtot = info(10)
      factors%nrlnec = info(11) - rspace
      factors%nirnec = info(12)

      ainfo%opsa   = rinfo(1)
      ainfo%opse   = rinfo(2)
      ainfo%flag   = info(1)
      if (info(1) == -18) ainfo%flag   = -10
      ainfo%more   = info(2)
      ainfo%oor    = info(3)
      ainfo%dup    = info(4)
      ainfo%nrladu = info(5)
      ainfo%niradu = info(6)
      ainfo%maxfrt = info(7)
      ainfo%nsteps  = info(8)
      ainfo%nrltot = info(9)
      ainfo%nirtot = info(10)
      ainfo%nrlnec = info(11)
      ainfo%nirnec = info(12)
      ainfo%ncmpa  = info(13)
      ainfo%ordering = info(36)

      deallocate (iw1, stat=stat)
      if (stat/=0) go to 100
      return

  100 if (control%ldiag>0 .and. control%lp>0 ) &
          write (control%lp,'(/a/a,i5)') &
         'Error return from MA57_ANALYSE: flag = -3', &
         'Allocate or deallocate failed with STAT=',stat
       ainfo%flag = -3
       ainfo%stat = stat

   end subroutine ma57_analyse

   subroutine ma57_factorize(matrix,factors,control,finfo)
      type(zd11_type), intent(in) :: matrix
      type(ma57_factors), intent(inout) :: factors
      type(ma57_control), intent(in) :: control
      type(ma57_finfo), intent(out) :: finfo

      integer :: la,liw,lkeep,oldla,oldliw
      integer stat  ! stat value in allocate statements

      integer icntl(20),info(40),n,exla,expne
      real(wp) cntl(5),rinfo(20)

      integer, allocatable :: iwork(:)
      real(wp), allocatable :: temp(:)
      integer, allocatable :: itemp(:)
      integer :: hold,iwpos,istk,nfront,apos,astk

      n = matrix%n
      lkeep = 5*n+matrix%ne+max(n,matrix%ne)+42
      allocate (iwork(n),stat=stat)
      if (stat/=0) go to 100

      if(factors%n/=matrix%n) then
         if (control%ldiag>0 .and. control%lp>0 ) &
         write (control%lp,'(/a/a,i12,a,i12)') &
         'Error return from MA57_FACTORIZE: flag = -1', &
         'MATRIX%N has the value', &
         matrix%n,' instead of',factors%n
       finfo%flag = -1
       finfo%more = factors%n
       return
      end if

      cntl(1) = control%u
      cntl(2) = control%tolerance
      cntl(4) = control%static_tolerance
      cntl(5) =  control%static_level
      icntl(1) = control%lp
      icntl(2) = control%wp
      icntl(3) = control%mp
      icntl(4) = control%sp
      icntl(5) = control%ldiag
      icntl(7) = control%pivoting
      factors%pivoting = control%pivoting
      factors%scaling = control%scaling
      icntl(8) = 1
      icntl(11) = control%factorblocking
      icntl(12) = control%nemin
      icntl(15) = control%scaling
      icntl(16) =  control%rank_deficient
      stat = 0

      expne = factors%keep(matrix%n+2)

      la = control%la
      if(la<factors%nrlnec)then
         la = 0
         if(allocated(factors%val))la = size(factors%val)
         if(la>control%reduce*factors%nrltot) la = factors%nrltot
         if(la<factors%nrlnec) la = factors%nrltot
      end if

! Needed because explicitly removed from factors%nrlnec and factors%nrltot
         exla = 0
         if(control%pivoting == 4) exla =  factors%n + 5
         if(control%scaling == 1)  exla =  exla + factors%n
! Space for biga is already in la but not in exla
         la = max(la+exla,exla+expne+2)
         if(control%scaling == 1) la = max(la,exla+3*expne+3*factors%n+1)

      if(allocated(factors%val))then
         if(la/=size(factors%val))then
            deallocate(factors%val,stat=stat)
            if (stat/=0) go to 100
            allocate(factors%val(la),stat=stat)
            if (stat/=0) go to 100
         end if
      else
         allocate(factors%val(la),stat=stat)
         if (stat/=0) go to 100
      end if

      liw = control%liw
      if(liw<factors%nirnec)then
         liw = 0
         if(allocated(factors%iw))liw = size(factors%iw)
         if(liw>control%reduce*factors%nirnec) liw = factors%nirtot
         if(liw<factors%nirnec) liw = factors%nirtot
      end if

      if(control%scaling == 1) liw = max(liw,3*expne+5*factors%n+1)

      if(allocated(factors%iw))then
         if(liw/=size(factors%iw))then
            deallocate(factors%iw,stat=stat)
            if (stat/=0) go to 100
            allocate(factors%iw(liw),stat=stat)
            if (stat/=0) go to 100
         end if
      else
         allocate(factors%iw(liw),stat=stat)
         if (stat/=0) go to 100
      end if

      do
         call ma57bd(matrix%n,matrix%ne,matrix%val,factors%val, &
                  la,factors%iw,liw,lkeep,factors%keep,iwork,icntl,cntl, &
                  info,rinfo)
         finfo%flag = info(1)

         if (info(1)==11) then
            oldliw = liw
            liw = max(liw, int(control%multiplier*liw))
            if (liw>control%maxliw) then
               if (control%ldiag>0 .and. control%lp>0 ) &
                 write (control%lp,'(/a/a,i10)') &
                 'Error return from MA57_FACTORIZE: iflag = -8 ', &
                 'Main integer array needs to be bigger than', control%maxliw
               finfo%flag = -8
               return
            end if

! Expand integer space
            hold   = matrix%n + 3
            iwpos  = factors%keep(hold+7)
            istk   = factors%keep(hold+14) + 1
            nfront = factors%keep(hold+23)
            allocate (itemp(oldliw),stat=stat)
            if (stat/=0) go to 100
! Changed because middle space could be undefined
!           itemp(1:oldliw) = factors%iw(1:oldliw)
            itemp(1:iwpos+nfront-1) = factors%iw(1:iwpos+nfront-1)
            itemp(istk:oldliw)      = factors%iw(istk:oldliw)
            deallocate (factors%iw,stat=stat)
            if (stat/=0) go to 100
            allocate (factors%iw(liw),stat=stat)
            if (stat/=0) go to 100
            call ma57ed &
                  (matrix%n,1,factors%keep,factors%val,la,factors%val,la, &
                   itemp,oldliw,factors%iw,liw,info)
            deallocate (itemp,stat=stat)
            if (stat/=0) go to 100

         else if (info(1)==10) then
            oldla = la
            la = max(la, int(control%multiplier*la))
            if (la>control%maxla) then
               if (control%ldiag>0 .and. control%lp>0 ) &
                 write (control%lp,'(/a/a,i10)') &
                 'Error return from MA57_FACTORIZE: flag = -7 ', &
                 'Main real array needs to be bigger than', control%maxla
               finfo%flag = -7
               return
            end if

! Expand real space
            hold = matrix%n + 3
            apos = factors%keep(hold+9) - 1
            astk = factors%keep(hold+15) + 1
            allocate (temp(oldla),stat=stat)
            if (stat/=0) go to 100
! Changed because middle portion may be undefined
!           temp(1:oldla) = factors%val(1:oldla)
            temp(1:apos)     = factors%val(1:apos)
            temp(astk:oldla) = factors%val(astk:oldla)
            deallocate (factors%val,stat=stat)
            if (stat/=0) go to 100
            allocate (factors%val(la),stat=stat)
            if (stat/=0) go to 100
            call ma57ed &
                 (matrix%n,0,factors%keep,temp,oldla,factors%val,la, &
                 factors%iw,liw,factors%iw,liw,info)
            deallocate (temp,stat=stat)
            if (stat/=0) go to 100

         else
            exit

         end if

      end do

      deallocate (iwork,stat=stat)
      if (stat/=0) go to 100


      finfo%more = info(2)
      if (info(1)>=0) then
        finfo%nebdu  = info(14)
        finfo%nrlbdu = info(15)
        finfo%nirbdu = info(16)
        finfo%nrlnec = info(17)
        finfo%nirnec = info(18)
        finfo%nrltot = info(19)
        finfo%nirtot = info(20)
        finfo%maxfrt = info(21)
        finfo%ntwo   = info(22)
        finfo%delay  = info(23)
        finfo%neig   = info(24)
        finfo%rank   = info(25)
        finfo%signc  = info(26)
        finfo%static = info(35)
        factors%rank = info(25)
        factors%nirbdu = info(16)
        factors%nebdu = info(14)
        factors%static = 0
        if (finfo%static > 0) factors%static = 1
        finfo%modstep= info(27)
        finfo%ncmpbr = info(28)
        finfo%ncmpbi = info(29)
        finfo%opsa   = rinfo(3)
        finfo%opse   = rinfo(4)
        finfo%opsb   = rinfo(5)
        finfo%smin   = rinfo(16)
        finfo%smax   = rinfo(17)
        if (finfo%modstep > 0) finfo%maxchange = rinfo(14)
      end if

      return

  100 if (control%ldiag>0 .and. control%lp>0 ) &
         write (control%lp,'(/a/a,i5)') &
         'Error return from MA57_FACTORIZE: flag = -3', &
         'Allocate or deallocate failed with STAT=',stat
       finfo%flag = -3
       finfo%stat = stat

   end subroutine ma57_factorize

   subroutine ma57_solve2(matrix,factors,x,control,sinfo,rhs,iter,cond)
! Solve subroutine for multiple right-hand sides
      type(zd11_type), intent(in) :: matrix
      type(ma57_factors), intent(in) :: factors
      real(wp), intent(inout) :: x(:,:)
      type(ma57_control), intent(in) :: control
      type(ma57_sinfo), intent(out) :: sinfo
      real(wp), optional, intent(in) :: rhs(:,:)
      integer, optional, intent(in) :: iter
      integer, optional, intent(in) :: cond
      integer icntl(20),info(40),job,stat
      real(wp) cntl(5),rinfo(20),zero
      integer i,lw,n,nrhs

      integer, allocatable :: iwork(:)
      real(wp), allocatable :: work(:),resid(:),start(:,:)

      parameter (zero=0.0d0)
      n = matrix%n

      nrhs = size(x,2)

      stat = 0
      sinfo%flag = 0

! If rhs is present, then ma57dd will be called
! If iter is also present, then ADD algorithm is used
! If cond is also present, then condition number and error estimated

! lw is length of array work
      lw = n*nrhs
      if (present(rhs))  lw = n
      if (present(iter)) lw = 3*n
      if (factors%static == 1) lw = 3*n
      if (present(cond)) lw = 4*n

      allocate (iwork(n),work(lw),resid(n),stat=stat)
      if (stat/=0) go to 100

      if (factors%static == 1 .and. .not. present(rhs)) &
          allocate(start(n,nrhs),stat=stat)
      if (stat/=0) go to 100

      icntl(1) = control%lp
      icntl(2) = control%wp
      icntl(3) = control%mp
      icntl(4) = control%sp
      icntl(5) = control%ldiag
      icntl(13) = control%solveblocking
      icntl(15) = control%scaling

      cntl(3) = control%convergence

      if (present(iter)) then
! If iter is present, then rhs must be also, and user must set x.
        icntl(9)=100
        icntl(10)=0
        if (present(cond)) icntl(10)=1
        job = 2
        do i = 1,nrhs
        call ma57dd(job,matrix%n,matrix%ne,matrix%val,matrix%row,matrix%col, &
          factors%val,size(factors%val),factors%iw, &
          size(factors%iw),rhs(:,i), &
          x(:,i),resid,work,iwork,icntl,cntl,info,rinfo)
        enddo
        if (present(cond)) then
          sinfo%cond  = rinfo(11)
          sinfo%cond2 = rinfo(12)
          sinfo%berr  = rinfo(6)
          sinfo%berr2 = rinfo(7)
          sinfo%error = rinfo(13)
        endif
      else
        if(present(rhs)) then
          icntl(9) = 1
          icntl(10) = 0
          job = 2
          do i = 1,nrhs
          call ma57dd(job,matrix%n,matrix%ne,matrix%val,matrix%row, &
            matrix%col,factors%val,size(factors%val),factors%iw, &
            size(factors%iw),rhs(:,i), &
            x(:,i),resid,work,iwork,icntl,cntl,info,rinfo)
          enddo
        else
          if (factors%static == 1) then
            icntl(9) = 1
            icntl(10) = 0
            job = 2
            start = zero
            do i = 1,nrhs
            call ma57dd(job,matrix%n,matrix%ne,matrix%val,matrix%row, &
              matrix%col,factors%val,size(factors%val),factors%iw, &
              size(factors%iw),x(:,i), &
              start(:,i),resid,work,iwork,icntl,cntl,info,rinfo)
            enddo
            x = start
          else
            job=1
            call ma57cd(job,factors%n,factors%val,size(factors%val), &
              factors%iw,size(factors%iw),nrhs,x,size(x,1),   &
              work,nrhs*matrix%n,iwork,icntl,info)
          end if
        end if
      endif

      if (info(1) == -8 .or. info(1) == -14) sinfo%flag = -11

      deallocate (iwork,work,resid,stat=stat)
      if (factors%static == 1 .and. .not. present(rhs)) &
          deallocate (start,stat=stat)
      if (stat==0) return

  100 if (control%ldiag>0 .and. control%lp>0 )  &
          write (control%lp,'(/a/a,i5)') &
         'Error return from MA57_ANALYSE: flag = -3', &
         'Allocate or deallocate failed with STAT=',stat
      sinfo%flag = -3
      sinfo%stat = stat

   end subroutine ma57_solve2

   subroutine ma57_solve1(matrix,factors,x,control,sinfo,rhs,iter,cond)
      type(zd11_type), intent(in) :: matrix
      type(ma57_factors), intent(in) :: factors
      real(wp), intent(inout) :: x(:)
      type(ma57_control), intent(in) :: control
      type(ma57_sinfo), intent(out) :: sinfo
      real(wp), optional, intent(in) :: rhs(:)
      integer, optional, intent(in) :: iter
      integer, optional, intent(in) :: cond
      integer icntl(20),info(40),job,stat
      real(wp) cntl(5),rinfo(20),zero
      integer n,nrhs,lw

      integer, allocatable :: iwork(:)
      real(wp), allocatable :: work(:),resid(:),start(:)

      parameter (zero=0.0d0)

      n = matrix%n

      nrhs = 1

      stat = 0
      sinfo%flag = 0

! If rhs is present, then ma57dd will be called
! If iter is also present, then ADD algorithm is used
! If cond is also present, then condition number and error estimated

! lw is length of array work
      lw = n
      if (present(rhs))  lw = n
      if (present(iter)) lw = 3*n
      if (factors%static == 1) lw = 3*n
      if (present(cond)) lw = 4*n

      allocate (iwork(n),work(lw),resid(n),stat=stat)
      if (stat/=0) go to 100

      if (factors%static == 1 .and. .not. present(rhs)) &
          allocate(start(n),stat=stat)
      if (stat/=0) go to 100

      icntl(1) = control%lp
      icntl(2) = control%wp
      icntl(3) = control%mp
      icntl(4) = control%sp
      icntl(5) = control%ldiag
      icntl(13) = control%solveblocking
      icntl(15) = control%scaling

      cntl(3) = control%convergence

      stat = 0

      if (present(iter)) then
! If iter is present, then rhs must be also, and user must set x.
        icntl(9)=100
        icntl(10)=0
        if (present(cond)) icntl(10)=1
        job = 2
        call ma57dd(job,matrix%n,matrix%ne,matrix%val,matrix%row, &
          matrix%col,factors%val,size(factors%val),factors%iw, &
          size(factors%iw),rhs, &
          x,resid,work,iwork,icntl,cntl,info,rinfo)
        if (present(cond)) then
          sinfo%cond  = rinfo(11)
          sinfo%cond2 = rinfo(12)
          sinfo%berr  = rinfo(6)
          sinfo%berr2 = rinfo(7)
          sinfo%error = rinfo(13)
        endif
      else
        if(present(rhs)) then
          icntl(9) = 1
          icntl(10) = 0
          job = 2
          call ma57dd(job,matrix%n,matrix%ne,matrix%val,matrix%row, &
            matrix%col,factors%val,size(factors%val),factors%iw, &
            size(factors%iw),rhs, &
            x,resid,work,iwork,icntl,cntl,info,rinfo)
        else
          if (factors%static == 1) then
            icntl(9) = 1
            icntl(10) = 0
            job = 2
            start = zero
            call ma57dd(job,matrix%n,matrix%ne,matrix%val,matrix%row, &
              matrix%col,factors%val,size(factors%val),factors%iw, &
              size(factors%iw),x, &
              start,resid,work,iwork,icntl,cntl,info,rinfo)
            x = start
          else
            job=1
            call ma57cd(job,factors%n,factors%val,size(factors%val), &
              factors%iw,size(factors%iw),nrhs,x,size(x,1),   &
              work,nrhs*matrix%n,iwork,icntl,info)
          end if
        end if
      endif

      deallocate (iwork,work,resid,stat=stat)
      if (factors%static == 1 .and. .not. present(rhs))  &
          deallocate (start,stat=stat)
      if (stat==0) return

  100 if (control%ldiag>0 .and. control%lp>0 )  &
         write (control%lp,'(/a/a,i5)')  &
         'Error return from MA57_ANALYSE: flag = -3', &
         'Allocate or deallocate failed with STAT=',stat
      sinfo%flag = -3
      sinfo%stat = stat

   end subroutine ma57_solve1

   subroutine ma57_finalize(factors,control,info)
      type(ma57_factors), intent(inout) :: factors
      type(ma57_control), intent(in) :: control
      integer, intent(out) :: info
      integer :: inf

      info = 0
      inf = 0
      if (allocated(factors%keep)) deallocate(factors%keep,stat=inf)
      if (inf/=0) info = inf
      if (allocated(factors%iw)) deallocate(factors%iw,stat=inf)
      if (inf/=0) info = inf
      if (allocated(factors%val)) deallocate(factors%val,stat=inf)
      if (inf/=0) info = inf
      if (info==0) return

      if (control%ldiag>0 .and. control%lp>0 ) &
         write (control%lp,'(/a/a,i5)') &
         'Error return from MA57_FINALIZE:', &
         'Deallocate failed with STAT=',info

    end subroutine ma57_finalize

   subroutine ma57_enquire(factors,perm,pivots,d,perturbation,scaling)
      type(ma57_factors), intent(in) :: factors
      integer, intent(out), optional :: perm(factors%n),pivots(factors%n)
      real(wp), intent(out), optional :: d(2,factors%n)
      real(wp), intent(out), optional :: perturbation(factors%n)
      real(wp), intent(out), optional :: scaling(factors%n)

      real(wp) one,zero
      parameter (one=1.0d0,zero=0.0d0)

      integer block ! block index
      integer i     ! row index within the block
      integer ka    ! posn in array factors%val
      integer k2    ! posn in for off-diagonals of 2 x 2 pivots in factors%val
      integer kd    ! index in array d
      integer kp    ! index in array pivots
      integer kw    ! posn in array factors%iw
      integer ncols ! number of columns in the block
      integer nrows ! number of rows in the block
      logical two   ! flag for two by two pivots

      if(present(perm)) then
        perm = factors%keep(1:factors%n)
      endif

      if (present(perturbation)) then
        if (factors%pivoting == 4) then
          if (factors%scaling == 1) then
            perturbation = factors%val(size(factors%val) -  2*factors%n : &
                                       size(factors%val) -  factors%n -1)
          else
            perturbation = factors%val(size(factors%val) -  factors%n : &
                                       size(factors%val) - 1)
          endif
        else
          perturbation = zero
        end if
      endif


      if (present(scaling)) then
        if (factors%scaling == 1) then
          scaling = factors%val(size(factors%val) - factors%n : &
                                size(factors%val) - 1)
        else
          scaling = one
        end if
      endif

      if(present(pivots).or.present(d)) then
        ka = 1
        k2 = factors%iw(1)
        kd = 0
        kp = 0
        kw = 4
        if(present(d)) d = 0
          do block = 1, abs(factors%iw(3))
            ncols = factors%iw(kw)
            nrows = factors%iw(kw+1)
            if(present(pivots)) then
              pivots(kp+1:kp+nrows) = factors%iw(kw+2:kw+nrows+1)
              kp = kp + nrows
            end if
            if(present(d)) then
            two = .false.
            do i = 1,nrows
              kd = kd + 1
              d(1,kd) = factors%val(ka)
              if(factors%iw(kw+1+i)<0) two = .not.two
              if (two) then
                d(2,kd) = factors%val(k2)
                k2 = k2 + 1
              endif
              ka = ka + nrows + 1 - i
            end do
            ka = ka + nrows*(ncols-nrows)
          end if
          kw = kw + ncols + 2
        end do
      endif

   end subroutine ma57_enquire


   subroutine ma57_alter_d(factors,d,info)
      type(ma57_factors), intent(inout) :: factors
      real(wp), intent(in) :: d(2,factors%n)
      integer, intent(out) :: info

      integer block ! block index
      integer i     ! row index within the block
      integer ka    ! posn in array factors%val
      integer k2    ! posn in for off-diagonals of 2 x 2 pivots in factors%val
      integer kd    ! index in array d
      integer kw    ! posn in array factors%iw
      integer ncols ! number of columns in the block
      integer nrows ! number of rows in the block
      logical two   ! flag for two by two pivots

      info = 0
      ka = 1
      k2 = factors%iw(1)
      kd = 0
      kw = 4
      do block = 1, abs(factors%iw(3))
        ncols = factors%iw(kw)
        nrows = factors%iw(kw+1)
        two = .false.
        do i = 1,nrows
          kd = kd + 1
          factors%val(ka) = d(1,kd)
          if(factors%iw(kw+1+i)<0) two = .not.two
          if (two) then
            factors%val(k2) = d(2,kd)
            k2 = k2 + 1
          else
            if (d(2,kd) /= 0) info = kd
          end if
          ka = ka + nrows + 1 - i
        end do
        ka = ka + nrows*(ncols-nrows)
        kw = kw + ncols + 2
      end do

     end subroutine ma57_alter_d


   subroutine ma57_part_solve2(factors,control,part,x,info)
      type(ma57_factors), intent(in) :: factors
      type(ma57_control), intent(in) :: control
      character, intent(in) :: part
      real(wp), intent(inout) :: x(:,:)
      integer, intent(out) :: info
      integer icntl(20),n,nrhs,inf(40)

      integer, allocatable :: iwork(:)
      real(wp), allocatable :: work(:)

      n = factors%n

      nrhs = size(x,2)

      allocate (iwork(n),work(nrhs*n),stat=info)
      if (info/=0) go to 100

      icntl(1) = control%lp
      icntl(2) = control%wp
      icntl(3) = control%mp
      icntl(4) = control%sp
      icntl(5) = control%ldiag
      icntl(13) = control%solveblocking
      icntl(15) = control%scaling
      info = 0

      if(part=='L')then
         call ma57cd(2,factors%n,factors%val,size(factors%val),factors%iw, &
           size(factors%iw),nrhs,x,size(x,1),   &
           work,nrhs*factors%n,iwork,icntl,inf)
      else if(part=='D') then
         call ma57cd(3,factors%n,factors%val,size(factors%val),factors%iw, &
           size(factors%iw),nrhs,x,size(x,1),   &
           work,nrhs*factors%n,iwork,icntl,inf)
      else if(part=='U') then
         call ma57cd(4,factors%n,factors%val,size(factors%val),factors%iw, &
           size(factors%iw),nrhs,x,size(x,1),   &
           work,nrhs*factors%n,iwork,icntl,inf)
      end if


      deallocate (iwork,work,stat=info)
      if (info==0) return

  100 if (control%ldiag>0 .and. control%lp>0 )  &
          write (control%lp,'(/a/a,i5)') &
         'Error return from MA57_ANALYSE: flag = -3', &
         'Allocate or deallocate failed with STAT=',info

end subroutine ma57_part_solve2

   subroutine ma57_part_solve1(factors,control,part,x,info)
      type(ma57_factors), intent(in) :: factors
      type(ma57_control), intent(in) :: control
      character, intent(in) :: part
      real(wp), intent(inout) :: x(:)
      integer, intent(out) :: info
      integer inf(40),icntl(20),n,nrhs

      integer, allocatable :: iwork(:)
      real(wp), allocatable :: work(:)

      n = factors%n

      nrhs = 1

      allocate (iwork(n),work(n),stat=info)
      if (info/=0) go to 100

      icntl(1) = control%lp
      icntl(2) = control%wp
      icntl(3) = control%mp
      icntl(4) = control%sp
      icntl(5) = control%ldiag
      icntl(13) = control%solveblocking
      icntl(15) = control%scaling
      info = 0

      if(part=='L')then
         call ma57cd(2,factors%n,factors%val,size(factors%val),factors%iw, &
           size(factors%iw),nrhs,x,size(x,1),   &
           work,nrhs*factors%n,iwork,icntl,inf)
      else if(part=='D') then
         call ma57cd(3,factors%n,factors%val,size(factors%val),factors%iw, &
           size(factors%iw),nrhs,x,size(x,1),   &
           work,nrhs*factors%n,iwork,icntl,inf)
      else if(part=='U') then
         call ma57cd(4,factors%n,factors%val,size(factors%val),factors%iw, &
           size(factors%iw),nrhs,x,size(x,1),   &
           work,nrhs*factors%n,iwork,icntl,inf)
      end if


      deallocate (iwork,work,stat=info)
      if (info==0) return

  100 if (control%ldiag>0 .and. control%lp>0 ) &
          write (control%lp,'(/a/a,i5)') &
         'Error return from MA57_ANALYSE: flag = -3', &
         'Allocate or deallocate failed with STAT=',info

end subroutine ma57_part_solve1

   subroutine ma57_fredholm_alternative(factors,control,x,fredx,sinfo)
! This subroutine either returns the solution of the system if it is consistent
! or it returns a vector y in the null space of A that is non-orthogonal to the
! right-hand side. That is, Ay = 0 and y^T rhs ne 0.
      type(ma57_factors), intent(in) :: factors
      type(ma57_control), intent(in) :: control
! x     Right-hand side vector.  On exit it is solution if system consistent
!       otherwise it is solution of nonsingular part.
! fredx returns a Fredholm vector, y, if system is inconsistent.
!sinfo  is used for error return from allocate/deallocate and also to flag
!       whether system is consistent.
      real(wp), intent(inout) :: x(factors%n)
      real(wp), intent(out) :: fredx(factors%n)
      type(ma57_sinfo), intent(out) :: sinfo

! Local variables
      integer i,inf(40),info,icntl(20),k,n,rank,nrhs

! Work arrays
      integer, allocatable :: iwork(:),singular_rows(:)
      real(wp), allocatable :: work(:)

! Constants
      real(wp) zero
      parameter (zero=0.0d0)

      sinfo%flag = 0
      info = 0

      n = factors%n
      rank = factors%rank

      nrhs = 1

      allocate (iwork(n),work(n),singular_rows(n-rank),stat=info)
      if (info/=0) go to 100

      icntl(1) = control%lp
      icntl(2) = control%wp
      icntl(3) = control%mp
      icntl(4) = control%sp
      icntl(5) = control%ldiag
      icntl(13) = control%solveblocking
      icntl(15) = control%scaling

! Compute indices of singular rows
      k = factors%nirbdu
      do i = 1, n-rank
        singular_rows(i) = factors%iw(k)
        k = k - 3
      enddo

! Perform forward elimination

      call ma57cd(2,factors%n,factors%val,size(factors%val),factors%iw, &
                  size(factors%iw),nrhs,x,size(x),   &
                  work,nrhs*factors%n,iwork,icntl,inf)

! Check to see if system was consistent
      do k = 1,n-rank
        i = singular_rows(k)
        if (abs(x(i)).gt.control%consist) then
          sinfo%flag = 1
        endif
      enddo

      if (sinfo%flag.eq.1) then
! System is inconsistent
! Set up right-hand side for Fredholm vector
        fredx = zero
        do k = 1,n-rank
          i = singular_rows(k)
          fredx(i) = x(i)
        enddo
      endif

! Return solution.  In inconsistent case, this will be a minimum norm
! solution in the LL^T norm.

      call ma57cd(3,factors%n,factors%val,size(factors%val),factors%iw, &
                  size(factors%iw),nrhs,x,size(x),   &
                  work,nrhs*factors%n,iwork,icntl,inf)

      call ma57cd(4,factors%n,factors%val,size(factors%val),factors%iw, &
                  size(factors%iw),nrhs,x,size(x),   &
                  work,nrhs*factors%n,iwork,icntl,inf)

      if (sinfo%flag.eq.0) go to 50

! System is inconsistent.  Compute Fredholm vector
! Do not need to solve for D because this part of rhs is all zero.
! Solve for L transpose

      call ma57cd(4,factors%n,factors%val,size(factors%val),factors%iw, &
                  size(factors%iw),nrhs,fredx,size(fredx),   &
                  work,nrhs*factors%n,iwork,icntl,inf)

   50 deallocate (iwork,work,singular_rows,stat=info)
      if (info==0) return

  100 sinfo%flag = -3
      sinfo%stat = info
      if (control%ldiag>0 .and. control%lp>0 ) &
      write (control%lp,'(/a/a,i5)') &
      'Error return from ma57_fredholm_alternative', &
      'Allocate or deallocate failed with STAT=',info
      return

end subroutine ma57_fredholm_alternative

subroutine ma57_get_factors(factors,control,nzl,iptrl,lrows,lvals, &
                   nzd,iptrd,drows,dvals,perm,invperm,scale,sinfo)
! This subroutine returns the factors L and D in a standard format (column
!     pointer, row index).
      type(ma57_factors), intent(in) :: factors
      type(ma57_control), intent(in) :: control
! L is in arrays iptrl(n+1), lrows(nzl), lvals(nzl)
! D is in arrays iptrd(n+1), drows(nzd), dvals(nzd)
! The scaling vector is returned in scale(n)
! The permutation vectors in perm(n) and invperm(n)
      real(wp), intent(out) :: lvals(factors%nebdu),dvals(2*factors%n), &
                               scale(factors%n)
      integer, intent(out) :: nzl,nzd,iptrl(factors%n+1),lrows(factors%nebdu), &
                              iptrd(factors%n+1),drows(2*factors%n),  &
                              perm(factors%n),invperm(factors%n)
      type(ma57_sinfo), intent(out) :: sinfo

! Local variables
! i is a do loop index
! info returns possible error flag from allocate/deallocate failure
! iscal is used to point to scaling vector in factors%val
! n is set to order of matrix
! rank is set to computed rank of matrix
      integer i,info,iscal,n,rank

! Set constants
      real(wp) zero, one
      parameter (zero=0.0d0, one=1.0d0)

      sinfo%flag = 0

      n = factors%n
      rank = factors%rank

! Record scaling vector in array scale
      if (control%scaling.eq.1) then
        iscal = size(factors%val) - n
        DO i = 1,n
          scale(i) = factors%val(iscal + i -1)
        ENDDO
      else
        DO i = 1,N
          scale(i) = one
        ENDDO
      endif


      call ma57lfd(n,factors%val,size(factors%val),factors%iw,size(factors%iw),&
                   factors%nebdu,nzl,iptrl,lrows,lvals,nzd,iptrd,drows,dvals,  &
                   invperm,perm,rank,info)

      if (info.ne.0) then
        sinfo%flag = -3
        sinfo%stat = info
        if(control%ldiag>0 .and. control%lp>0 ) &
           write (control%lp,'(/a/a,i5)') &
          'Error return from ma57_get_factors', &
          'Allocate or deallocate failed with STAT=',info
      endif

      return

contains

      subroutine ma57lfd(n,fact,lfact,ifact,lifact,nebdu,nnz,ipl,irn,fl,  &
                        nnzd,ipd,id,d,ivp,iperm,rank,info)
! This subroutine extracts the factors stored in FACT/IFACT by MA57BD.
!
! n      must be set to the order of the matrix. It is not altered.
! fact   must be set to hold the real values corresponding to the
!        factors. This must be unchanged since the preceding call to
!        ma57_factorize or ma57_alter_d. It is not altered.
! lfact  length of array fact. It is not altered.
! ifact  holds the integer indexing information for the matrix factors
!        in fact. This must be unchanged since the preceding call to
!        ma57_factorize or ma57_alter_d. It is not altered.
! lifact length of array ifact. It is not altered.
! nnz    is the number of nonzeros in L
! ipl    on output will hold pointers to the start of the row entries
!        held by columns. The entries in JCN between ipl(j) and ipl(j+1)-1 are
!        the the row indices of column j in the L factor.
!        Integer array length N+1.
! irn    on output, will hold the row indices of the matrix L.
! fl     on output, will hold the values of the matrix L.
! nnzd   is the number of nonzeros in D (<2*N).
! ipd    on output will hold pointers to the start of the row entries
!        held by columns. The entries in id between ipd(j) and ipd(j+1)-1 are
!        the the row indices of column j in the D factor.
!        Integer array length N+1.
! d      on output, will hold the value of the pivots (1x1 and 2x2)
! id     on output, will hold the row coordinates of D
! ivp    need not be set on entry. On exit IVP(IPERM(I)) = I
! iperm  on output holds the permutation of the matrix A
! rank   is the calculated rank of the matrix
! info   is an integer flag to transmit any error in allocate/deallocate
      integer n,lfact,lifact,nebdu,nnz,nnzd
      real (wp) fact(lfact),FL(nebdu),D(2*n)
      integer ifaCt(lifact),irn(nebdu),id(2*n)
      integer ipl(n+1),ipd(n+1)
      integer iperm(n),ivp(n),rank,info
!
! icl59  is set to the control parameters of MC59.
! info59 error array for MC59.
!
      integer icl59(10),info59(10)
      EXTERNAL MC59AD

! Work arrays
! jcn will hold the column indices for entry to MC59 for the L factor
! jd  will hold the column indices for entry to MC59 for the D factor
         integer, allocatable :: jcn(:),jd(:),iw(:)

! Procedures
      INTRINSIC ABS

! Local variables
! apos  Current position in array FACT
! i     Temporary DO index
! iblk  Index of block pivot
! ii    Temporary index
! ip    Temporary index
! iwpos Position in IFACT of start of current index list
! j     Temporary DO index
! k     Temporary pointer to position in real array
! ncols Number of columns in the block pivot
! nrows Number of rows in the block pivot
! ind   Is an offset for indices
! kd    Counter variable
! apos2 Position in array FACT of D factors
! l     Temporary index
! c     Temporary variable
! temp  Temporary variable
!
      integer apos,i,iblk,ii,ip,iwpos,j,k,ncols,nrows,ind,kd,apos2,l
      real (wp) zero, one, c, temp
      parameter (zero=0.0d0,one=1.0d0)

      info = 0

!
! Compute IPERM and initialize the pointers

!
      iperm = 0
      APOS = 1
      apos2 = 1
      ind = 0
      iwpos = 4
      k = 0
! Run through each block in turn
      DO 270 iblk = 1,ifact(3)
! Find the number of rows and columns in the block
         ncols = ifact(iwpos)
         nrows = ifact(iwpos+1)
         iwpos = iwpos + 2
! Build the permutation
         do I = 1,nrows
            ii = ifact(iwpos+i-1)
            iperm(i+ind) = ii
         enddo
         iwpos = iwpos + ncols
         ind = ind + nrows
         apos2 = apos2 + (nrows*(nrows+1))/2 + nrows*(ncols-nrows)
 270  CONTINUE
! Build the inverse permutation
      do i = 1,N
         ivp(abs(iperm(I))) = i
      enddo
!
! Check if the rank = 0
!
      if (rank .EQ. 0) then
! Rank = 0 ==> D = 0 L = Identity
         do i=1,n
           irn(i) = i
           ipl(i) = i
           fl(i)  = one
           ipd(i) = i
           id(i)  = i
           d(i)   = zero
         enddo
         ipl(n+1) = n+1
         ipd(n+1) = n+1
         nnzd = n
         nnz  = n
         return
      endif

      allocate(jcn(nebdu),jd(2*n),iw(n+1),stat=info)
      if (info.ne.0) return

         apos = 1
         ind = 0
         iwpos = 4
         k = 1
         kd = 1
         DO 370 iblk = 1,ifact(3)
! Find the number of rows and columns in the block
           ncols = ifact(iwpos)
           nrows = ifact(iwpos+1)
           iwpos = iwpos + 2
! Treat diagonal block and build L and D^(-1)
           DO 30 j = 1,nrows
             DO 40 i = j,nrows
                if (i.eq.j) then
                  fl(k) = one
                  jcn(k) = j + ind
                  irn(k) = j + ind
!CC Why stop when singular
                  if ((i+ind).le.rank) then
                    d(kd) = fact(apos)
                    id(kd) = i + ind
                    jd(kd) = i + ind
                    kd = kd + 1
                  endif
! Build the matrix D^(-1) structure as a nonsymmetric matrix
                  if (iperm(i+ind).lt.0) then
! Two by two pivot
                    d(kd) = fact(apos2)
                    id(kd) = i + ind
                    jd(kd) = i + 1 + ind
                    d(kd+1) = fact(apos2)
                    id(kd+1) = jd(kd)
                    jd(kd+1) = id(kd)
                    kd = kd + 2
                    apos2 = apos2 + 1
                    iperm(i+ind) = abs(iperm(i+ind))
                    iperm(i+1+ind) = abs(iperm(i+1+ind))
                  endif
                else
! Save entries of L
                  fl(k) = fact(apos)
                  jcn(k) = j + ind
                  irn(k) = i + ind
                endif
                apos = apos + 1
                k = k + 1
 40         CONTINUE
 30       CONTINUE

! Treat off-diagonal block
           if (ncols.gt.nrows) then
             do 60 i = 1,nrows
               do 50 j = 1,ncols-nrows
                 ii = ifact(iwpos+j-1+nrows)
                 ip = ivp(abs(ii))
                 fl(k) = -fact(apos)
                 jcn(k) = i + ind
                 irn(k) = abs(ip)
                 apos = apos + 1
                 k = k + 1
 50            continue
 60          continue
           ENDIF
           iwpos = iwpos + ncols
           ind = ind + nrows
  370    CONTINUE
         nnz = k - 1
         nnzd = kd - 1

! Now sort entries
         icl59 = 0
         icl59(1) = 1
         icl59(2) = 1
         icl59(4) = -1
         icl59(5) = -1

         call mc59ad(icl59,n,n,nnz,irn,nnz,jcn,nnz,fl,n+1,ipl,n+1, &
                    iw,info59)
         if (info59(1) .eq. 0) then
            call mc59ad(icl59,n,n,nnzd,id,nnzd,jd,nnzd,d,n+1,ipd,n+1, &
                       iw,info59)
         endif

         if (info59(1) .ge. 0) then
!
!      Inversion of D
            ii = 0
            DO 580 k=1,n
               if ( ii .eq. 0 ) then
                 l = ipd(k+1) - ipd(k)
                 if (l .ne. 0) then
                   if (l .eq. 1) then
! Inverting 1 by 1 pivot
                     d(ipd(k)) = one / d(ipd(k))
                     ii = 0
                   else
! Inverting 2 by 2 pivot
                      i  = ipd(k)
                      c  = d(i)*d(i+3)-d(i+1)*d(i+2)
                      temp = d(i)/c
                      d(i) = d(i+3)/c
                      d(i+3)= temp
                      d(i+1) = -d(i+1)/c
                      d(i+2)= d(i+1)
                      ii = 1
                   endif
                 else
                   ii = 0
                 endif
               else
                  ii = 0
               endif
 580       CONTINUE
! Complete D if matrix is rank deficient
           k = nnzd
           DO i = rank+1,n
             k = k + 1
             ipd(i) = k
             d(k) = zero
             id(k) = I
           enddo
           nnzd = k
           ipd(n+1) = nnzd+1
!
!      Finished inversion of d
!
! I don't see how we can ever execute this.
! Would need out-of-range or duplicates to MC59.
         else
           do i=1,n
              ipd(i) = i
              id(i)  = i
              d(i)   = zero
           enddo
           ipd(n+1) = n+1
         endif
      deallocate(jcn,jd,iw,stat=info)
      return

      end subroutine MA57LFD

end subroutine ma57_get_factors

   subroutine ma57_sparse_lsolve(factors,control,nzrhs,irhs,nzsoln,isoln,  &
                                 rhs,sinfo)
! This subroutine performs forward elimination on a sparse right-hand side
      type(ma57_factors), intent(in) :: factors
      type(ma57_control), intent(in) :: control
      integer, intent(in) :: nzrhs
      integer, intent(in) :: irhs(nzrhs)
      integer, intent(out) :: nzsoln,isoln(*)
      real(wp), intent(inout) :: rhs(factors%n)
      type(ma57_sinfo), intent(out) :: sinfo

! Loacl variables
      integer blcntl,info,n,node,nodes,nstk

      integer, allocatable :: tree(:),flag(:)
      real(wp), allocatable :: scaling(:)

      n = factors%n

! Used in partitioning factors%keep array
      node = 2*n+43
      nstk = 3*n+43

      nodes = factors%iw(3)

      info = 0
      sinfo%flag = 0

      allocate (tree(nodes),flag(nodes),stat=info)
      if (info/=0) go to 100

! Generate tree pointers from children to parents
      call ma57_regen(nodes,factors%keep(nstk),tree,info)
      if (info.ne.0) go to 100

! Set flags for particular rhs
      call ma57_nflag(n,nodes,nzrhs,factors%keep(node),irhs,tree,flag)
      if (factors%scaling == 1) then
          allocate (scaling(n),stat= info)
          scaling = factors%val(size(factors%val) - factors%n : &
                                size(factors%val) - 1)
          rhs = scaling*rhs
      endif

      blcntl = control%solveblocking

! Call L solve
      call lsolve(n,factors%val,size(factors%val), factors%iw, &
           size(factors%iw),rhs,isoln,nzsoln,flag,blcntl,info)
      if (info.ne.0) go to 100
      call dsolve(n,factors%val,size(factors%val), &
           factors%iw,size(factors%iw), &
           rhs,flag,blcntl,info)
      if (info.ne.0) go to 100

      deallocate (flag,tree,stat=info)
      if (info==0) return

  100 sinfo%stat = info
      sinfo%flag = -3
      if (control%ldiag>0 .and. control%lp>0 ) &
          write (control%lp,'(/a/a,i5)') &
         'Error return from ma57_sparse_lsolve', &
         'Allocate or deallocate failed with STAT=',info
      return

contains
!CCCC
!BBBB

      subroutine lsolve(N,FACT,LFACT,IFACT,LIFACT,RHS,&
                        ISOLN,NZSOL,flag,blascntl,info)
! This subroutine performs forward elimination using the factors
!     stored in FACT/IFACT by MA57BD.
! It is designed for efficiency on a sparse right-hand side.
      INTEGER N,LFACT,NZSOL
      real(wp) FACT(LFACT)
      INTEGER LIFACT,IFACT(LIFACT)
      real(wp) RHS(N)
      INTEGER blascntl,info
      INTEGER FLAG(*),ISOLN(*)
      real (wp), allocatable :: w(:)
! N   must be set to the order of the matrix. It is not altered.
! FACT   must be set to hold the real values corresponding to the
!     factors. This must be unchanged since the preceding call to
!     MA57BD. It is not altered.
! LFACT  length of array FACT. It is not altered.
! IFACT  holds the integer indexing information for the matrix factors
!     in FACT. This must be unchanged since the preceding call to
!     MA57BD. It is not altered.
! LIFACT length of array IFACT. It is not altered.
! RHS on input, must be set to hold the right hand side vector.  On
!     return, it will hold the modified vector following forward
!     elimination.
! W   used as workspace to hold the components of the right hand
!     sides corresponding to current block pivotal rows.
!

! Procedures
      INTRINSIC ABS
      EXTERNAL DGEMV,DTPSV

! Constant
      real(wp) ONE
      PARAMETER (ONE=1.0D0)
!
! Local variables
! APOS  Current position in array FACT.
! I     Temporary DO index
! IBLK  Index of block pivot.
! II    Temporary index.
! IPIV  Pivot index.
! IRHS  RHS index.
! IWPOS Position in IFACT of start of current index list.
! J     Temporary DO index
! K     Temporary pointer to position in real array.
! J1    Position in IFACT of index of leading entry of row.
! J2    Position in IFACT of index of trailing entry of row.
! NCOLS Number of columns in the block pivot.
! NROWS Number of rows in the block pivot.
! W1    RHS value.
! W2    RHS value.
      INTEGER APOS,I,IBLK,II,IPIV,IRHS,IWPOS,J,J1,J2,K,K1,K2,NCOLS,NROWS
      real(wp) W1,W2

      allocate (w(n),stat=info)
      if (info.ne.0) return

      APOS = 1
      IWPOS = 4
      NZSOL = 0

      DO 270 IBLK = 1,IFACT(3)

! Find the number of rows and columns in the block.
        NCOLS = IFACT(IWPOS)
        NROWS = IFACT(IWPOS+1)

        IF (FLAG(IBLK) == 0) THEN
! Must set APOS and IWPOS for next block
          IWPOS = IWPOS + NCOLS + 2
          APOS = APOS + (NROWS* (NROWS+1))/2
          APOS = APOS + NROWS* (NCOLS-NROWS)
          GO TO 270
        ENDIF

        IWPOS = IWPOS + 2

        IF (NROWS.GT.4 .AND. NCOLS.GT.blascntl) THEN

! Perform operations using direct addressing.

! Load appropriate components of right-hand sides into W.
          DO 10 I = 1,NCOLS
            II = ABS(IFACT(IWPOS+I-1))
            IF (I.LE.NROWS) THEN
              NZSOL = NZSOL + 1
              ISOLN(NZSOL) = II
            ENDIF
            W(I) = RHS(II)
   10     CONTINUE


! Treat diagonal block (direct addressing)
          CALL DTPSV('L','N','U',NROWS,FACT(APOS),W,1)
          APOS = APOS + (NROWS* (NROWS+1))/2

! Treat off-diagonal block (direct addressing)
          IF (NCOLS.GT.NROWS) CALL DGEMV('N',NCOLS-NROWS,NROWS,  &
                                        ONE,FACT(APOS),NCOLS-NROWS,  &
                                        W,1,ONE,W(NROWS+1),1)
          APOS = APOS + NROWS* (NCOLS-NROWS)

! Reload W back into RHS.
          DO 35 I = 1,NCOLS
            II = ABS(IFACT(IWPOS+I-1))
            RHS(II) = W(I)
   35     CONTINUE

        ELSE

! Perform operations using indirect addressing.

        J1 = IWPOS
        J2 = IWPOS + NROWS - 1


! Treat diagonal block (indirect addressing)
        DO 130 IPIV = 1,NROWS
          APOS = APOS + 1
          I = ABS(IFACT(J1))
          NZSOL = NZSOL + 1
          ISOLN(NZSOL) = I
          W1 = RHS(I)
          K = APOS
          DO 100 J = J1+1,J2
            IRHS = ABS(IFACT(J))
            RHS(IRHS) = RHS(IRHS) - FACT(K)*W1
            K = K + 1
  100     CONTINUE
          APOS = K
          J1 = J1 + 1
  130   CONTINUE

! Loop unrolling
        J2 = IWPOS + NCOLS - 1
        DO 136 IPIV = 1,NROWS-1,2
          K1 = APOS
          K2 = APOS+NCOLS-NROWS
          W1 = RHS(ABS(IFACT(IWPOS+IPIV-1)))
          W2 = RHS(ABS(IFACT(IWPOS+IPIV)))
          DO 133 J = J1,J2
            IRHS = ABS(IFACT(J))
            RHS(IRHS) = RHS(IRHS) + W1*FACT(K1) + W2*FACT(K2)
            K1 = K1 + 1
            K2 = K2 + 1
  133     CONTINUE
          APOS = K2
  136   CONTINUE

        IF (MOD(NROWS,2).EQ.1) THEN
          K = APOS
          W1 = RHS(ABS(IFACT(IWPOS+IPIV-1)))
          DO 137 J = J1,J2
            IRHS = ABS(IFACT(J))
            RHS(IRHS) = RHS(IRHS) + W1*FACT(K)
            K = K + 1
  137     CONTINUE
          APOS = K
        ENDIF
      END IF

      IWPOS = IWPOS + NCOLS
  270 CONTINUE

      deallocate(w,stat=info)
      return

end subroutine lsolve

end subroutine ma57_sparse_lsolve

   subroutine ma57_sparse_dsolve(factors,control,nzrhs,irhs,rhs,sinfo)
      type(ma57_factors), intent(in) :: factors
      type(ma57_control), intent(in) :: control
      real(wp), intent(inout) :: rhs(factors%n)
      integer, intent(in) :: nzrhs
      integer, intent(in) :: irhs(nzrhs)
      type(ma57_sinfo), intent(out) :: sinfo

      integer, allocatable :: tree(:),flag(:)
      real(wp), allocatable :: scaling(:)

      integer blcntl,info,n,node,nodes,nstk

      n = factors%n

! Used in partitioning factors%keep array
      node = 2*n+43
      nstk = 3*n+43

      nodes = factors%iw(3)

      info = 0
      sinfo%flag = 0

      allocate (tree(nodes),flag(nodes),stat=info)
      if (info/=0) go to 100

! Generate tree pointers from children to parents
      call ma57_regen(nodes,factors%keep(nstk),tree,info)
      if (info.ne.0) go to 100

! Set flags for particular rhs
      call ma57_nflag(n,nodes,nzrhs,factors%keep(node),irhs,tree,flag)
      if (factors%scaling == 1) then
          allocate (scaling(n),stat= info)
          if (info.ne.0) go to 100
          scaling = factors%val(size(factors%val) - factors%n : &
                                size(factors%val) - 1)
!         rhs = scaling*rhs
      endif


      blcntl = control%solveblocking

! Call dsolve
      call dsolve(n,factors%val,size(factors%val), &
           factors%iw,size(factors%iw), &
           rhs,flag,blcntl,info)
      if (info.ne.0) return
       if (factors%scaling == 1) then
          rhs = rhs*scaling
          deallocate (scaling,stat=info)
          if (info.ne.0) go to 100
      endif

      deallocate (flag,tree,stat=info)
      if (info==0) return

  100 sinfo%stat = info
      sinfo%flag = -3
      if (control%ldiag>0 .and. control%lp>0 ) &
          write (control%lp,'(/a/a,i5)') &
         'Error return from MA57_SPARSE_LSOLVE: flag = -3', &
         'Allocate or deallocate failed with STAT=',info

!contains
!CCCCC
!BBBB
end subroutine ma57_sparse_dsolve

   subroutine ma57_lmultiply(factors,control,trans,x,y,sinfo)
! This subroutine multiplies a vetor (x) by L or L transpose.  The
! result is the vector y.
! If trans is 'N' (or 'n') then product is: S^(-1) P^T L S x otherwise
! it is: S L^T P^T S^(-1) x.
      type(ma57_factors), intent(in) :: factors
      type(ma57_control), intent(in) :: control
      character, intent(in) :: trans
      real(wp), intent(in) :: x(:)
      real(wp), intent(out) :: y(:)
      type(ma57_sinfo), intent(out) :: sinfo

! Work array is only used if scaling is on
      real(wp), allocatable :: work(:)

! Local variables
! i is do loop variable
! iscale is position of scaling in array factors%val
! n is local variable holding order of L
      integer i,info,iscale,n

      sinfo%flag = 0
      info = 0
      n = factors%n
      if (control%scaling.eq.1) then
        allocate (work(n),stat=info)
        if (info.ne.0) go to 100
      endif

      if (trans.eq.'N' .or. trans.eq.'n') then
! Multiply by L
! Scale input vector
        if (control%scaling.eq.1) then
          iscale = size(factors%val)-n
          do i = 1,n
            work(i) = x(i)*factors%val(iscale+i-1)
          enddo
          call lmult(factors%n,factors%val,size(factors%val),  &
                     factors%iw,size(factors%iw),work,y)
        else
          call lmult(factors%n,factors%val,size(factors%val),  &
                     factors%iw,size(factors%iw),x,y)
        endif
! Scale result
        if (control%scaling.eq.1) then
          do i = 1,n
            y(i) = y(i)/factors%val(iscale+i-1)
          enddo
        endif

      else
! Multiply by L transpose
! Scale input vector
        if (control%scaling.eq.1) then
          iscale = size(factors%val)-n
          do i = 1,n
            work(i) = x(i)/factors%val(iscale+i-1)
          enddo
          call ltmult(factors%n,factors%val,size(factors%val),  &
                     factors%iw,size(factors%iw),work,y)
        else
          call ltmult(factors%n,factors%val,size(factors%val),  &
                     factors%iw,size(factors%iw),x,y)
        endif
! Scale result
        if (control%scaling.eq.1) then
          do i = 1,n
            y(i) = y(i)*factors%val(iscale+i-1)
          enddo
        endif

      endif

      if (control%scaling.eq.1) deallocate (work,stat=info)
      if (info.eq.0) return

 100  sinfo%flag = -3
      sinfo%stat = info
      if (control%ldiag>0 .and. control%lp>0 ) &
          write (control%lp,'(/a/a,i5)') &
          'Error return from ma57_lmultiply', &
          'Allocate or deallocate failed with STAT=',info

      return

contains
      subroutine lmult(N,FACT,LFACT,IFACT,LIFACT,X,RESULT)
! This subroutine performs a multiplication by L using the factors
!     stored in FACT/IFACT by MA57BD.
      integer N,LFACT
      real (wp) FACT(LFACT)
      integer LIFACT,IFACT(LIFACT)
      real (wp) X(N),RESULT(N)
! N   must be set to the order of the matrix. It is not altered.
! FACT   must be set to hold the real values corresponding to the
!     factors. This must be unchanged since the preceding call to
!     MA57BD. It is not altered.
! LFACT  length of array FACT. It is not altered.
! IFACT  holds the integer indexing information for the matrix factors
!     in FACT. This must be unchanged since the preceding call to
!     MA57BD. It is not altered.
! LIFACT length of array IFACT. It is not altered.
! X on input, must be set to hold the right hand side vector.  On
!     return, it will hold the modified vector following forward
!     elimination.
! W   used as workspace to hold the components of the right hand
!     sides corresponding to current block pivotal rows.

! Procedures
      INTRINSIC ABS

! Constant
      real (wp) ONE,MONE,ZERO
      PARAMETER (ONE=1.0D0,ZERO=0.0D0,MONE=-1.0D0)
!
! Local variables
!

! APOS  Current position in array FACT.
! I     Temporary DO index
! IBLK  Index of block pivot.
! II    Temporary index.
! IWPOS Position in IFACT of start of current index list.
! NCOLS Number of columns in the block pivot.
! NROWS Number of rows in the block pivot.
      integer APOS,I,IBLK,II,IWPOS,NCOLS,NROWS

      real(wp), allocatable :: w(:)

      EXTERNAL DTPMV,DGEMV

! Could be allocated to maximum front size
      allocate (w(n))

      RESULT = 0

      APOS = 1
      IWPOS = 4

! Each pass through this loop performs the operations from one block
! of the factors.
      DO 270 IBLK = 1,IFACT(3)

! Find the number of rows and columns in the block.
        NCOLS = IFACT(IWPOS)
        NROWS = IFACT(IWPOS+1)
        IWPOS = IWPOS + 2

! Load appropriate components of right-hand sides into W.
        DO 10 I = 1,NROWS
          II = ABS(IFACT(IWPOS+I-1))
          W(I) = X(II)
   10   CONTINUE

! Treat diagonal block
        CALL DTPMV('L','N','U',NROWS,FACT(APOS),W,1)


! Update product vector
        DO I = 1,NROWS
          II = ABS(IFACT(IWPOS+I-1))
          RESULT(II) = RESULT(II) + W(I)
        ENDDO


        APOS = APOS + (NROWS* (NROWS+1))/2

        IF (NCOLS.GT.NROWS) THEN
! Treat off-diagonal block

! Reset W
          DO I = 1,NROWS
            II = ABS(IFACT(IWPOS+I-1))
            W(I) = X(II)
          ENDDO

          CALL DGEMV('N',NCOLS-NROWS,NROWS,MONE,FACT(APOS),NCOLS-NROWS,  &
                     W,1,ZERO,W(NROWS+1),1)

          APOS = APOS + NROWS* (NCOLS-NROWS)

! Update product vector
          DO I = NROWS+1,NCOLS
            II = ABS(IFACT(IWPOS+I-1))
            RESULT(II) = RESULT(II) + W(I)
          ENDDO

        ENDIF

        IWPOS = IWPOS + NCOLS

  270 CONTINUE

      deallocate (w)

      end subroutine lmult

      subroutine ltmult(N,FACT,LFACT,IFACT,LIFACT,X,RESULT)
! This subroutine performs multiplication by L transpose using the
!      factors stored in FACT/IFACT by MA57BD.
      integer N,LFACT
      real (wp) FACT(LFACT)
      integer LIFACT,IFACT(LIFACT)
      real (wp) X(N),RESULT(N)
! N   must be set to the order of the matrix. It is not altered.
! FACT   must be set to hold the real values corresponding to the
!     factors. This must be unchanged since the preceding call to
!     MA57BD. It is not altered.
! LFACT  length of array FACT. It is not altered.
! IFACT  holds the integer indexing information for the matrix factors
!     in FACT. This must be unchanged since the preceding call to
!     MA57BD. It is not altered.
! LIFACT length of array IFACT. It is not altered.
! X on input, must be set to hold the right hand side vector.  On
!     return, it will hold the modified vector following forward
!     elimination.
! W   used as workspace to hold the components of the right hand
!     sides corresponding to current block pivotal rows.
!

! Procedures
      INTRINSIC ABS

! Constant
      real (wp) ONE,MONE,ZERO
      PARAMETER (ONE=1.0D0,ZERO=0.0D0,MONE=-1.0D0)
!
! Local variables
! APOS  Current position in array FACT.
! I     Temporary DO index
! IBLK  Index of block pivot.
! II    Temporary index.
! IWPOS Position in IFACT of start of current index list.
! NCOLS Number of columns in the block pivot.
! NROWS Number of rows in the block pivot.
!
      integer APOS,I,IBLK,II,IWPOS,NCOLS,NROWS


      real(wp), allocatable :: w(:)

      EXTERNAL DTPMV,DGEMV

      allocate (w(n))

      RESULT = 0

      APOS = 1
      IWPOS = 4

! Each pass through this loop performs the operations from one block
! of the factors.
      DO 270 IBLK = 1,IFACT(3)

! Find the number of rows and columns in the block.
        NCOLS = IFACT(IWPOS)
        NROWS = IFACT(IWPOS+1)
        IWPOS = IWPOS + 2

! Load appropriate components of right-hand sides into W.
        DO 10 I = 1,NROWS
          II = ABS(IFACT(IWPOS+I-1))
          W(I) = X(II)
   10   CONTINUE

! Treat diagonal block
        CALL DTPMV('L','T','U',NROWS,FACT(APOS),W,1)

! Update product vector
        DO I = 1,NROWS
          II = ABS(IFACT(IWPOS+I-1))
          RESULT(II) = RESULT(II) + W(I)
        ENDDO

        APOS = APOS + (NROWS* (NROWS+1))/2

        IF (NCOLS.GT.NROWS) THEN
! Treat off-diagonal block

! Reset W
          DO I = 1,NCOLS-NROWS
            II = ABS(IFACT(IWPOS+NROWS+I-1))
            W(NROWS+I) = X(II)
          ENDDO

          CALL DGEMV('T',NCOLS-NROWS,NROWS,MONE,FACT(APOS),NCOLS-NROWS,  &
                     W(NROWS+1),1,ZERO,W,1)

          APOS = APOS + NROWS* (NCOLS-NROWS)

! Update product vector
          DO I = 1,NROWS
            II = ABS(IFACT(IWPOS+I-1))
            RESULT(II) = RESULT(II) + W(I)
          ENDDO

        ENDIF

        IWPOS = IWPOS + NCOLS

  270 CONTINUE

      deallocate (w)

      end subroutine ltmult

end subroutine ma57_lmultiply
      subroutine ma57_regen(nodes,nstk,tree,info)
! This subroutine regenerates child to parent pointers for the tree
! nodes is the number of nodes in the tree
! nstk is an array of length nodes giving number of assemblies at
!     each tree node
      integer nodes
! tree (output) gives child to parent pointers for the tree
! info is set to 1 if there is a failure in allocate or deallocate
      integer, intent(in)  :: nstk(nodes)
      integer, intent(out) :: tree(nodes),info

! Local variables
! i,j are do loop variables
! k is pointer to last/top entry in stack
! stack is work array holding stack

      INTEGER i,j,k,nass

      integer, allocatable :: stack(:)

      info = 0
      allocate(stack(nodes),stat=info)
      if (info.ne.0) return

! The tree pointer is zero at the roots
      tree = 0

      k = 0
      do i = 1,nodes
        nass = nstk(i)
! The previous nass nodes on the stack are now assembled
        if (nass .ne. 0) then
          do j = 1,nass
! Set child to parent pointers and remove children from stack
            tree(stack(k)) = i
            k = k - 1
          enddo
        endif
! Put node on stack
        k = k + 1
        stack(k) = i
      enddo

      deallocate(stack,stat=info)
      return
      end subroutine ma57_regen

      subroutine ma57_nflag(n,nodes,nzrhs,node,rhs,tree,flag)
! This subroutine sets flags for nodes that are required for forward
!     elimination for sparse right-hand sides with pattern given by rhs.
! n is the order of the matrix
! nodes is the number of nodes in the tree
! nzrhs is the number of entries in the right-hand side
! node is the number of the node at which variable is fully summed
! rhs is a list of the position of the entries in the sparse
!     right-hand side
! tree (input) gives child to parent pointers for the tree
! flag (output) is set to 1 only if node is used in solution of
!     sparse right-hand side
      integer n,nodes,nzrhs
      integer, intent(in)  :: node(n),rhs(nzrhs),tree(nodes)
      integer, intent(out) :: flag(nodes)

!local variables
! i,j are do loop variables
! inode holds the node number
! irhs holds a position in the sparse right-hand side
      INTEGER i,j,inode,irhs

      flag = 0

! Run through entries in sparse right-hand side
      do i = 1,nzrhs
        irhs = rhs(i)
! Identify node at which entry in rhs is fully-summed and flag it
        inode = node(irhs)
        if (flag(inode) == 1) exit
        flag(inode) = 1
! Set flag for all nodes from inode to root, stopping if existing
!     flag encountered
        do j = 1,nodes
          inode = tree(inode)
! Stop if at root
          if (inode == 0 ) exit
          if (flag(inode) == 1) exit
          flag(inode) = 1
        enddo
      enddo
      return
      end subroutine ma57_nflag

      subroutine dsolve(N,FACT,LFACT,IFACT,LIFACT,RHS,FLAG,blascntl,info)
! This subroutine divides a vector by the block diagonal matrix of
!     the matrix factors using factor entries stored in FACT/IFACT
!      by MA57BD.
      INTEGER N,LFACT
      real (wp) FACT(LFACT)
      INTEGER LIFACT,IFACT(LIFACT)
      real (wp) RHS(N)
      INTEGER blascntl,info
      INTEGER FLAG(*)
! FACT    must be set to hold the real values corresponding to the
!      factors. This must be unchanged since the
!      preceding call to MA57BD. It is not altered.
! LFACT   length of array FACT. It is not altered.
! IFACT   holds the integer indexing information for the matrix factors
!      in FACT. This must be unchanged since the preceding call to
!      MA57BD. It is not altered.
! LIFACT  length of array IFACT. It is not altered.
! RHS  on entry, must be set to hold the right hand side modified by
!      the forward substitution operations. On exit, holds the
!      solution vector.
! W    used as workspace to hold the components of the right hand
!      sides corresponding to current block pivotal rows.
! blascntl Threshold on number of columns in a block for direct
!       addressing using Level 2 and Level 3 BLAS.

      real (wp), allocatable :: w(:)
! Procedures
      INTRINSIC ABS
      EXTERNAL DGEMM,DTPSV

!
! Local variables.
      INTEGER APOS,APOS2,I,IBLK,II,IPIV,IRHS,IWPOS,JPIV,NCOLS,NROWS
! APOS  Current position in array FACT.
! APOS2 Current position in array FACT for off-diagonal entry of 2x2
!       pivot.
! I     Temporary DO index
! IBLK  Index of block pivot.
! II    Temporary index.
! IPIV  Pivot index.
! IRHS  RHS index.
! IWPOS Position in IFACT of start of current index list.
! JPIV  Has the value 1 for the first row of a 2 by 2 pivot and -1 for
!       the second.
! K     Temporary pointer to position in real array.
! NCOLS Number of columns in the block pivot.
! NROWS Number of rows in the block pivot.
!
      info = 0

      allocate (w(n),stat=info)
      if (info.ne.0) return

      APOS = 1
      APOS2 = IFACT(1)
      IWPOS = 4

! Run through block pivot rows
      DO 380 IBLK = 1,IFACT(3)

! Find the number of rows and columns in the block.
        NCOLS = IFACT(IWPOS)
        NROWS = IFACT(IWPOS+1)
        IWPOS = IWPOS + 2

        IF (FLAG(IBLK) == 0) THEN
! We can skip this block just resetting IWPOS and APOS
          IWPOS = IWPOS + NCOLS
          APOS  = APOS + (NROWS*(NROWS+1))/2
          APOS  = APOS + NROWS*(NCOLS-NROWS)
          GO TO 380
        ENDIF

        IF (NROWS.GT.4 .AND. NCOLS.GT.blascntl) THEN

! Perform operations using direct addressing.

! Multiply by the diagonal matrix (direct addressing)
          DO 10 IPIV = 1,NROWS
            IRHS = ABS(IFACT(IWPOS+IPIV-1))
            W(IPIV) = RHS(IRHS)*SQRT(FACT(APOS))
            APOS = APOS + (NROWS+1-IPIV)
   10     CONTINUE

! Reload W back into RHS.
          DO 60 I = 1,NROWS
            II = ABS(IFACT(IWPOS+I-1))
            RHS(II) = W(I)
   60     CONTINUE

        ELSE
!
! Perform operations using indirect addressing.

! Multiply by the diagonal matrix (indirect addressing)
          JPIV = 1
          DO 210 IPIV = 1,NROWS
            IRHS = IFACT(IWPOS+IPIV-1)

! 1 by 1 pivot.
            RHS(IRHS) = RHS(IRHS)*SQRT(FACT(APOS))

            APOS = APOS + NROWS - IPIV + 1

  210     CONTINUE

        END IF

        IWPOS = IWPOS + NCOLS
        APOS = APOS + NROWS*(NCOLS-NROWS)

  380 CONTINUE
        deallocate(w,stat=info)
        return

      end subroutine dsolve

end module hsl_ma57_double
