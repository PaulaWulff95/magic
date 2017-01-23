module finite_differences
   !
   ! This module is used to calculate the radial grid when finite differences
   ! are used
   !

   use precision_mod
   use parallel_mod, only: rank
   use constants, only: one, two
   use useful, only: logWrite
   use mem_alloc, only: bytes_allocated
   use radial_scheme, only: type_rscheme

   implicit none

   private

   type, public, extends(type_rscheme) :: type_fd
   contains
      procedure :: initialize
      procedure :: finalize
      procedure, private :: nullify_epsilon
      procedure, private :: get_FD_coeffs
      procedure :: get_der_mat
      procedure :: get_grid => get_FD_grid
   end type type_fd

contains

   subroutine initialize(this,n_r_max,order,order_boundary)
      !
      ! This subroutine allocates the arrays used when finite difference are used
      !

      class(type_fd) :: this
      integer, intent(in) :: n_r_max ! Number of radial grid points
      integer, intent(in) :: order   ! FD order
      integer, intent(in) :: order_boundary   ! FD order on the boundary

      this%order = order
      this%order_boundary = order_boundary
      this%rnorm = one
      this%n_max = n_r_max
      this%boundary_fac = one
      this%version = 'fd'

      allocate( this%dr(this%n_max,0:order) )
      allocate( this%ddr(this%n_max,0:order) )
      allocate( this%dddr(this%n_max,0:order+2) )
      allocate( this%dr_top(order/2,0:order_boundary) )
      allocate( this%dr_bot(order/2,0:order_boundary) )
      allocate( this%ddr_top(order/2,0:order_boundary+1) )
      allocate( this%ddr_bot(order/2,0:order_boundary+1) )
      allocate( this%dddr_top(order/2+1,0:order_boundary+2) )
      allocate( this%dddr_bot(order/2+1,0:order_boundary+2) )

      bytes_allocated=bytes_allocated+(this%n_max*(3*order+5)+         &
      &               order/2*(4*order_boundary+6)+                    &
      &               (order/2+1)*(2*order_boundary+6))*SIZEOF_DEF_REAL

      allocate( this%rMat(this%n_max,this%n_max) )
      allocate( this%drMat(this%n_max,this%n_max) )
      allocate( this%d2rMat(this%n_max,this%n_max) )
      allocate( this%d3rMat(this%n_max,this%n_max) )

      bytes_allocated=bytes_allocated+4*this%n_max*this%n_max*SIZEOF_DEF_REAL

   end subroutine initialize
!---------------------------------------------------------------------------
   subroutine finalize(this)
      !
      ! This subroutine deallocates the arrays used in FD
      !

      class(type_fd) :: this

      deallocate( this%dr, this%ddr, this%dddr )
      deallocate( this%dr_top, this%dr_bot )
      deallocate( this%ddr_top, this%ddr_bot )
      deallocate( this%dddr_top, this%dddr_bot )
      deallocate( this%rMat, this%drMat, this%d2rMat, this%d3rMat )

   end subroutine finalize
!---------------------------------------------------------------------------
   subroutine nullify_epsilon(this)
      !
      ! Nullify epsilon values
      !

      class(type_fd) :: this

      !-- Local variables
      real(cp) :: eps
      integer :: n_r, od

      eps = 100.0_cp*epsilon(one)

      !-- Bulk points for 1st and 2nd derivative
      do od=0,this%order
         do n_r=1,this%n_max
            if ( abs(this%dr(n_r,od)) < eps ) this%dr(n_r,od)=0.0_cp
            if ( abs(this%ddr(n_r,od)) < eps ) this%ddr(n_r,od)=0.0_cp
         end do
      end do

      !-- Bulk points for 3rd derivative
      do od=0,this%order+2
         do n_r=1,this%n_max
            if ( abs(this%dddr(n_r,od)) < eps ) this%dddr(n_r,od)=0.0_cp
         end do
      end do

      !-- Boundary points for 1st derivative
      do od=0,this%order_boundary
         do n_r=1,this%order/2
            if ( abs(this%dr_bot(n_r,od)) < eps ) this%dr_bot(n_r,od)=0.0_cp
            if ( abs(this%dr_top(n_r,od)) < eps ) this%dr_top(n_r,od)=0.0_cp
         end do
      end do

      !-- Boundary points for 2nd derivative
      do od=0,this%order_boundary+1
         do n_r=1,this%order/2
            if ( abs(this%ddr_bot(n_r,od)) < eps ) this%ddr_bot(n_r,od)=0.0_cp
            if ( abs(this%ddr_top(n_r,od)) < eps ) this%ddr_top(n_r,od)=0.0_cp
         end do
      end do

      !-- Boundary points for 3rd derivative
      do od=0,this%order_boundary+2
         do n_r=1,this%order/2+1
            if ( abs(this%dddr_bot(n_r,od)) < eps ) this%dddr_bot(n_r,od)=0.0_cp
            if ( abs(this%dddr_top(n_r,od)) < eps ) this%dddr_top(n_r,od)=0.0_cp
         end do
      end do

   end subroutine nullify_epsilon
!---------------------------------------------------------------------------
   subroutine get_FD_grid(this, n_r_max, ricb, rcmb, ratio1, ratio2, r)
      !
      ! This subroutine constructs the radial grid
      !

      class(type_fd) :: this

      !-- Input quantities:
      integer,  intent(in) :: n_r_max    ! Number of grid points
      real(cp), intent(inout) :: ratio1  ! Nboudary/Nbulk
      real(cp), intent(in) :: ratio2     ! drMin/drMax
      real(cp), intent(in) :: ricb       ! inner boundary
      real(cp), intent(in) :: rcmb       ! outer boundary

      !-- Output quantities:
      real(cp), intent(out) :: r(n_r_max) ! radius

      !-- Local quantities:
      real(cp) :: dr_before, dr_after
      character(len=80) :: message
      integer :: n_r
      integer :: n_bulk_points, n_boundary_points

      r(1)=rcmb ! start with the outer boundary

      if ( ratio2 == one ) then ! Regular grid

         dr_before = one/(real(n_r_max,cp)-one)
         dr_after = dr_before
         do n_r=2,n_r_max
            r(n_r)=r(n_r-1)-dr_before
         end do

      else  ! irregular grid

         n_boundary_points = int( real(n_r_max-1,cp)/(two*(one+ratio1)) )
         ratio1 = real(n_r_max-1, cp)/real(two*n_boundary_points)-one

         n_bulk_points = n_r_max-1-2*n_boundary_points
         dr_after  = exp( log(ratio2)/real(n_boundary_points,cp) )
         dr_before = one

         do n_r=1,n_boundary_points
            dr_before=dr_before*dr_after
         end do
         dr_before=one/(real(n_bulk_points,cp)+two*dr_after*((one-dr_before) &
         &         /(one-dr_after)))

         write(message,'(''!      drMax='',ES16.6)') dr_before
         call logWrite(message)

         do n_r=1,n_boundary_points
            dr_before = dr_before*dr_after
         end do

         write(message,'(''!      drMax='',ES16.6)') dr_before
         call logWrite(message)

         do n_r=2,n_boundary_points+1
            r(n_r) = r(n_r-1)-dr_before
            dr_before = dr_before/dr_after
         end do

         do n_r=1,n_bulk_points
            r(n_r+n_boundary_points+1)=r(n_r+n_boundary_points)-dr_before
         end do

         do n_r=1,n_boundary_points
            dr_before = dr_before*dr_after
            r(n_r+n_boundary_points+n_bulk_points+1)=         &
            &        r(n_r+n_boundary_points+n_bulk_points)-dr_before
         end do

      end if

      if ( abs(r(n_r_max)-ricb) > 10.0_cp*epsilon(one) ) then
         if ( rank == 0 ) then
            write(*,*) 'Wrong internal radius!'
            stop
         end if
      end if

      call this%get_FD_coeffs(r)

   end subroutine get_FD_grid
!---------------------------------------------------------------------------
   subroutine get_FD_coeffs(this, r)

      class(type_fd) :: this

      !-- Input quantities:
      real(cp), intent(in) :: r(:) ! Radius

      !-- Local quantities:
      real(cp), allocatable :: dr_spacing(:)
      real(cp), allocatable :: taylor_exp(:,:)
      integer :: n_r, od

      !
      !-- Step 1: First and 2nd derivatives in the bulk
      !
      allocate( dr_spacing(this%order+1) )
      allocate( taylor_exp(0:this%order,0:this%order) )

      do n_r=1+this%order/2,this%n_max-this%order/2
         do od=0,this%order
            dr_spacing(od+1)=r(n_r-this%order/2+od)-r(n_r)
         end do

         call populate_fd_weights(0.0_cp,dr_spacing,this%order, &
              &                   this%order,taylor_exp)

         do od=0,this%order
            this%dr(n_r,od) =taylor_exp(od,1)
            this%ddr(n_r,od)=taylor_exp(od,2)
         end do
      end do

      deallocate( dr_spacing, taylor_exp )

      !
      !-- Step 2: First derivative for the outer points
      !
      allocate( dr_spacing(this%order_boundary+1) )
      allocate( taylor_exp(0:this%order_boundary,0:this%order_boundary) )

      do n_r=1,this%order/2
         do od=0,this%order_boundary
            dr_spacing(od+1)=r(od+1)-r(n_r)
         end do

         call populate_fd_weights(0.0_cp,dr_spacing,this%order_boundary, &
              &                   this%order_boundary,taylor_exp)

         do od=0,this%order_boundary
            this%dr_top(n_r,od) =taylor_exp(od,1)
         end do
      end do

      !
      !-- Step 3: First derivative for the inner points
      !
      do n_r=1,this%order/2
         do od=0,this%order_boundary
            dr_spacing(od+1)=r(this%n_max-od)-r(this%n_max-n_r+1)
         end do

         call populate_fd_weights(0.0_cp,dr_spacing,this%order_boundary, &
              &                   this%order_boundary,taylor_exp)

         do od=0,this%order_boundary
            this%dr_bot(n_r,od) =taylor_exp(od,1)
         end do
      end do

      deallocate( dr_spacing, taylor_exp )

      !
      !-- Step 4: 2nd derivative for the outer points
      !
      allocate( dr_spacing(this%order_boundary+2) )
      allocate( taylor_exp(0:this%order_boundary+1,0:this%order_boundary+1) )

      do n_r=1,this%order/2
         do od=0,this%order_boundary+1
            dr_spacing(od+1)=r(od+1)-r(n_r)
         end do

         call populate_fd_weights(0.0_cp,dr_spacing,this%order_boundary+1, &
              &                   this%order_boundary+1,taylor_exp)

         do od=0,this%order_boundary+1
            this%ddr_top(n_r,od) =taylor_exp(od,2)
         end do
      end do

      !
      !-- Step 5: 2nd derivative for the inner points
      !
      do n_r=1,this%order/2
         do od=0,this%order_boundary+1
            dr_spacing(od+1)=r(this%n_max-od)-r(this%n_max-n_r+1)
         end do

         call populate_fd_weights(0.0_cp,dr_spacing,this%order_boundary+1, &
              &                   this%order_boundary+1,taylor_exp)

         do od=0,this%order_boundary+1
            this%ddr_bot(n_r,od) =taylor_exp(od,2)
         end do
      end do

      deallocate( dr_spacing, taylor_exp )

      !
      !-- Step 6: 3rd derivative in the bulk
      !
      allocate( dr_spacing(this%order+3) )
      allocate( taylor_exp(0:this%order+2,0:this%order+2) )

      do n_r=2+this%order/2,this%n_max-this%order/2-1
         do od=0,this%order+2
            dr_spacing(od+1)=r(n_r-this%order/2-1+od)-r(n_r)
         end do

         call populate_fd_weights(0.0_cp,dr_spacing,this%order+2,this%order+2, &
              &                   taylor_exp)
         do od=0,this%order+2
            this%dddr(n_r,od) =taylor_exp(od,3)
         end do
      end do

      deallocate( dr_spacing, taylor_exp )

      !
      !-- Step 7: 3rd derivative for the outer points
      !
      allocate( dr_spacing(this%order_boundary+3) )
      allocate( taylor_exp(0:this%order_boundary+2,0:this%order_boundary+2) )

      do n_r=1,this%order/2+1
         do od=0,this%order_boundary+2
            dr_spacing(od+1)=r(od+1)-r(n_r)
         end do

         call populate_fd_weights(0.0_cp,dr_spacing,this%order_boundary+2, &
              &                   this%order_boundary+2,taylor_exp)

         do od=0,this%order_boundary+2
            this%dddr_top(n_r,od) =taylor_exp(od,3)
         end do
      end do

      !
      !-- Step 8: 3rd derivative for the inner points
      !
      do n_r=1,this%order/2+1
         do od=0,this%order_boundary+2
            dr_spacing(od+1)=r(this%n_max-od)-r(this%n_max-n_r+1)
         end do

         call populate_fd_weights(0.0_cp,dr_spacing,this%order_boundary+2, &
              &                   this%order_boundary+2,taylor_exp)

         do od=0,this%order_boundary+2
            this%dddr_bot(n_r,od) =taylor_exp(od,3)
         end do
      end do

      deallocate( dr_spacing, taylor_exp )

   end subroutine get_FD_coeffs
!---------------------------------------------------------------------------
   subroutine get_der_mat(this, n_r_max)

      class(type_fd) :: this

      !-- Input variable
      integer,  intent(in) :: n_r_max

      !-- Local variables
      integer :: i, j, n_r

      !-- Set everything to zero
      do j=1,n_r_max
         do i=1,n_r_max
            this%drMat(i,j) =0.0_cp
            this%d2rMat(i,j)=0.0_cp
            this%d3rMat(i,j)=0.0_cp
            this%rMat(i,j)  =0.0_cp
         end do
         this%rMat(j,j)=1.0_cp
      end do

      !-- Bulk points for 1st and 2nd derivatives
      do n_r=this%order/2+1,n_r_max-this%order/2
         this%drMat(n_r,n_r-this%order/2:n_r+this%order/2) = this%dr(n_r,:)
         this%d2rMat(n_r,n_r-this%order/2:n_r+this%order/2)=this%ddr(n_r,:)
      end do

      !-- Bulk points for 3rd derivative
      do n_r=2+this%order/2,n_r_max-this%order/2-1
         this%d3rMat(n_r,n_r-this%order/2-1:n_r+this%order/2+1)=this%dddr(n_r,:)
      end do

      !-- Boundary points for 1st derivative
      do n_r=1,this%order/2
         this%drMat(n_r,1:this%order_boundary+1)                         =this%dr_top(n_r,:)
         this%drMat(n_r_max-n_r+1,n_r_max:n_r_max-this%order_boundary:-1)=this%dr_bot(n_r,:)
      end do

      !-- Boundary points for 2nd derivative
      do n_r=1,this%order/2
         this%d2rMat(n_r,1:this%order_boundary+2)                           =this%ddr_top(n_r,:)
         this%d2rMat(n_r_max-n_r+1,n_r_max:n_r_max-this%order_boundary-1:-1)=this%ddr_bot(n_r,:)
      end do

      !-- Boundary points for 3rd derivative
      do n_r=1,this%order/2+1
         this%d3rMat(n_r,1:this%order_boundary+3)                           =this%dddr_top(n_r,:)
         this%d3rMat(n_r_max-n_r+1,n_r_max:n_r_max-this%order_boundary-2:-1)=this%dddr_bot(n_r,:)
      end do

   end subroutine get_der_mat
!----------------------------------------------------------------------------
   subroutine populate_fd_weights(z, x, nd, m, c)
      !
      !  Input Parameters
      !    z            -  location where approximations are to be
      !                    accurate
      !    x(0:nd)      -  grid point locations, found in x(0:n)
      !    nd           -  dimension of x- and c-arrays in calling
      !                    program x(0:nd) and c(0:nd, 0:m), respectively
      !    m            -  highest derivative for which weights are
      !                    sought
      !
      !  Output Parameter
      !    c(0:nd,0:m)  -  weights at grid locations x(0:n) for
      !                    derivatives of order 0:m, found in c(0:nd, 0:m)
      !
      !  Reference:
      !      Generation of Finite Difference Formulas on Arbitrarily
      !          Spaced Grids, Bengt Fornberg,
      !          Mathematics of compuation, 51, 184, 1988, 699-706

      !-- Input quantities:
      real(cp), intent(in) :: z
      integer,  intent(in) :: nd, m
      real(cp), intent(in) :: x(0:nd)

      !-- Output:
      real(cp), intent(out) :: c(0:nd, 0:m)

      !-- Local variables
      real(cp) :: c1, c2, c3, c4, c5
      integer :: i, j, k, mn

      c1 = 1.0_cp
      c4 = x(0) - z
      c(:,:)  = 0.0_cp
      c(0, 0) = 1.0_cp
      do i=1, nd
         mn = min(i, m)
         c2 = 1.0_cp
         c5 = c4
         c4 = x(i) - z
         do j=0, i-1
            c3 = x(i) - x(j)
            c2 = c2*c3
            if (j == i-1) then
               do k = mn, 1, -1
                   c(i, k) = c1*(k*c(i-1, k-1) - c5*c(i-1, k))/c2
               end do
               c(i, 0) = -c1*c5*c(i-1, 0)/c2
            endif
            do k=mn, 1, -1
               c(j, k) = (c4*c(j, k) - k*c(j, k-1))/c3
            end do
            c(j, 0) = c4*c(j, 0)/c3
         end do
         c1 = c2
      end do

   end subroutine populate_fd_weights
!------------------------------------------------------------------------------
end module finite_differences
