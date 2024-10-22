! Generated automatically.  DO NOT EDIT!

  include 'fftw3.f03'

  integer(C_INTPTR_T), parameter :: FFTW_MPI_DEFAULT_BLOCK = 0
  integer(C_INT), parameter :: FFTW_MPI_SCRAMBLED_IN = 134217728
  integer(C_INT), parameter :: FFTW_MPI_SCRAMBLED_OUT = 268435456
  integer(C_INT), parameter :: FFTW_MPI_TRANSPOSED_IN = 536870912
  integer(C_INT), parameter :: FFTW_MPI_TRANSPOSED_OUT = 1073741824

  type, bind(C) :: fftw_mpi_ddim
     integer(C_INTPTR_T) n, ib, ob
  end type fftw_mpi_ddim

  interface
    subroutine fftw_mpi_init() bind(C, name='fftw_mpi_init')
      import
    end subroutine fftw_mpi_init
    
    subroutine fftw_mpi_cleanup() bind(C, name='fftw_mpi_cleanup')
      import
    end subroutine fftw_mpi_cleanup
    
    integer(C_INTPTR_T) function fftw_mpi_local_size_many_transposed(rnk,n,howmany,block0,block1,comm,local_n0,local_0_start, &
                                                                     local_n1,local_1_start) &
                                 bind(C, name='fftw_mpi_local_size_many_transposed_f03')
      import
      integer(C_INT), value :: rnk
      integer(C_INTPTR_T), dimension(*), intent(in) :: n
      integer(C_INTPTR_T), value :: howmany
      integer(C_INTPTR_T), value :: block0
      integer(C_INTPTR_T), value :: block1
      integer(C_INT32_T), value :: comm
      integer(C_INTPTR_T), intent(out) :: local_n0
      integer(C_INTPTR_T), intent(out) :: local_0_start
      integer(C_INTPTR_T), intent(out) :: local_n1
      integer(C_INTPTR_T), intent(out) :: local_1_start
    end function fftw_mpi_local_size_many_transposed
    
    integer(C_INTPTR_T) function fftw_mpi_local_size_many(rnk,n,howmany,block0,comm,local_n0,local_0_start) &
                                 bind(C, name='fftw_mpi_local_size_many_f03')
      import
      integer(C_INT), value :: rnk
      integer(C_INTPTR_T), dimension(*), intent(in) :: n
      integer(C_INTPTR_T), value :: howmany
      integer(C_INTPTR_T), value :: block0
      integer(C_INT32_T), value :: comm
      integer(C_INTPTR_T), intent(out) :: local_n0
      integer(C_INTPTR_T), intent(out) :: local_0_start
    end function fftw_mpi_local_size_many
    
    integer(C_INTPTR_T) function fftw_mpi_local_size_transposed(rnk,n,comm,local_n0,local_0_start,local_n1,local_1_start) &
                                 bind(C, name='fftw_mpi_local_size_transposed_f03')
      import
      integer(C_INT), value :: rnk
      integer(C_INTPTR_T), dimension(*), intent(in) :: n
      integer(C_INT32_T), value :: comm
      integer(C_INTPTR_T), intent(out) :: local_n0
      integer(C_INTPTR_T), intent(out) :: local_0_start
      integer(C_INTPTR_T), intent(out) :: local_n1
      integer(C_INTPTR_T), intent(out) :: local_1_start
    end function fftw_mpi_local_size_transposed
    
    integer(C_INTPTR_T) function fftw_mpi_local_size(rnk,n,comm,local_n0,local_0_start) bind(C, name='fftw_mpi_local_size_f03')
      import
      integer(C_INT), value :: rnk
      integer(C_INTPTR_T), dimension(*), intent(in) :: n
      integer(C_INT32_T), value :: comm
      integer(C_INTPTR_T), intent(out) :: local_n0
      integer(C_INTPTR_T), intent(out) :: local_0_start
    end function fftw_mpi_local_size
    
    integer(C_INTPTR_T) function fftw_mpi_local_size_many_1d(n0,howmany,comm,sign,flags,local_ni,local_i_start,local_no, &
                                                             local_o_start) bind(C, name='fftw_mpi_local_size_many_1d_f03')
      import
      integer(C_INTPTR_T), value :: n0
      integer(C_INTPTR_T), value :: howmany
      integer(C_INT32_T), value :: comm
      integer(C_INT), value :: sign
      integer(C_INT), value :: flags
      integer(C_INTPTR_T), intent(out) :: local_ni
      integer(C_INTPTR_T), intent(out) :: local_i_start
      integer(C_INTPTR_T), intent(out) :: local_no
      integer(C_INTPTR_T), intent(out) :: local_o_start
    end function fftw_mpi_local_size_many_1d
    
    integer(C_INTPTR_T) function fftw_mpi_local_size_1d(n0,comm,sign,flags,local_ni,local_i_start,local_no,local_o_start) &
                                 bind(C, name='fftw_mpi_local_size_1d_f03')
      import
      integer(C_INTPTR_T), value :: n0
      integer(C_INT32_T), value :: comm
      integer(C_INT), value :: sign
      integer(C_INT), value :: flags
      integer(C_INTPTR_T), intent(out) :: local_ni
      integer(C_INTPTR_T), intent(out) :: local_i_start
      integer(C_INTPTR_T), intent(out) :: local_no
      integer(C_INTPTR_T), intent(out) :: local_o_start
    end function fftw_mpi_local_size_1d
    
    integer(C_INTPTR_T) function fftw_mpi_local_size_2d(n0,n1,comm,local_n0,local_0_start) &
                                 bind(C, name='fftw_mpi_local_size_2d_f03')
      import
      integer(C_INTPTR_T), value :: n0
      integer(C_INTPTR_T), value :: n1
      integer(C_INT32_T), value :: comm
      integer(C_INTPTR_T), intent(out) :: local_n0
      integer(C_INTPTR_T), intent(out) :: local_0_start
    end function fftw_mpi_local_size_2d
    
    integer(C_INTPTR_T) function fftw_mpi_local_size_2d_transposed(n0,n1,comm,local_n0,local_0_start,local_n1,local_1_start) &
                                 bind(C, name='fftw_mpi_local_size_2d_transposed_f03')
      import
      integer(C_INTPTR_T), value :: n0
      integer(C_INTPTR_T), value :: n1
      integer(C_INT32_T), value :: comm
      integer(C_INTPTR_T), intent(out) :: local_n0
      integer(C_INTPTR_T), intent(out) :: local_0_start
      integer(C_INTPTR_T), intent(out) :: local_n1
      integer(C_INTPTR_T), intent(out) :: local_1_start
    end function fftw_mpi_local_size_2d_transposed
    
    integer(C_INTPTR_T) function fftw_mpi_local_size_3d(n0,n1,n2,comm,local_n0,local_0_start) &
                                 bind(C, name='fftw_mpi_local_size_3d_f03')
      import
      integer(C_INTPTR_T), value :: n0
      integer(C_INTPTR_T), value :: n1
      integer(C_INTPTR_T), value :: n2
      integer(C_INT32_T), value :: comm
      integer(C_INTPTR_T), intent(out) :: local_n0
      integer(C_INTPTR_T), intent(out) :: local_0_start
    end function fftw_mpi_local_size_3d
    
    integer(C_INTPTR_T) function fftw_mpi_local_size_3d_transposed(n0,n1,n2,comm,local_n0,local_0_start,local_n1,local_1_start) &
                                 bind(C, name='fftw_mpi_local_size_3d_transposed_f03')
      import
      integer(C_INTPTR_T), value :: n0
      integer(C_INTPTR_T), value :: n1
      integer(C_INTPTR_T), value :: n2
      integer(C_INT32_T), value :: comm
      integer(C_INTPTR_T), intent(out) :: local_n0
      integer(C_INTPTR_T), intent(out) :: local_0_start
      integer(C_INTPTR_T), intent(out) :: local_n1
      integer(C_INTPTR_T), intent(out) :: local_1_start
    end function fftw_mpi_local_size_3d_transposed
    
    type(C_PTR) function fftw_mpi_plan_many_transpose(n0,n1,howmany,block0,block1,in,out,comm,flags) &
                         bind(C, name='fftw_mpi_plan_many_transpose_f03')
      import
      integer(C_INTPTR_T), value :: n0
      integer(C_INTPTR_T), value :: n1
      integer(C_INTPTR_T), value :: howmany
      integer(C_INTPTR_T), value :: block0
      integer(C_INTPTR_T), value :: block1
      real(C_DOUBLE), dimension(*), intent(out) :: in
      real(C_DOUBLE), dimension(*), intent(out) :: out
      integer(C_INT32_T), value :: comm
      integer(C_INT), value :: flags
    end function fftw_mpi_plan_many_transpose
    
    type(C_PTR) function fftw_mpi_plan_transpose(n0,n1,in,out,comm,flags) bind(C, name='fftw_mpi_plan_transpose_f03')
      import
      integer(C_INTPTR_T), value :: n0
      integer(C_INTPTR_T), value :: n1
      real(C_DOUBLE), dimension(*), intent(out) :: in
      real(C_DOUBLE), dimension(*), intent(out) :: out
      integer(C_INT32_T), value :: comm
      integer(C_INT), value :: flags
    end function fftw_mpi_plan_transpose
    
    type(C_PTR) function fftw_mpi_plan_many_dft(rnk,n,howmany,block,tblock,in,out,comm,sign,flags) &
                         bind(C, name='fftw_mpi_plan_many_dft_f03')
      import
      integer(C_INT), value :: rnk
      integer(C_INTPTR_T), dimension(*), intent(in) :: n
      integer(C_INTPTR_T), value :: howmany
      integer(C_INTPTR_T), value :: block
      integer(C_INTPTR_T), value :: tblock
      complex(C_DOUBLE_COMPLEX), dimension(*), intent(out) :: in
      complex(C_DOUBLE_COMPLEX), dimension(*), intent(out) :: out
      integer(C_INT32_T), value :: comm
      integer(C_INT), value :: sign
      integer(C_INT), value :: flags
    end function fftw_mpi_plan_many_dft
    
    type(C_PTR) function fftw_mpi_plan_dft(rnk,n,in,out,comm,sign,flags) bind(C, name='fftw_mpi_plan_dft_f03')
      import
      integer(C_INT), value :: rnk
      integer(C_INTPTR_T), dimension(*), intent(in) :: n
      complex(C_DOUBLE_COMPLEX), dimension(*), intent(out) :: in
      complex(C_DOUBLE_COMPLEX), dimension(*), intent(out) :: out
      integer(C_INT32_T), value :: comm
      integer(C_INT), value :: sign
      integer(C_INT), value :: flags
    end function fftw_mpi_plan_dft
    
    type(C_PTR) function fftw_mpi_plan_dft_1d(n0,in,out,comm,sign,flags) bind(C, name='fftw_mpi_plan_dft_1d_f03')
      import
      integer(C_INTPTR_T), value :: n0
      complex(C_DOUBLE_COMPLEX), dimension(*), intent(out) :: in
      complex(C_DOUBLE_COMPLEX), dimension(*), intent(out) :: out
      integer(C_INT32_T), value :: comm
      integer(C_INT), value :: sign
      integer(C_INT), value :: flags
    end function fftw_mpi_plan_dft_1d
    
    type(C_PTR) function fftw_mpi_plan_dft_2d(n0,n1,in,out,comm,sign,flags) bind(C, name='fftw_mpi_plan_dft_2d_f03')
      import
      integer(C_INTPTR_T), value :: n0
      integer(C_INTPTR_T), value :: n1
      complex(C_DOUBLE_COMPLEX), dimension(*), intent(out) :: in
      complex(C_DOUBLE_COMPLEX), dimension(*), intent(out) :: out
      integer(C_INT32_T), value :: comm
      integer(C_INT), value :: sign
      integer(C_INT), value :: flags
    end function fftw_mpi_plan_dft_2d
    
    type(C_PTR) function fftw_mpi_plan_dft_3d(n0,n1,n2,in,out,comm,sign,flags) bind(C, name='fftw_mpi_plan_dft_3d_f03')
      import
      integer(C_INTPTR_T), value :: n0
      integer(C_INTPTR_T), value :: n1
      integer(C_INTPTR_T), value :: n2
      complex(C_DOUBLE_COMPLEX), dimension(*), intent(out) :: in
      complex(C_DOUBLE_COMPLEX), dimension(*), intent(out) :: out
      integer(C_INT32_T), value :: comm
      integer(C_INT), value :: sign
      integer(C_INT), value :: flags
    end function fftw_mpi_plan_dft_3d
    
    type(C_PTR) function fftw_mpi_plan_many_r2r(rnk,n,howmany,iblock,oblock,in,out,comm,kind,flags) &
                         bind(C, name='fftw_mpi_plan_many_r2r_f03')
      import
      integer(C_INT), value :: rnk
      integer(C_INTPTR_T), dimension(*), intent(in) :: n
      integer(C_INTPTR_T), value :: howmany
      integer(C_INTPTR_T), value :: iblock
      integer(C_INTPTR_T), value :: oblock
      real(C_DOUBLE), dimension(*), intent(out) :: in
      real(C_DOUBLE), dimension(*), intent(out) :: out
      integer(C_INT32_T), value :: comm
      integer(C_FFTW_R2R_KIND), dimension(*), intent(in) :: kind
      integer(C_INT), value :: flags
    end function fftw_mpi_plan_many_r2r
    
    type(C_PTR) function fftw_mpi_plan_r2r(rnk,n,in,out,comm,kind,flags) bind(C, name='fftw_mpi_plan_r2r_f03')
      import
      integer(C_INT), value :: rnk
      integer(C_INTPTR_T), dimension(*), intent(in) :: n
      real(C_DOUBLE), dimension(*), intent(out) :: in
      real(C_DOUBLE), dimension(*), intent(out) :: out
      integer(C_INT32_T), value :: comm
      integer(C_FFTW_R2R_KIND), dimension(*), intent(in) :: kind
      integer(C_INT), value :: flags
    end function fftw_mpi_plan_r2r
    
    type(C_PTR) function fftw_mpi_plan_r2r_2d(n0,n1,in,out,comm,kind0,kind1,flags) bind(C, name='fftw_mpi_plan_r2r_2d_f03')
      import
      integer(C_INTPTR_T), value :: n0
      integer(C_INTPTR_T), value :: n1
      real(C_DOUBLE), dimension(*), intent(out) :: in
      real(C_DOUBLE), dimension(*), intent(out) :: out
      integer(C_INT32_T), value :: comm
      integer(C_FFTW_R2R_KIND), value :: kind0
      integer(C_FFTW_R2R_KIND), value :: kind1
      integer(C_INT), value :: flags
    end function fftw_mpi_plan_r2r_2d
    
    type(C_PTR) function fftw_mpi_plan_r2r_3d(n0,n1,n2,in,out,comm,kind0,kind1,kind2,flags) bind(C, name='fftw_mpi_plan_r2r_3d_f03')
      import
      integer(C_INTPTR_T), value :: n0
      integer(C_INTPTR_T), value :: n1
      integer(C_INTPTR_T), value :: n2
      real(C_DOUBLE), dimension(*), intent(out) :: in
      real(C_DOUBLE), dimension(*), intent(out) :: out
      integer(C_INT32_T), value :: comm
      integer(C_FFTW_R2R_KIND), value :: kind0
      integer(C_FFTW_R2R_KIND), value :: kind1
      integer(C_FFTW_R2R_KIND), value :: kind2
      integer(C_INT), value :: flags
    end function fftw_mpi_plan_r2r_3d
    
    type(C_PTR) function fftw_mpi_plan_many_dft_r2c(rnk,n,howmany,iblock,oblock,in,out,comm,flags) &
                         bind(C, name='fftw_mpi_plan_many_dft_r2c_f03')
      import
      integer(C_INT), value :: rnk
      integer(C_INTPTR_T), dimension(*), intent(in) :: n
      integer(C_INTPTR_T), value :: howmany
      integer(C_INTPTR_T), value :: iblock
      integer(C_INTPTR_T), value :: oblock
      real(C_DOUBLE), dimension(*), intent(out) :: in
      complex(C_DOUBLE_COMPLEX), dimension(*), intent(out) :: out
      integer(C_INT32_T), value :: comm
      integer(C_INT), value :: flags
    end function fftw_mpi_plan_many_dft_r2c
    
    type(C_PTR) function fftw_mpi_plan_dft_r2c(rnk,n,in,out,comm,flags) bind(C, name='fftw_mpi_plan_dft_r2c_f03')
      import
      integer(C_INT), value :: rnk
      integer(C_INTPTR_T), dimension(*), intent(in) :: n
      real(C_DOUBLE), dimension(*), intent(out) :: in
      complex(C_DOUBLE_COMPLEX), dimension(*), intent(out) :: out
      integer(C_INT32_T), value :: comm
      integer(C_INT), value :: flags
    end function fftw_mpi_plan_dft_r2c
    
    type(C_PTR) function fftw_mpi_plan_dft_r2c_2d(n0,n1,in,out,comm,flags) bind(C, name='fftw_mpi_plan_dft_r2c_2d_f03')
      import
      integer(C_INTPTR_T), value :: n0
      integer(C_INTPTR_T), value :: n1
      real(C_DOUBLE), dimension(*), intent(out) :: in
      complex(C_DOUBLE_COMPLEX), dimension(*), intent(out) :: out
      integer(C_INT32_T), value :: comm
      integer(C_INT), value :: flags
    end function fftw_mpi_plan_dft_r2c_2d
    
    type(C_PTR) function fftw_mpi_plan_dft_r2c_3d(n0,n1,n2,in,out,comm,flags) bind(C, name='fftw_mpi_plan_dft_r2c_3d_f03')
      import
      integer(C_INTPTR_T), value :: n0
      integer(C_INTPTR_T), value :: n1
      integer(C_INTPTR_T), value :: n2
      real(C_DOUBLE), dimension(*), intent(out) :: in
      complex(C_DOUBLE_COMPLEX), dimension(*), intent(out) :: out
      integer(C_INT32_T), value :: comm
      integer(C_INT), value :: flags
    end function fftw_mpi_plan_dft_r2c_3d
    
    type(C_PTR) function fftw_mpi_plan_many_dft_c2r(rnk,n,howmany,iblock,oblock,in,out,comm,flags) &
                         bind(C, name='fftw_mpi_plan_many_dft_c2r_f03')
      import
      integer(C_INT), value :: rnk
      integer(C_INTPTR_T), dimension(*), intent(in) :: n
      integer(C_INTPTR_T), value :: howmany
      integer(C_INTPTR_T), value :: iblock
      integer(C_INTPTR_T), value :: oblock
      complex(C_DOUBLE_COMPLEX), dimension(*), intent(out) :: in
      real(C_DOUBLE), dimension(*), intent(out) :: out
      integer(C_INT32_T), value :: comm
      integer(C_INT), value :: flags
    end function fftw_mpi_plan_many_dft_c2r
    
    type(C_PTR) function fftw_mpi_plan_dft_c2r(rnk,n,in,out,comm,flags) bind(C, name='fftw_mpi_plan_dft_c2r_f03')
      import
      integer(C_INT), value :: rnk
      integer(C_INTPTR_T), dimension(*), intent(in) :: n
      complex(C_DOUBLE_COMPLEX), dimension(*), intent(out) :: in
      real(C_DOUBLE), dimension(*), intent(out) :: out
      integer(C_INT32_T), value :: comm
      integer(C_INT), value :: flags
    end function fftw_mpi_plan_dft_c2r
    
    type(C_PTR) function fftw_mpi_plan_dft_c2r_2d(n0,n1,in,out,comm,flags) bind(C, name='fftw_mpi_plan_dft_c2r_2d_f03')
      import
      integer(C_INTPTR_T), value :: n0
      integer(C_INTPTR_T), value :: n1
      complex(C_DOUBLE_COMPLEX), dimension(*), intent(out) :: in
      real(C_DOUBLE), dimension(*), intent(out) :: out
      integer(C_INT32_T), value :: comm
      integer(C_INT), value :: flags
    end function fftw_mpi_plan_dft_c2r_2d
    
    type(C_PTR) function fftw_mpi_plan_dft_c2r_3d(n0,n1,n2,in,out,comm,flags) bind(C, name='fftw_mpi_plan_dft_c2r_3d_f03')
      import
      integer(C_INTPTR_T), value :: n0
      integer(C_INTPTR_T), value :: n1
      integer(C_INTPTR_T), value :: n2
      complex(C_DOUBLE_COMPLEX), dimension(*), intent(out) :: in
      real(C_DOUBLE), dimension(*), intent(out) :: out
      integer(C_INT32_T), value :: comm
      integer(C_INT), value :: flags
    end function fftw_mpi_plan_dft_c2r_3d
    
    subroutine fftw_mpi_gather_wisdom(comm_) bind(C, name='fftw_mpi_gather_wisdom_f03')
      import
      integer(C_INT32_T), value :: comm_
    end subroutine fftw_mpi_gather_wisdom
    
    subroutine fftw_mpi_broadcast_wisdom(comm_) bind(C, name='fftw_mpi_broadcast_wisdom_f03')
      import
      integer(C_INT32_T), value :: comm_
    end subroutine fftw_mpi_broadcast_wisdom
    
    subroutine fftw_mpi_execute_dft(p,in,out) bind(C, name='fftw_mpi_execute_dft')
      import
      type(C_PTR), value :: p
      complex(C_DOUBLE_COMPLEX), dimension(*), intent(inout) :: in
      complex(C_DOUBLE_COMPLEX), dimension(*), intent(out) :: out
    end subroutine fftw_mpi_execute_dft
    
    subroutine fftw_mpi_execute_dft_r2c(p,in,out) bind(C, name='fftw_mpi_execute_dft_r2c')
      import
      type(C_PTR), value :: p
      real(C_DOUBLE), dimension(*), intent(inout) :: in
      complex(C_DOUBLE_COMPLEX), dimension(*), intent(out) :: out
    end subroutine fftw_mpi_execute_dft_r2c
    
    subroutine fftw_mpi_execute_dft_c2r(p,in,out) bind(C, name='fftw_mpi_execute_dft_c2r')
      import
      type(C_PTR), value :: p
      complex(C_DOUBLE_COMPLEX), dimension(*), intent(inout) :: in
      real(C_DOUBLE), dimension(*), intent(out) :: out
    end subroutine fftw_mpi_execute_dft_c2r
    
    subroutine fftw_mpi_execute_r2r(p,in,out) bind(C, name='fftw_mpi_execute_r2r')
      import
      type(C_PTR), value :: p
      real(C_DOUBLE), dimension(*), intent(inout) :: in
      real(C_DOUBLE), dimension(*), intent(out) :: out
    end subroutine fftw_mpi_execute_r2r
    
  end interface

  type, bind(C) :: fftwf_mpi_ddim
     integer(C_INTPTR_T) n, ib, ob
  end type fftwf_mpi_ddim

  interface
    subroutine fftwf_mpi_init() bind(C, name='fftwf_mpi_init')
      import
    end subroutine fftwf_mpi_init
    
    subroutine fftwf_mpi_cleanup() bind(C, name='fftwf_mpi_cleanup')
      import
    end subroutine fftwf_mpi_cleanup
    
    integer(C_INTPTR_T) function fftwf_mpi_local_size_many_transposed(rnk,n,howmany,block0,block1,comm,local_n0,local_0_start, &
                                                                      local_n1,local_1_start) &
                                 bind(C, name='fftwf_mpi_local_size_many_transposed_f03')
      import
      integer(C_INT), value :: rnk
      integer(C_INTPTR_T), dimension(*), intent(in) :: n
      integer(C_INTPTR_T), value :: howmany
      integer(C_INTPTR_T), value :: block0
      integer(C_INTPTR_T), value :: block1
      integer(C_INT32_T), value :: comm
      integer(C_INTPTR_T), intent(out) :: local_n0
      integer(C_INTPTR_T), intent(out) :: local_0_start
      integer(C_INTPTR_T), intent(out) :: local_n1
      integer(C_INTPTR_T), intent(out) :: local_1_start
    end function fftwf_mpi_local_size_many_transposed
    
    integer(C_INTPTR_T) function fftwf_mpi_local_size_many(rnk,n,howmany,block0,comm,local_n0,local_0_start) &
                                 bind(C, name='fftwf_mpi_local_size_many_f03')
      import
      integer(C_INT), value :: rnk
      integer(C_INTPTR_T), dimension(*), intent(in) :: n
      integer(C_INTPTR_T), value :: howmany
      integer(C_INTPTR_T), value :: block0
      integer(C_INT32_T), value :: comm
      integer(C_INTPTR_T), intent(out) :: local_n0
      integer(C_INTPTR_T), intent(out) :: local_0_start
    end function fftwf_mpi_local_size_many
    
    integer(C_INTPTR_T) function fftwf_mpi_local_size_transposed(rnk,n,comm,local_n0,local_0_start,local_n1,local_1_start) &
                                 bind(C, name='fftwf_mpi_local_size_transposed_f03')
      import
      integer(C_INT), value :: rnk
      integer(C_INTPTR_T), dimension(*), intent(in) :: n
      integer(C_INT32_T), value :: comm
      integer(C_INTPTR_T), intent(out) :: local_n0
      integer(C_INTPTR_T), intent(out) :: local_0_start
      integer(C_INTPTR_T), intent(out) :: local_n1
      integer(C_INTPTR_T), intent(out) :: local_1_start
    end function fftwf_mpi_local_size_transposed
    
    integer(C_INTPTR_T) function fftwf_mpi_local_size(rnk,n,comm,local_n0,local_0_start) bind(C, name='fftwf_mpi_local_size_f03')
      import
      integer(C_INT), value :: rnk
      integer(C_INTPTR_T), dimension(*), intent(in) :: n
      integer(C_INT32_T), value :: comm
      integer(C_INTPTR_T), intent(out) :: local_n0
      integer(C_INTPTR_T), intent(out) :: local_0_start
    end function fftwf_mpi_local_size
    
    integer(C_INTPTR_T) function fftwf_mpi_local_size_many_1d(n0,howmany,comm,sign,flags,local_ni,local_i_start,local_no, &
                                                              local_o_start) bind(C, name='fftwf_mpi_local_size_many_1d_f03')
      import
      integer(C_INTPTR_T), value :: n0
      integer(C_INTPTR_T), value :: howmany
      integer(C_INT32_T), value :: comm
      integer(C_INT), value :: sign
      integer(C_INT), value :: flags
      integer(C_INTPTR_T), intent(out) :: local_ni
      integer(C_INTPTR_T), intent(out) :: local_i_start
      integer(C_INTPTR_T), intent(out) :: local_no
      integer(C_INTPTR_T), intent(out) :: local_o_start
    end function fftwf_mpi_local_size_many_1d
    
    integer(C_INTPTR_T) function fftwf_mpi_local_size_1d(n0,comm,sign,flags,local_ni,local_i_start,local_no,local_o_start) &
                                 bind(C, name='fftwf_mpi_local_size_1d_f03')
      import
      integer(C_INTPTR_T), value :: n0
      integer(C_INT32_T), value :: comm
      integer(C_INT), value :: sign
      integer(C_INT), value :: flags
      integer(C_INTPTR_T), intent(out) :: local_ni
      integer(C_INTPTR_T), intent(out) :: local_i_start
      integer(C_INTPTR_T), intent(out) :: local_no
      integer(C_INTPTR_T), intent(out) :: local_o_start
    end function fftwf_mpi_local_size_1d
    
    integer(C_INTPTR_T) function fftwf_mpi_local_size_2d(n0,n1,comm,local_n0,local_0_start) &
                                 bind(C, name='fftwf_mpi_local_size_2d_f03')
      import
      integer(C_INTPTR_T), value :: n0
      integer(C_INTPTR_T), value :: n1
      integer(C_INT32_T), value :: comm
      integer(C_INTPTR_T), intent(out) :: local_n0
      integer(C_INTPTR_T), intent(out) :: local_0_start
    end function fftwf_mpi_local_size_2d
    
    integer(C_INTPTR_T) function fftwf_mpi_local_size_2d_transposed(n0,n1,comm,local_n0,local_0_start,local_n1,local_1_start) &
                                 bind(C, name='fftwf_mpi_local_size_2d_transposed_f03')
      import
      integer(C_INTPTR_T), value :: n0
      integer(C_INTPTR_T), value :: n1
      integer(C_INT32_T), value :: comm
      integer(C_INTPTR_T), intent(out) :: local_n0
      integer(C_INTPTR_T), intent(out) :: local_0_start
      integer(C_INTPTR_T), intent(out) :: local_n1
      integer(C_INTPTR_T), intent(out) :: local_1_start
    end function fftwf_mpi_local_size_2d_transposed
    
    integer(C_INTPTR_T) function fftwf_mpi_local_size_3d(n0,n1,n2,comm,local_n0,local_0_start) &
                                 bind(C, name='fftwf_mpi_local_size_3d_f03')
      import
      integer(C_INTPTR_T), value :: n0
      integer(C_INTPTR_T), value :: n1
      integer(C_INTPTR_T), value :: n2
      integer(C_INT32_T), value :: comm
      integer(C_INTPTR_T), intent(out) :: local_n0
      integer(C_INTPTR_T), intent(out) :: local_0_start
    end function fftwf_mpi_local_size_3d
    
    integer(C_INTPTR_T) function fftwf_mpi_local_size_3d_transposed(n0,n1,n2,comm,local_n0,local_0_start,local_n1,local_1_start) &
                                 bind(C, name='fftwf_mpi_local_size_3d_transposed_f03')
      import
      integer(C_INTPTR_T), value :: n0
      integer(C_INTPTR_T), value :: n1
      integer(C_INTPTR_T), value :: n2
      integer(C_INT32_T), value :: comm
      integer(C_INTPTR_T), intent(out) :: local_n0
      integer(C_INTPTR_T), intent(out) :: local_0_start
      integer(C_INTPTR_T), intent(out) :: local_n1
      integer(C_INTPTR_T), intent(out) :: local_1_start
    end function fftwf_mpi_local_size_3d_transposed
    
    type(C_PTR) function fftwf_mpi_plan_many_transpose(n0,n1,howmany,block0,block1,in,out,comm,flags) &
                         bind(C, name='fftwf_mpi_plan_many_transpose_f03')
      import
      integer(C_INTPTR_T), value :: n0
      integer(C_INTPTR_T), value :: n1
      integer(C_INTPTR_T), value :: howmany
      integer(C_INTPTR_T), value :: block0
      integer(C_INTPTR_T), value :: block1
      real(C_FLOAT), dimension(*), intent(out) :: in
      real(C_FLOAT), dimension(*), intent(out) :: out
      integer(C_INT32_T), value :: comm
      integer(C_INT), value :: flags
    end function fftwf_mpi_plan_many_transpose
    
    type(C_PTR) function fftwf_mpi_plan_transpose(n0,n1,in,out,comm,flags) bind(C, name='fftwf_mpi_plan_transpose_f03')
      import
      integer(C_INTPTR_T), value :: n0
      integer(C_INTPTR_T), value :: n1
      real(C_FLOAT), dimension(*), intent(out) :: in
      real(C_FLOAT), dimension(*), intent(out) :: out
      integer(C_INT32_T), value :: comm
      integer(C_INT), value :: flags
    end function fftwf_mpi_plan_transpose
    
    type(C_PTR) function fftwf_mpi_plan_many_dft(rnk,n,howmany,block,tblock,in,out,comm,sign,flags) &
                         bind(C, name='fftwf_mpi_plan_many_dft_f03')
      import
      integer(C_INT), value :: rnk
      integer(C_INTPTR_T), dimension(*), intent(in) :: n
      integer(C_INTPTR_T), value :: howmany
      integer(C_INTPTR_T), value :: block
      integer(C_INTPTR_T), value :: tblock
      complex(C_FLOAT_COMPLEX), dimension(*), intent(out) :: in
      complex(C_FLOAT_COMPLEX), dimension(*), intent(out) :: out
      integer(C_INT32_T), value :: comm
      integer(C_INT), value :: sign
      integer(C_INT), value :: flags
    end function fftwf_mpi_plan_many_dft
    
    type(C_PTR) function fftwf_mpi_plan_dft(rnk,n,in,out,comm,sign,flags) bind(C, name='fftwf_mpi_plan_dft_f03')
      import
      integer(C_INT), value :: rnk
      integer(C_INTPTR_T), dimension(*), intent(in) :: n
      complex(C_FLOAT_COMPLEX), dimension(*), intent(out) :: in
      complex(C_FLOAT_COMPLEX), dimension(*), intent(out) :: out
      integer(C_INT32_T), value :: comm
      integer(C_INT), value :: sign
      integer(C_INT), value :: flags
    end function fftwf_mpi_plan_dft
    
    type(C_PTR) function fftwf_mpi_plan_dft_1d(n0,in,out,comm,sign,flags) bind(C, name='fftwf_mpi_plan_dft_1d_f03')
      import
      integer(C_INTPTR_T), value :: n0
      complex(C_FLOAT_COMPLEX), dimension(*), intent(out) :: in
      complex(C_FLOAT_COMPLEX), dimension(*), intent(out) :: out
      integer(C_INT32_T), value :: comm
      integer(C_INT), value :: sign
      integer(C_INT), value :: flags
    end function fftwf_mpi_plan_dft_1d
    
    type(C_PTR) function fftwf_mpi_plan_dft_2d(n0,n1,in,out,comm,sign,flags) bind(C, name='fftwf_mpi_plan_dft_2d_f03')
      import
      integer(C_INTPTR_T), value :: n0
      integer(C_INTPTR_T), value :: n1
      complex(C_FLOAT_COMPLEX), dimension(*), intent(out) :: in
      complex(C_FLOAT_COMPLEX), dimension(*), intent(out) :: out
      integer(C_INT32_T), value :: comm
      integer(C_INT), value :: sign
      integer(C_INT), value :: flags
    end function fftwf_mpi_plan_dft_2d
    
    type(C_PTR) function fftwf_mpi_plan_dft_3d(n0,n1,n2,in,out,comm,sign,flags) bind(C, name='fftwf_mpi_plan_dft_3d_f03')
      import
      integer(C_INTPTR_T), value :: n0
      integer(C_INTPTR_T), value :: n1
      integer(C_INTPTR_T), value :: n2
      complex(C_FLOAT_COMPLEX), dimension(*), intent(out) :: in
      complex(C_FLOAT_COMPLEX), dimension(*), intent(out) :: out
      integer(C_INT32_T), value :: comm
      integer(C_INT), value :: sign
      integer(C_INT), value :: flags
    end function fftwf_mpi_plan_dft_3d
    
    type(C_PTR) function fftwf_mpi_plan_many_r2r(rnk,n,howmany,iblock,oblock,in,out,comm,kind,flags) &
                         bind(C, name='fftwf_mpi_plan_many_r2r_f03')
      import
      integer(C_INT), value :: rnk
      integer(C_INTPTR_T), dimension(*), intent(in) :: n
      integer(C_INTPTR_T), value :: howmany
      integer(C_INTPTR_T), value :: iblock
      integer(C_INTPTR_T), value :: oblock
      real(C_FLOAT), dimension(*), intent(out) :: in
      real(C_FLOAT), dimension(*), intent(out) :: out
      integer(C_INT32_T), value :: comm
      integer(C_FFTW_R2R_KIND), dimension(*), intent(in) :: kind
      integer(C_INT), value :: flags
    end function fftwf_mpi_plan_many_r2r
    
    type(C_PTR) function fftwf_mpi_plan_r2r(rnk,n,in,out,comm,kind,flags) bind(C, name='fftwf_mpi_plan_r2r_f03')
      import
      integer(C_INT), value :: rnk
      integer(C_INTPTR_T), dimension(*), intent(in) :: n
      real(C_FLOAT), dimension(*), intent(out) :: in
      real(C_FLOAT), dimension(*), intent(out) :: out
      integer(C_INT32_T), value :: comm
      integer(C_FFTW_R2R_KIND), dimension(*), intent(in) :: kind
      integer(C_INT), value :: flags
    end function fftwf_mpi_plan_r2r
    
    type(C_PTR) function fftwf_mpi_plan_r2r_2d(n0,n1,in,out,comm,kind0,kind1,flags) bind(C, name='fftwf_mpi_plan_r2r_2d_f03')
      import
      integer(C_INTPTR_T), value :: n0
      integer(C_INTPTR_T), value :: n1
      real(C_FLOAT), dimension(*), intent(out) :: in
      real(C_FLOAT), dimension(*), intent(out) :: out
      integer(C_INT32_T), value :: comm
      integer(C_FFTW_R2R_KIND), value :: kind0
      integer(C_FFTW_R2R_KIND), value :: kind1
      integer(C_INT), value :: flags
    end function fftwf_mpi_plan_r2r_2d
    
    type(C_PTR) function fftwf_mpi_plan_r2r_3d(n0,n1,n2,in,out,comm,kind0,kind1,kind2,flags) &
                         bind(C, name='fftwf_mpi_plan_r2r_3d_f03')
      import
      integer(C_INTPTR_T), value :: n0
      integer(C_INTPTR_T), value :: n1
      integer(C_INTPTR_T), value :: n2
      real(C_FLOAT), dimension(*), intent(out) :: in
      real(C_FLOAT), dimension(*), intent(out) :: out
      integer(C_INT32_T), value :: comm
      integer(C_FFTW_R2R_KIND), value :: kind0
      integer(C_FFTW_R2R_KIND), value :: kind1
      integer(C_FFTW_R2R_KIND), value :: kind2
      integer(C_INT), value :: flags
    end function fftwf_mpi_plan_r2r_3d
    
    type(C_PTR) function fftwf_mpi_plan_many_dft_r2c(rnk,n,howmany,iblock,oblock,in,out,comm,flags) &
                         bind(C, name='fftwf_mpi_plan_many_dft_r2c_f03')
      import
      integer(C_INT), value :: rnk
      integer(C_INTPTR_T), dimension(*), intent(in) :: n
      integer(C_INTPTR_T), value :: howmany
      integer(C_INTPTR_T), value :: iblock
      integer(C_INTPTR_T), value :: oblock
      real(C_FLOAT), dimension(*), intent(out) :: in
      complex(C_FLOAT_COMPLEX), dimension(*), intent(out) :: out
      integer(C_INT32_T), value :: comm
      integer(C_INT), value :: flags
    end function fftwf_mpi_plan_many_dft_r2c
    
    type(C_PTR) function fftwf_mpi_plan_dft_r2c(rnk,n,in,out,comm,flags) bind(C, name='fftwf_mpi_plan_dft_r2c_f03')
      import
      integer(C_INT), value :: rnk
      integer(C_INTPTR_T), dimension(*), intent(in) :: n
      real(C_FLOAT), dimension(*), intent(out) :: in
      complex(C_FLOAT_COMPLEX), dimension(*), intent(out) :: out
      integer(C_INT32_T), value :: comm
      integer(C_INT), value :: flags
    end function fftwf_mpi_plan_dft_r2c
    
    type(C_PTR) function fftwf_mpi_plan_dft_r2c_2d(n0,n1,in,out,comm,flags) bind(C, name='fftwf_mpi_plan_dft_r2c_2d_f03')
      import
      integer(C_INTPTR_T), value :: n0
      integer(C_INTPTR_T), value :: n1
      real(C_FLOAT), dimension(*), intent(out) :: in
      complex(C_FLOAT_COMPLEX), dimension(*), intent(out) :: out
      integer(C_INT32_T), value :: comm
      integer(C_INT), value :: flags
    end function fftwf_mpi_plan_dft_r2c_2d
    
    type(C_PTR) function fftwf_mpi_plan_dft_r2c_3d(n0,n1,n2,in,out,comm,flags) bind(C, name='fftwf_mpi_plan_dft_r2c_3d_f03')
      import
      integer(C_INTPTR_T), value :: n0
      integer(C_INTPTR_T), value :: n1
      integer(C_INTPTR_T), value :: n2
      real(C_FLOAT), dimension(*), intent(out) :: in
      complex(C_FLOAT_COMPLEX), dimension(*), intent(out) :: out
      integer(C_INT32_T), value :: comm
      integer(C_INT), value :: flags
    end function fftwf_mpi_plan_dft_r2c_3d
    
    type(C_PTR) function fftwf_mpi_plan_many_dft_c2r(rnk,n,howmany,iblock,oblock,in,out,comm,flags) &
                         bind(C, name='fftwf_mpi_plan_many_dft_c2r_f03')
      import
      integer(C_INT), value :: rnk
      integer(C_INTPTR_T), dimension(*), intent(in) :: n
      integer(C_INTPTR_T), value :: howmany
      integer(C_INTPTR_T), value :: iblock
      integer(C_INTPTR_T), value :: oblock
      complex(C_FLOAT_COMPLEX), dimension(*), intent(out) :: in
      real(C_FLOAT), dimension(*), intent(out) :: out
      integer(C_INT32_T), value :: comm
      integer(C_INT), value :: flags
    end function fftwf_mpi_plan_many_dft_c2r
    
    type(C_PTR) function fftwf_mpi_plan_dft_c2r(rnk,n,in,out,comm,flags) bind(C, name='fftwf_mpi_plan_dft_c2r_f03')
      import
      integer(C_INT), value :: rnk
      integer(C_INTPTR_T), dimension(*), intent(in) :: n
      complex(C_FLOAT_COMPLEX), dimension(*), intent(out) :: in
      real(C_FLOAT), dimension(*), intent(out) :: out
      integer(C_INT32_T), value :: comm
      integer(C_INT), value :: flags
    end function fftwf_mpi_plan_dft_c2r
    
    type(C_PTR) function fftwf_mpi_plan_dft_c2r_2d(n0,n1,in,out,comm,flags) bind(C, name='fftwf_mpi_plan_dft_c2r_2d_f03')
      import
      integer(C_INTPTR_T), value :: n0
      integer(C_INTPTR_T), value :: n1
      complex(C_FLOAT_COMPLEX), dimension(*), intent(out) :: in
      real(C_FLOAT), dimension(*), intent(out) :: out
      integer(C_INT32_T), value :: comm
      integer(C_INT), value :: flags
    end function fftwf_mpi_plan_dft_c2r_2d
    
    type(C_PTR) function fftwf_mpi_plan_dft_c2r_3d(n0,n1,n2,in,out,comm,flags) bind(C, name='fftwf_mpi_plan_dft_c2r_3d_f03')
      import
      integer(C_INTPTR_T), value :: n0
      integer(C_INTPTR_T), value :: n1
      integer(C_INTPTR_T), value :: n2
      complex(C_FLOAT_COMPLEX), dimension(*), intent(out) :: in
      real(C_FLOAT), dimension(*), intent(out) :: out
      integer(C_INT32_T), value :: comm
      integer(C_INT), value :: flags
    end function fftwf_mpi_plan_dft_c2r_3d
    
    subroutine fftwf_mpi_gather_wisdom(comm_) bind(C, name='fftwf_mpi_gather_wisdom_f03')
      import
      integer(C_INT32_T), value :: comm_
    end subroutine fftwf_mpi_gather_wisdom
    
    subroutine fftwf_mpi_broadcast_wisdom(comm_) bind(C, name='fftwf_mpi_broadcast_wisdom_f03')
      import
      integer(C_INT32_T), value :: comm_
    end subroutine fftwf_mpi_broadcast_wisdom
    
    subroutine fftwf_mpi_execute_dft(p,in,out) bind(C, name='fftwf_mpi_execute_dft')
      import
      type(C_PTR), value :: p
      complex(C_FLOAT_COMPLEX), dimension(*), intent(inout) :: in
      complex(C_FLOAT_COMPLEX), dimension(*), intent(out) :: out
    end subroutine fftwf_mpi_execute_dft
    
    subroutine fftwf_mpi_execute_dft_r2c(p,in,out) bind(C, name='fftwf_mpi_execute_dft_r2c')
      import
      type(C_PTR), value :: p
      real(C_FLOAT), dimension(*), intent(inout) :: in
      complex(C_FLOAT_COMPLEX), dimension(*), intent(out) :: out
    end subroutine fftwf_mpi_execute_dft_r2c
    
    subroutine fftwf_mpi_execute_dft_c2r(p,in,out) bind(C, name='fftwf_mpi_execute_dft_c2r')
      import
      type(C_PTR), value :: p
      complex(C_FLOAT_COMPLEX), dimension(*), intent(inout) :: in
      real(C_FLOAT), dimension(*), intent(out) :: out
    end subroutine fftwf_mpi_execute_dft_c2r
    
    subroutine fftwf_mpi_execute_r2r(p,in,out) bind(C, name='fftwf_mpi_execute_r2r')
      import
      type(C_PTR), value :: p
      real(C_FLOAT), dimension(*), intent(inout) :: in
      real(C_FLOAT), dimension(*), intent(out) :: out
    end subroutine fftwf_mpi_execute_r2r
    
  end interface

  type, bind(C) :: fftwl_mpi_ddim
     integer(C_INTPTR_T) n, ib, ob
  end type fftwl_mpi_ddim

  interface
    subroutine fftwl_mpi_init() bind(C, name='fftwl_mpi_init')
      import
    end subroutine fftwl_mpi_init
    
    subroutine fftwl_mpi_cleanup() bind(C, name='fftwl_mpi_cleanup')
      import
    end subroutine fftwl_mpi_cleanup
    
    integer(C_INTPTR_T) function fftwl_mpi_local_size_many_transposed(rnk,n,howmany,block0,block1,comm,local_n0,local_0_start, &
                                                                      local_n1,local_1_start) &
                                 bind(C, name='fftwl_mpi_local_size_many_transposed_f03')
      import
      integer(C_INT), value :: rnk
      integer(C_INTPTR_T), dimension(*), intent(in) :: n
      integer(C_INTPTR_T), value :: howmany
      integer(C_INTPTR_T), value :: block0
      integer(C_INTPTR_T), value :: block1
      integer(C_INT32_T), value :: comm
      integer(C_INTPTR_T), intent(out) :: local_n0
      integer(C_INTPTR_T), intent(out) :: local_0_start
      integer(C_INTPTR_T), intent(out) :: local_n1
      integer(C_INTPTR_T), intent(out) :: local_1_start
    end function fftwl_mpi_local_size_many_transposed
    
    integer(C_INTPTR_T) function fftwl_mpi_local_size_many(rnk,n,howmany,block0,comm,local_n0,local_0_start) &
                                 bind(C, name='fftwl_mpi_local_size_many_f03')
      import
      integer(C_INT), value :: rnk
      integer(C_INTPTR_T), dimension(*), intent(in) :: n
      integer(C_INTPTR_T), value :: howmany
      integer(C_INTPTR_T), value :: block0
      integer(C_INT32_T), value :: comm
      integer(C_INTPTR_T), intent(out) :: local_n0
      integer(C_INTPTR_T), intent(out) :: local_0_start
    end function fftwl_mpi_local_size_many
    
    integer(C_INTPTR_T) function fftwl_mpi_local_size_transposed(rnk,n,comm,local_n0,local_0_start,local_n1,local_1_start) &
                                 bind(C, name='fftwl_mpi_local_size_transposed_f03')
      import
      integer(C_INT), value :: rnk
      integer(C_INTPTR_T), dimension(*), intent(in) :: n
      integer(C_INT32_T), value :: comm
      integer(C_INTPTR_T), intent(out) :: local_n0
      integer(C_INTPTR_T), intent(out) :: local_0_start
      integer(C_INTPTR_T), intent(out) :: local_n1
      integer(C_INTPTR_T), intent(out) :: local_1_start
    end function fftwl_mpi_local_size_transposed
    
    integer(C_INTPTR_T) function fftwl_mpi_local_size(rnk,n,comm,local_n0,local_0_start) bind(C, name='fftwl_mpi_local_size_f03')
      import
      integer(C_INT), value :: rnk
      integer(C_INTPTR_T), dimension(*), intent(in) :: n
      integer(C_INT32_T), value :: comm
      integer(C_INTPTR_T), intent(out) :: local_n0
      integer(C_INTPTR_T), intent(out) :: local_0_start
    end function fftwl_mpi_local_size
    
    integer(C_INTPTR_T) function fftwl_mpi_local_size_many_1d(n0,howmany,comm,sign,flags,local_ni,local_i_start,local_no, &
                                                              local_o_start) bind(C, name='fftwl_mpi_local_size_many_1d_f03')
      import
      integer(C_INTPTR_T), value :: n0
      integer(C_INTPTR_T), value :: howmany
      integer(C_INT32_T), value :: comm
      integer(C_INT), value :: sign
      integer(C_INT), value :: flags
      integer(C_INTPTR_T), intent(out) :: local_ni
      integer(C_INTPTR_T), intent(out) :: local_i_start
      integer(C_INTPTR_T), intent(out) :: local_no
      integer(C_INTPTR_T), intent(out) :: local_o_start
    end function fftwl_mpi_local_size_many_1d
    
    integer(C_INTPTR_T) function fftwl_mpi_local_size_1d(n0,comm,sign,flags,local_ni,local_i_start,local_no,local_o_start) &
                                 bind(C, name='fftwl_mpi_local_size_1d_f03')
      import
      integer(C_INTPTR_T), value :: n0
      integer(C_INT32_T), value :: comm
      integer(C_INT), value :: sign
      integer(C_INT), value :: flags
      integer(C_INTPTR_T), intent(out) :: local_ni
      integer(C_INTPTR_T), intent(out) :: local_i_start
      integer(C_INTPTR_T), intent(out) :: local_no
      integer(C_INTPTR_T), intent(out) :: local_o_start
    end function fftwl_mpi_local_size_1d
    
    integer(C_INTPTR_T) function fftwl_mpi_local_size_2d(n0,n1,comm,local_n0,local_0_start) &
                                 bind(C, name='fftwl_mpi_local_size_2d_f03')
      import
      integer(C_INTPTR_T), value :: n0
      integer(C_INTPTR_T), value :: n1
      integer(C_INT32_T), value :: comm
      integer(C_INTPTR_T), intent(out) :: local_n0
      integer(C_INTPTR_T), intent(out) :: local_0_start
    end function fftwl_mpi_local_size_2d
    
    integer(C_INTPTR_T) function fftwl_mpi_local_size_2d_transposed(n0,n1,comm,local_n0,local_0_start,local_n1,local_1_start) &
                                 bind(C, name='fftwl_mpi_local_size_2d_transposed_f03')
      import
      integer(C_INTPTR_T), value :: n0
      integer(C_INTPTR_T), value :: n1
      integer(C_INT32_T), value :: comm
      integer(C_INTPTR_T), intent(out) :: local_n0
      integer(C_INTPTR_T), intent(out) :: local_0_start
      integer(C_INTPTR_T), intent(out) :: local_n1
      integer(C_INTPTR_T), intent(out) :: local_1_start
    end function fftwl_mpi_local_size_2d_transposed
    
    integer(C_INTPTR_T) function fftwl_mpi_local_size_3d(n0,n1,n2,comm,local_n0,local_0_start) &
                                 bind(C, name='fftwl_mpi_local_size_3d_f03')
      import
      integer(C_INTPTR_T), value :: n0
      integer(C_INTPTR_T), value :: n1
      integer(C_INTPTR_T), value :: n2
      integer(C_INT32_T), value :: comm
      integer(C_INTPTR_T), intent(out) :: local_n0
      integer(C_INTPTR_T), intent(out) :: local_0_start
    end function fftwl_mpi_local_size_3d
    
    integer(C_INTPTR_T) function fftwl_mpi_local_size_3d_transposed(n0,n1,n2,comm,local_n0,local_0_start,local_n1,local_1_start) &
                                 bind(C, name='fftwl_mpi_local_size_3d_transposed_f03')
      import
      integer(C_INTPTR_T), value :: n0
      integer(C_INTPTR_T), value :: n1
      integer(C_INTPTR_T), value :: n2
      integer(C_INT32_T), value :: comm
      integer(C_INTPTR_T), intent(out) :: local_n0
      integer(C_INTPTR_T), intent(out) :: local_0_start
      integer(C_INTPTR_T), intent(out) :: local_n1
      integer(C_INTPTR_T), intent(out) :: local_1_start
    end function fftwl_mpi_local_size_3d_transposed
    
    type(C_PTR) function fftwl_mpi_plan_many_transpose(n0,n1,howmany,block0,block1,in,out,comm,flags) &
                         bind(C, name='fftwl_mpi_plan_many_transpose_f03')
      import
      integer(C_INTPTR_T), value :: n0
      integer(C_INTPTR_T), value :: n1
      integer(C_INTPTR_T), value :: howmany
      integer(C_INTPTR_T), value :: block0
      integer(C_INTPTR_T), value :: block1
      real(C_LONG_DOUBLE), dimension(*), intent(out) :: in
      real(C_LONG_DOUBLE), dimension(*), intent(out) :: out
      integer(C_INT32_T), value :: comm
      integer(C_INT), value :: flags
    end function fftwl_mpi_plan_many_transpose
    
    type(C_PTR) function fftwl_mpi_plan_transpose(n0,n1,in,out,comm,flags) bind(C, name='fftwl_mpi_plan_transpose_f03')
      import
      integer(C_INTPTR_T), value :: n0
      integer(C_INTPTR_T), value :: n1
      real(C_LONG_DOUBLE), dimension(*), intent(out) :: in
      real(C_LONG_DOUBLE), dimension(*), intent(out) :: out
      integer(C_INT32_T), value :: comm
      integer(C_INT), value :: flags
    end function fftwl_mpi_plan_transpose
    
    type(C_PTR) function fftwl_mpi_plan_many_dft(rnk,n,howmany,block,tblock,in,out,comm,sign,flags) &
                         bind(C, name='fftwl_mpi_plan_many_dft_f03')
      import
      integer(C_INT), value :: rnk
      integer(C_INTPTR_T), dimension(*), intent(in) :: n
      integer(C_INTPTR_T), value :: howmany
      integer(C_INTPTR_T), value :: block
      integer(C_INTPTR_T), value :: tblock
      complex(C_LONG_DOUBLE_COMPLEX), dimension(*), intent(out) :: in
      complex(C_LONG_DOUBLE_COMPLEX), dimension(*), intent(out) :: out
      integer(C_INT32_T), value :: comm
      integer(C_INT), value :: sign
      integer(C_INT), value :: flags
    end function fftwl_mpi_plan_many_dft
    
    type(C_PTR) function fftwl_mpi_plan_dft(rnk,n,in,out,comm,sign,flags) bind(C, name='fftwl_mpi_plan_dft_f03')
      import
      integer(C_INT), value :: rnk
      integer(C_INTPTR_T), dimension(*), intent(in) :: n
      complex(C_LONG_DOUBLE_COMPLEX), dimension(*), intent(out) :: in
      complex(C_LONG_DOUBLE_COMPLEX), dimension(*), intent(out) :: out
      integer(C_INT32_T), value :: comm
      integer(C_INT), value :: sign
      integer(C_INT), value :: flags
    end function fftwl_mpi_plan_dft
    
    type(C_PTR) function fftwl_mpi_plan_dft_1d(n0,in,out,comm,sign,flags) bind(C, name='fftwl_mpi_plan_dft_1d_f03')
      import
      integer(C_INTPTR_T), value :: n0
      complex(C_LONG_DOUBLE_COMPLEX), dimension(*), intent(out) :: in
      complex(C_LONG_DOUBLE_COMPLEX), dimension(*), intent(out) :: out
      integer(C_INT32_T), value :: comm
      integer(C_INT), value :: sign
      integer(C_INT), value :: flags
    end function fftwl_mpi_plan_dft_1d
    
    type(C_PTR) function fftwl_mpi_plan_dft_2d(n0,n1,in,out,comm,sign,flags) bind(C, name='fftwl_mpi_plan_dft_2d_f03')
      import
      integer(C_INTPTR_T), value :: n0
      integer(C_INTPTR_T), value :: n1
      complex(C_LONG_DOUBLE_COMPLEX), dimension(*), intent(out) :: in
      complex(C_LONG_DOUBLE_COMPLEX), dimension(*), intent(out) :: out
      integer(C_INT32_T), value :: comm
      integer(C_INT), value :: sign
      integer(C_INT), value :: flags
    end function fftwl_mpi_plan_dft_2d
    
    type(C_PTR) function fftwl_mpi_plan_dft_3d(n0,n1,n2,in,out,comm,sign,flags) bind(C, name='fftwl_mpi_plan_dft_3d_f03')
      import
      integer(C_INTPTR_T), value :: n0
      integer(C_INTPTR_T), value :: n1
      integer(C_INTPTR_T), value :: n2
      complex(C_LONG_DOUBLE_COMPLEX), dimension(*), intent(out) :: in
      complex(C_LONG_DOUBLE_COMPLEX), dimension(*), intent(out) :: out
      integer(C_INT32_T), value :: comm
      integer(C_INT), value :: sign
      integer(C_INT), value :: flags
    end function fftwl_mpi_plan_dft_3d
    
    type(C_PTR) function fftwl_mpi_plan_many_r2r(rnk,n,howmany,iblock,oblock,in,out,comm,kind,flags) &
                         bind(C, name='fftwl_mpi_plan_many_r2r_f03')
      import
      integer(C_INT), value :: rnk
      integer(C_INTPTR_T), dimension(*), intent(in) :: n
      integer(C_INTPTR_T), value :: howmany
      integer(C_INTPTR_T), value :: iblock
      integer(C_INTPTR_T), value :: oblock
      real(C_LONG_DOUBLE), dimension(*), intent(out) :: in
      real(C_LONG_DOUBLE), dimension(*), intent(out) :: out
      integer(C_INT32_T), value :: comm
      integer(C_FFTW_R2R_KIND), dimension(*), intent(in) :: kind
      integer(C_INT), value :: flags
    end function fftwl_mpi_plan_many_r2r
    
    type(C_PTR) function fftwl_mpi_plan_r2r(rnk,n,in,out,comm,kind,flags) bind(C, name='fftwl_mpi_plan_r2r_f03')
      import
      integer(C_INT), value :: rnk
      integer(C_INTPTR_T), dimension(*), intent(in) :: n
      real(C_LONG_DOUBLE), dimension(*), intent(out) :: in
      real(C_LONG_DOUBLE), dimension(*), intent(out) :: out
      integer(C_INT32_T), value :: comm
      integer(C_FFTW_R2R_KIND), dimension(*), intent(in) :: kind
      integer(C_INT), value :: flags
    end function fftwl_mpi_plan_r2r
    
    type(C_PTR) function fftwl_mpi_plan_r2r_2d(n0,n1,in,out,comm,kind0,kind1,flags) bind(C, name='fftwl_mpi_plan_r2r_2d_f03')
      import
      integer(C_INTPTR_T), value :: n0
      integer(C_INTPTR_T), value :: n1
      real(C_LONG_DOUBLE), dimension(*), intent(out) :: in
      real(C_LONG_DOUBLE), dimension(*), intent(out) :: out
      integer(C_INT32_T), value :: comm
      integer(C_FFTW_R2R_KIND), value :: kind0
      integer(C_FFTW_R2R_KIND), value :: kind1
      integer(C_INT), value :: flags
    end function fftwl_mpi_plan_r2r_2d
    
    type(C_PTR) function fftwl_mpi_plan_r2r_3d(n0,n1,n2,in,out,comm,kind0,kind1,kind2,flags) &
                         bind(C, name='fftwl_mpi_plan_r2r_3d_f03')
      import
      integer(C_INTPTR_T), value :: n0
      integer(C_INTPTR_T), value :: n1
      integer(C_INTPTR_T), value :: n2
      real(C_LONG_DOUBLE), dimension(*), intent(out) :: in
      real(C_LONG_DOUBLE), dimension(*), intent(out) :: out
      integer(C_INT32_T), value :: comm
      integer(C_FFTW_R2R_KIND), value :: kind0
      integer(C_FFTW_R2R_KIND), value :: kind1
      integer(C_FFTW_R2R_KIND), value :: kind2
      integer(C_INT), value :: flags
    end function fftwl_mpi_plan_r2r_3d
    
    type(C_PTR) function fftwl_mpi_plan_many_dft_r2c(rnk,n,howmany,iblock,oblock,in,out,comm,flags) &
                         bind(C, name='fftwl_mpi_plan_many_dft_r2c_f03')
      import
      integer(C_INT), value :: rnk
      integer(C_INTPTR_T), dimension(*), intent(in) :: n
      integer(C_INTPTR_T), value :: howmany
      integer(C_INTPTR_T), value :: iblock
      integer(C_INTPTR_T), value :: oblock
      real(C_LONG_DOUBLE), dimension(*), intent(out) :: in
      complex(C_LONG_DOUBLE_COMPLEX), dimension(*), intent(out) :: out
      integer(C_INT32_T), value :: comm
      integer(C_INT), value :: flags
    end function fftwl_mpi_plan_many_dft_r2c
    
    type(C_PTR) function fftwl_mpi_plan_dft_r2c(rnk,n,in,out,comm,flags) bind(C, name='fftwl_mpi_plan_dft_r2c_f03')
      import
      integer(C_INT), value :: rnk
      integer(C_INTPTR_T), dimension(*), intent(in) :: n
      real(C_LONG_DOUBLE), dimension(*), intent(out) :: in
      complex(C_LONG_DOUBLE_COMPLEX), dimension(*), intent(out) :: out
      integer(C_INT32_T), value :: comm
      integer(C_INT), value :: flags
    end function fftwl_mpi_plan_dft_r2c
    
    type(C_PTR) function fftwl_mpi_plan_dft_r2c_2d(n0,n1,in,out,comm,flags) bind(C, name='fftwl_mpi_plan_dft_r2c_2d_f03')
      import
      integer(C_INTPTR_T), value :: n0
      integer(C_INTPTR_T), value :: n1
      real(C_LONG_DOUBLE), dimension(*), intent(out) :: in
      complex(C_LONG_DOUBLE_COMPLEX), dimension(*), intent(out) :: out
      integer(C_INT32_T), value :: comm
      integer(C_INT), value :: flags
    end function fftwl_mpi_plan_dft_r2c_2d
    
    type(C_PTR) function fftwl_mpi_plan_dft_r2c_3d(n0,n1,n2,in,out,comm,flags) bind(C, name='fftwl_mpi_plan_dft_r2c_3d_f03')
      import
      integer(C_INTPTR_T), value :: n0
      integer(C_INTPTR_T), value :: n1
      integer(C_INTPTR_T), value :: n2
      real(C_LONG_DOUBLE), dimension(*), intent(out) :: in
      complex(C_LONG_DOUBLE_COMPLEX), dimension(*), intent(out) :: out
      integer(C_INT32_T), value :: comm
      integer(C_INT), value :: flags
    end function fftwl_mpi_plan_dft_r2c_3d
    
    type(C_PTR) function fftwl_mpi_plan_many_dft_c2r(rnk,n,howmany,iblock,oblock,in,out,comm,flags) &
                         bind(C, name='fftwl_mpi_plan_many_dft_c2r_f03')
      import
      integer(C_INT), value :: rnk
      integer(C_INTPTR_T), dimension(*), intent(in) :: n
      integer(C_INTPTR_T), value :: howmany
      integer(C_INTPTR_T), value :: iblock
      integer(C_INTPTR_T), value :: oblock
      complex(C_LONG_DOUBLE_COMPLEX), dimension(*), intent(out) :: in
      real(C_LONG_DOUBLE), dimension(*), intent(out) :: out
      integer(C_INT32_T), value :: comm
      integer(C_INT), value :: flags
    end function fftwl_mpi_plan_many_dft_c2r
    
    type(C_PTR) function fftwl_mpi_plan_dft_c2r(rnk,n,in,out,comm,flags) bind(C, name='fftwl_mpi_plan_dft_c2r_f03')
      import
      integer(C_INT), value :: rnk
      integer(C_INTPTR_T), dimension(*), intent(in) :: n
      complex(C_LONG_DOUBLE_COMPLEX), dimension(*), intent(out) :: in
      real(C_LONG_DOUBLE), dimension(*), intent(out) :: out
      integer(C_INT32_T), value :: comm
      integer(C_INT), value :: flags
    end function fftwl_mpi_plan_dft_c2r
    
    type(C_PTR) function fftwl_mpi_plan_dft_c2r_2d(n0,n1,in,out,comm,flags) bind(C, name='fftwl_mpi_plan_dft_c2r_2d_f03')
      import
      integer(C_INTPTR_T), value :: n0
      integer(C_INTPTR_T), value :: n1
      complex(C_LONG_DOUBLE_COMPLEX), dimension(*), intent(out) :: in
      real(C_LONG_DOUBLE), dimension(*), intent(out) :: out
      integer(C_INT32_T), value :: comm
      integer(C_INT), value :: flags
    end function fftwl_mpi_plan_dft_c2r_2d
    
    type(C_PTR) function fftwl_mpi_plan_dft_c2r_3d(n0,n1,n2,in,out,comm,flags) bind(C, name='fftwl_mpi_plan_dft_c2r_3d_f03')
      import
      integer(C_INTPTR_T), value :: n0
      integer(C_INTPTR_T), value :: n1
      integer(C_INTPTR_T), value :: n2
      complex(C_LONG_DOUBLE_COMPLEX), dimension(*), intent(out) :: in
      real(C_LONG_DOUBLE), dimension(*), intent(out) :: out
      integer(C_INT32_T), value :: comm
      integer(C_INT), value :: flags
    end function fftwl_mpi_plan_dft_c2r_3d
    
    subroutine fftwl_mpi_gather_wisdom(comm_) bind(C, name='fftwl_mpi_gather_wisdom_f03')
      import
      integer(C_INT32_T), value :: comm_
    end subroutine fftwl_mpi_gather_wisdom
    
    subroutine fftwl_mpi_broadcast_wisdom(comm_) bind(C, name='fftwl_mpi_broadcast_wisdom_f03')
      import
      integer(C_INT32_T), value :: comm_
    end subroutine fftwl_mpi_broadcast_wisdom
    
    subroutine fftwl_mpi_execute_dft(p,in,out) bind(C, name='fftwl_mpi_execute_dft')
      import
      type(C_PTR), value :: p
      complex(C_LONG_DOUBLE_COMPLEX), dimension(*), intent(inout) :: in
      complex(C_LONG_DOUBLE_COMPLEX), dimension(*), intent(out) :: out
    end subroutine fftwl_mpi_execute_dft
    
    subroutine fftwl_mpi_execute_dft_r2c(p,in,out) bind(C, name='fftwl_mpi_execute_dft_r2c')
      import
      type(C_PTR), value :: p
      real(C_LONG_DOUBLE), dimension(*), intent(inout) :: in
      complex(C_LONG_DOUBLE_COMPLEX), dimension(*), intent(out) :: out
    end subroutine fftwl_mpi_execute_dft_r2c
    
    subroutine fftwl_mpi_execute_dft_c2r(p,in,out) bind(C, name='fftwl_mpi_execute_dft_c2r')
      import
      type(C_PTR), value :: p
      complex(C_LONG_DOUBLE_COMPLEX), dimension(*), intent(inout) :: in
      real(C_LONG_DOUBLE), dimension(*), intent(out) :: out
    end subroutine fftwl_mpi_execute_dft_c2r
    
    subroutine fftwl_mpi_execute_r2r(p,in,out) bind(C, name='fftwl_mpi_execute_r2r')
      import
      type(C_PTR), value :: p
      real(C_LONG_DOUBLE), dimension(*), intent(inout) :: in
      real(C_LONG_DOUBLE), dimension(*), intent(out) :: out
    end subroutine fftwl_mpi_execute_r2r
    
  end interface
