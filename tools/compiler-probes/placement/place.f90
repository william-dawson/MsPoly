program place
  use, intrinsic :: iso_c_binding
  use mpi_f08
  implicit none
  interface
     integer(c_int) function place_init(r) bind(C)
       import :: c_int
       integer(c_int), value :: r
     end function
     real(c_double) function place_read(x, n) bind(C)
       import :: c_double, c_long
       real(c_double), dimension(*) :: x
       integer(c_long), value :: n
     end function
     real(c_double) function place_gpu_write(x, n) bind(C)
       import :: c_double, c_long
       real(c_double), dimension(*) :: x
       integer(c_long), value :: n
     end function
     real(c_double) function place_prefetch(x, n) bind(C)
       import :: c_double, c_long
       real(c_double), dimension(*) :: x
       integer(c_long), value :: n
     end function
     real(c_double) function place_copy_from_device(x, n) bind(C)
       import :: c_double, c_long
       real(c_double), dimension(*) :: x
       integer(c_long), value :: n
     end function
  end interface
  integer(c_long), parameter :: n = 2_c_long**25
  real(c_double), allocatable :: a(:), b(:), c(:)
  real(c_double) :: t1, t2
  integer :: rank, device
  character(len=200) :: line
  call MPI_Init()
  call MPI_Comm_rank(MPI_COMM_WORLD, rank)
  device = place_init(rank)
  allocate(a(n), b(n), c(n))
  t1 = place_gpu_write(a, n); t2 = place_read(a, n)
  write (line, '(A,I0,A,I0,A,F8.2,A,F8.1,A)') "rank ", rank, " gpu ", device, &
       " (a) GPU first touch: write ", 1d3*t1, " ms, read ", 8d0*n/t2/1d9, " GB/s"
  print '(A)', trim(line)
  b = 1
  t1 = place_prefetch(b, n); t2 = place_read(b, n)
  write (line, '(A,I0,A,F8.2,A,F8.1,A)') "rank ", rank, &
       " (b) CPU fill + prefetch: prefetch ", 1d3*t1, " ms, read ", 8d0*n/t2/1d9, " GB/s"
  print '(A)', trim(line)
  t1 = place_copy_from_device(c, n); t2 = place_read(c, n)
  write (line, '(A,I0,A,F8.2,A,F8.1,A)') "rank ", rank, &
       " (c) copy from device: copy ", 1d3*t1, " ms, read ", 8d0*n/t2/1d9, " GB/s"
  print '(A)', trim(line)
  call MPI_Finalize()
end program
