program LinearEquationSolver
  implicit none
  
  ! Matrix and vector declaration
  integer, parameter :: n = 3       
  integer :: i, j, k
  real(8), allocatable :: a(:,:), b(:), x(:)
  
  ! factor multiplier and 
  real(8) :: factor, sum_val
    
  ! a program that shows the stability and reliability of the Gaussian
  ! elimination to solve systems of linear equations   

  ! Initialization of the 
  allocate(a(n,n), b(n), x(n))
  
  a(1, :) = (/2, 0, 1/)
  a(2, :) = (/0, 4, 1/)
  a(3, :) = (/1, 1, 2/)
  b = (/30, 40, 30/)
  x = 0

  ! Initial System Visualization
  write(*, '(a)') 'Initial System A*x = b:'
  do i = 1, n
      write(*, '(3(f10.4), " | ", f10.4)') (a(i, j), j=1, n), b(i)
  enddo

  ! --- Gaussian Elimination (Forward Pass for upper triangular form) ---
  do k = 1, n - 1 ! pivot k 
      do i = k + 1, n ! row selection  
          if (abs(a(k, k)) < 1.0d-12) then 
              write(*, '(a)') 'Error: Near-zero pivot detected.'
              goto 100
          endif
          
          factor = a(i, k) / a(k, k) ! m(i,k) elimination factor 
          b(i) = b(i) - factor * b(k) ! simultaneous numerical swap to b(i)' in order to update the full augmented matrix

   ! This loop subtracts the factored pivot row(k) from row(i) across all columns(j) 
          do j = k, n 
              a(i, j) = a(i, j) - factor * a(k, j)! row transformation 
          enddo
      enddo
  enddo
  
  ! --- Back-Substitution ---
  do i = n, 1, -1 
      sum_val = 0 ! clears the temporary sum for each new row calculation 

      do j = i + 1, n
          sum_val = sum_val + a(i, j) * x(j)
      enddo
      
      if (abs(a(i, i)) < 1.0d-12) then
          write(*, '(a)') 'Error: Zero diagonal element.'
          goto 100
      endif

      x(i) = (b(i) - sum_val) / a(i, i)
  enddo
  
  ! --- Output Results ---
  write(*, '(a)') 'Final Upper Triangular Matrix A (U):'
  do i = 1, n
      write(*, '(3f10.4)') (a(i, j), j=1, n)
  enddo
  write(*, *)

  write(*, '(a)') 'Solution Vector x:'
  do i = 1, n
      write(*, '(a, i1, a, f15.8)') 'x(', i, ') = ', x(i)
  enddo
  
  ! --- Deallocation and Error Stop ---
  100 continue 
  if (allocated(a)) deallocate(a, b, x)
  
end program LinearEquationSolver
