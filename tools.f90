module tools
    use types
    implicit none

    public :: trapintng,gauss

    contains

  function trapintng(a,b,step,f)result(volume)

    use types
    use functions

    real(kind=rk), intent(in) :: a,b
    real(kind=rk) :: volume,dx,value1,value2
    integer, intent(in) :: step
    integer :: i

     interface
        function f(x) result(value)
        use types
         real(kind=rk), intent(in) :: x
         real(kind=rk) :: value
        end function f
     end interface

    volume = 0_rk
    value1 = 0_rk
    value2 = 0_rk

    dx = abs((b-a)/step)

    Do i=0,step-1
            value1 = value1 + f(a+(i*dx))
    End do

    Do i=1,step
            value2 = value2 + f(a+(i*dx))
    End do

        volume = volume + (dx/2) * (value1+value2)

    print*, "the volume of integral is:",volume

  end function trapintng


  subroutine gauss(a,b,x)
    implicit none
    real (kind=rk) :: c

    real (kind=rk), dimension (:,:), intent (in),allocatable :: a
    real (kind=rk), dimension (:), intent (out),allocatable :: x
    real (kind=rk), dimension (:), intent (in),allocatable :: b

    real (kind=rk), dimension (:,:), allocatable :: acopy
    real (kind=rk), dimension (:), allocatable :: bcopy

    integer :: i, j, k, n

    if (.not. allocated(a) .or. .not.(allocated(b))) then
        print *, "inputs to gem are not allocated..."
        return
    end if

    if ((ubound(a,1) /= ubound(a,2))) then
        print*, "matrix a is not square..."
        return
    end if

    if ((ubound(a,1) /= ubound(b,1))) then
        print*, "vector b and matrix a have incompatible lengths..."
        return
    end if

    if (.not. allocated (x)) then
        allocate(x(ubound(a,1)))
    end if

    allocate(acopy(ubound(a,1), ubound(a,1)))
    acopy(:, 1:ubound(a,1)) = a
    !acopy(:, ubound(a,1)+1) = b
    allocate(bcopy(ubound(b,1)))
    bcopy(1:ubound(b,1)) = b

    print*,"Original matrix"
    print*,"-------------------------------------------------"

    call printmatrix(acopy,(ubound(a,1)),(ubound(a,1)))
    print*,"-------------------------------------------------"

    n = ubound(a,1)

    !forward elimination

!    do k=1, n-1
!        do i=k+1,n
!            c=acopy(i,k)/acopy(k,k)
!            acopy(i,k) = 0.0
!            bcopy(i)=bcopy(i)- c*bcopy(k)
!                do j=k+1,n
!                    acopy(i,j) = acopy(i,j)-c*acopy(k,j)
!                end do
!        end do
!    end do
        do j=1, n-1
            pos = maxloc(abs(acopy(j:n,j)),1) + j-1
            rowtemp = acopy(pos,:)
            acopy(pos,:) = acopy(j,:)
            acopy(j,:) = rowtemp

!            do i=j=1, n
!                acopy(i,1) = acopy(i,1) * acopy(j,1) / acopy(j,j) * acopy(i,j)
!            end do

        end do

        do j=2, ubound(a,1)
          do i = j, ubound(a,1)
            acopy(i, :) = acopy(i,:) - acopy(j-1,:)/acopy(j-1,j-1)*acopy(i,j-1)
          end do
        end do

    !back substitution

    x(n) = bcopy(n)/acopy(n,n)
    do i=n-1,1,-1
        c=0.0
        do j=i+1,n
            c= c + acopy(i,j)*x(j)
        end do
        x(i) = (bcopy(i)- c)/acopy(i,i)
    end do

    print*,"Transformed matrix"
    print*,"-------------------------------------------------"

    call printmatrix(acopy,(ubound(a,1)),(ubound(a,1)))
    print*,"-------------------------------------------------"


  end subroutine gauss

  subroutine printmatrix(a,n_rows,n_cols)
    implicit none
    integer, intent(in) :: n_rows,n_cols
    real (kind=rk), intent(in)    :: a(n_rows,n_cols)
    integer rows,cols

    do rows=1,n_rows
        do cols=1,n_cols
            write(*,1000,advance='no')a(rows,cols)
        end do
        write(*,1010)
    end do

    1000  format(f10.6)
    1010  format()

    return
  end subroutine printmatrix

end module tools
