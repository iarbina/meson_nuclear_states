!***********************************************************************
! Linear spline interpolation
!
! - Takes as input a grid (x, y)_i for i in 1 to n, and the point
!   where the interpolation is required z
!
! - Returns the s(z) being s the linear spline
!***********************************************************************
subroutine spls (x, y, ndim, z, s)

    implicit none

    ! define double precision
    integer, parameter :: dp = selected_real_kind (15, 307)

    ! input variables
    integer, intent(in) :: ndim
    real(dp), dimension(ndim), intent(in) :: x, y
    real(dp), intent(in) :: z

    ! output variable
    real(dp), intent(out) :: s

    ! internal variables
    integer :: i, gdim         ! dimension of the grid

    !-------------------------------------------------------------------

    ! Throw error if size of the arrays of the grid are different
    if (size(x) /= size(y)) then
        print *, 'size(x) =', size(x), 'size(y) =', size(y)
        print *, 'ERROR (spls): Sizes of the input arrays are different!'
        STOP
    else
        gdim = size(x)
    end if

    do i = 1, gdim

        if (x(i) > x(i+1) .and. i < gdim) then
            print *, 'ERROR (spls): Independent variable array not ordered in increasing values', x(i), x(i+1)
            STOP
        end if
        
        if (x(i) >= z) go to 876
    
    end do
876 continue

    if (i == 1) then
        print *, 'ERROR (spls): Point out of grid range. Extrapolation required'
PRINT *, 'Input =', z, 'First array point =', x(i), 'Index =', i
        !STOP
    else
        s = ( y(i-1) * ( x(i) - z ) + y(i) * ( z - x(i-1) ) ) / ( x(i) - x(i-1) )   
    end if

end subroutine spls
