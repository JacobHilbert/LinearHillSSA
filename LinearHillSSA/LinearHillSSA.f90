! f2py --f90flags="-Ofast -fbackslash" -c LinearHillSSA.f90 -m ssa 
module LinearHillSSA
    implicit none
    !integer, parameter :: dp = 4
    character*1, parameter :: creturn = achar(13)
    
    contains
    
    ! hill rate function v*(x^n)/(x^n+alpha^n)
    ! use negative n for repression
    pure elemental function hill(X,v,alpha,n) result(r)
        real(8), intent(in) :: X,v,alpha,n
        real(8) :: r
        r = v*(x**n)/(x**n+alpha**n)
    end function hill
    
    ! propensity function with 
    ! * linear relation matrix A, 
    ! * constant term b, 
    ! * Hill relation matrix M,
    ! * Hill constants, multipliers, and coefficeints v, alpha, n
    pure function propensity(x,A,b,M,v,alpha,n) result(p)
        real(8), intent(in) :: x(:),A(:,:),b(:),M(:,:),v(:),alpha(:),n(:)
        real(8) :: p(size(b))
        p = matmul(A,x) + b + hill(matmul(M,x),v,alpha,n)
    end function propensity
    
    pure function jacobian(t,x,A,b,M,v,alpha,n,R) result(dx)
        real(8), intent(in) :: t,x(:),A(:,:),b(:),M(:,:),v(:),alpha(:),n(:),R(:,:)
        real(8) :: dx(size(x))
        dx = matmul(transpose(R),propensity(x,A,b,M,v,alpha,n))
    end function jacobian
    
    ! cummulative sum of real array
    pure function cumsum(a) result(r)
        real(8), intent(in) :: a(:)
        real(8) :: r(size(a))
        integer :: i
        r = [(sum(a(1:i)), i=1,size(a))]
    end function cumsum
    
    ! returns the index of where element should be inserted into array,
    ! such that the resultant array is sorted.
    ! (Assumes array is already sorted)
    pure function searchsortedfirst(array,element) result(i)
        real(8), intent(in) :: array(:), element 
        integer :: i, n
        n = size(array)
        do i = 1,n 
            if (element < array(i)) then 
                exit
            end if
        end do
    end function searchsortedfirst
    
    ! Performs one step of the SSA. 
    ! MODIFIES x.
    ! returns the time elapsed until a reaction happened.
    function SSA_step(x,A,b,M,v,alpha,n,R) result(tau)
        real(8), intent(inout) :: x(:)
        real(8), intent(in) :: A(:,:),b(:),M(:,:),v(:),alpha(:),n(:),R(:,:)
        real(8) :: tau, p(size(b)), p_cumsum(size(b)), p_sum, u
        integer :: i
        p = propensity(x,A,b,M,v,alpha,n)
        p_sum = sum(p)
        p_cumsum = cumsum(p)
        call random_number(u)
        i = searchsortedfirst(p_cumsum,u*p_sum)
        x = x + R(i,:)
        tau = -log(u)/p_sum
    end function SSA_step
    
    ! gives a raw path of the jump process
    ! mainly for purposes of studing the resulting times distribution
    ! MODIFIES x0, X, and times.
    subroutine SSA_path(x0,X,times,A,b,M,v,alpha,n,R)
        real(8), intent(inout) :: x0(:), X(:,:), times(:)
        real(8), intent(in) :: A(:,:),b(:),M(:,:),v(:),alpha(:),n(:),R(:,:)
        real(8) :: t
        integer :: i 
        t = 0.0
        do i=1,size(times)
            times(i) = t
            X(:,i) = x0
            t = t + SSA_step(x0,A,b,M,v,alpha,n,R)
        end do
    end subroutine SSA_path
    
    ! realizes the jump stochastic process, and measures its value  at specified times.
    ! MODIFIES x0 and X.
    ! X(:,i) is the resultant state vector at time times(i).
    subroutine SSA_measure(x0,times,X,A,b,M,v,alpha,n,R)
        real(8), intent(in) :: times(:), A(:,:),b(:),M(:,:),v(:),alpha(:),n(:),R(:,:)
        real(8), intent(inout) :: x0(:),X(:,:)
        real(8) :: t
        integer :: i
        t = 0.0
        do i=1,size(times)
            do while (t<times(i))
                t = t+SSA_step(x0,A,b,M,v,alpha,n,R)
            end do 
            X(:,i) = x0 
        end do
    end subroutine SSA_measure
    
    ! performs multiple realizations of the jump stochastic process
    ! system j is started from state vector x0(:,j)
    ! MODIFIES x0 and X
    ! X(:,i,j) is the resultant state vector for system j and time times(i).
    ! this subroutine prints 
    subroutine SSA_ensemble(x0,times,X,A,b,M,v,alpha,n,R)
        real(8), intent(in) :: times(:), A(:,:),b(:),M(:,:),v(:),alpha(:),n(:),R(:,:)
        real(8), intent(inout) :: x0(:,:),X(:,:,:)
        integer j,j_max
        j_max = size(x0,dim=2)
        do j = 1,j_max
            write(*,'(a,i5)', advance='no') creturn, j
            call SSA_measure(x0(:,j),times,X(:,:,j),A,b,M,v,alpha,n,R)
        end do
        print *,""
    end subroutine SSA_ensemble
    
end module LinearHillSSA