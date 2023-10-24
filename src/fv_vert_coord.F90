!=================================================================================|
! Set up sigma and more, to be continued....                                                       !
!								                                                  !
! case(1): SIGMA LEVELS                                                           !
! Sigma levels are determined by a formula of                                     !
!                       sigma(k)=-[(k-1)/(kb-1)]^pow_sig                          !
!    pow_sig=1: uniform sigma layers                                              !
!    pow_sig=2: layers satisfying a parabolic function with high                  !
!               vertical resolution near the surface and bottom.                  !
!    pow_sig can be used any real number                                          !
!									                                              !
! case(2): GENERAL VERTICAL LEVEL                                                 !
! Vertical levels are determined by the formula                                   !
!                                                                                 !
!          tanh[(lev_bot+lev_sur)*(nsigma-k)/(nsigma-1)-lev_bot]+tanh(lev_bot)    !
! sigma(k)= ------------------------------------------------------------------ - 1!
!                      tanh(lev_bot) + tanh(lev_sur)                              !
!                                                                                 !
!                                                                                 !
! case(3): CONSTANT LAYER TRANSFORMATION                                          !
! It is not a real sigma, this is useful approach for the sediment solutions      !
! Near bottom or/and near surface you have layers with constant depth.            !
! If the depth is too large, than you have case of ind_sigma=2. This model        !
! is useful for some task in frame sediment movement simulation.                  !
! For this case you need prescribe the number of these constant layers            !
! and its thickness. Also you prescribe the critical depth, from which you will   !
! shift to usual sigma distribution (ind_sigma=2).If the depth is less than sum   !
! of the depths of 'constant layers', than the uniform sigma layers               !
! will be prescribed. Useful refernce is:                                         !
! Pietrzak, Jakobson, Burchard, Vested, Petersen, 2002. A three-dimensional       !
! hydrostatic model for coastal and ocean modelling using a generalised topography!
! following co-ordinate system. Ocean Modelling 4, 173-205                        !
!=================================================================================|
subroutine set_deepsigma(coeff, d)
  USE o_MESH
  USE o_ARRAYS
  USE o_PARAM

  USE g_parsup
   !subrountine to get depth-sigma levels distribution for fang-eddy case
   !coeff - coefficient to reduce all resolutions, used on second step to adapt to requested nsigma
   ! coeff =1 for nsigma=352
   !ssize - preliminary number of sigma levels
   !d - depth distribution
   !dnum - final number of sigma
   !--parameters: r1,r2,r3 - resolutions, d1,d2,d3 - corresponding depths d3 - depth
   ! example of use:
   !     integer :: dnum
   !     real(wp) :: d(1000)
   !     call set_deepsigma(1.0_wp,1000, d,dnum)
   !     call set_deepsigma(dnum/(nsigma-3),nsigma, d,dnum)
   !     sigma = 1-d(1:nsigma)/4400

   implicit none
   real(WP),intent(in) :: coeff
   integer,parameter :: ssize=1000
!   integer,intent(out) :: dnum
   real(WP),intent(out) :: d(ssize)
   integer :: i, id1, nn
   real(WP) :: r1,r2,r3,d1,d2,d3,a,b,rn3,x1,x2,x3,lb,ls
   r1=0.5_WP
   r2=4.0_WP*coeff
   r3=250.0_WP*coeff
   d1=120.0_WP
   d2=150.0_WP
   d3=4400.0_WP
   lb=0.3_WP
   ls=2.0_WP
   d(1)=0.0_WP
   i=1
   do while (d(i)<d1)
      i=i+1
      if (i>ssize) stop "set_deepsigma err 1, ssize< sigmalev"
      d(i)=d(i-1)+r1
   end do
!   a=(r2-r1)/(d2-d(i))
!   b=r1-a*d(i)
!   do while (d(i)<d2)
!      i=i+1
!      if (i>ssize) stop "set_deepsigma err 2, ssize< sigmalev"
!      d(i)=d(i-1)+a*d(i-1)+b
!   end do
   id1=i
   nn=nsigma-id1
   if (mype==0) print *,"id1=",id1,d(id1),nn
!   do i=1,nn
!        x1=tanh((lb+ls)*(nn-i)/(nn) - lb)
!        x2=tanh(lb)
!        x3=x2+tanh(ls)
!        d(i+id1)=d(id1)+(d3-d(id1))*(1-(x1+x2)/x3)
!        if (mype==0) print *,"i,d=",i,i+id1,d(i+id1)
!   enddo
   if (mod(nn,2) == 0) then
      d(id1+1)=d(id1)+r1
      id1=id1+1
   end if
   nn=nsigma-id1
   if (mype==0) print *,"id1=",id1,d(id1),nn,(nn+1)/2

   do i=1,(nn+1)/2
      d(i+id1)=r1+d(id1)+(d3-d(id1))*&
                 ( ((i-1)/float((nn+1)/2-1))**pow_sig*0.5_WP   )
       !if (mype==0) print *,"FIRSTPART", i, i+id1-1, d(i+id1-1)
   enddo
   do i=(nn+1)/2+1,nn
      d(i+id1)=r1+d(id1)+(d3-d(id1))*&
                 ( 1 -((nn-i)/float((nn+1)/2-1))**pow_sig*0.5_WP   )
       !if (mype==0) print *,"SECONDPART", i, i+id1-1, d(i+id1-1)
       !if (mype==0) print *, i, nn-i,((nn-i)/float((nn+1)/2-1)),( 1 -((nn-i)/float((nn+1)/2-1))**pow_sig*0.5_WP   )
   enddo

   if (mype==0) then
      !print *, "N rest levels: ",nn
   end if
   d(nsigma) = d3

!   a=(r3-r2)/(d3-d(i))
!   b=r2-a*d(i)
!   do while (d(i)<d3)
!      i=i+1
!      if (i>ssize) then
!         print *, r1,r2,r3
!         print *, "set_deepsigma err 3, ssize< sigmalev", d(i-1),i,ssize
!         print *, d(1:10)
!         print *, d(90:100)
!         stop
!      end if
!      d(i)=d(i-1)+a*d(i-1)+b
!   end do
!   d(i) = d3
!   dnum = nsigma

end subroutine set_deepsigma

subroutine SET_SIGMA

  !=================================================================================|
  USE o_MESH
  USE o_ARRAYS
  USE o_PARAM

  USE g_parsup

  implicit none

  integer :: K,KK,lb,ls
  integer :: I
  real(WP):: X1,X2,X3
  real(WP):: DR,CL,CU, CB, CS, D_ref_min
  real(WP), allocatable   :: ZKBw(:), ZKSw(:)

  integer :: dnum
  real(wp) :: d(1000),coeff


  !==============================================================================|

  if (mype==0) print *,'Entering SET_SIGMA'
  select case (ind_sigma)
  case(4)
     ! special case for experiment fang_eddy
     coeff = 352.0_wp/float(nsigma)
     call set_deepsigma(coeff, d)
     sigma = 1-d(1:nsigma)/4400.0_wp

!     call set_deepsigma(1.0_wp,1000, d,dnum)
!     coeff = float(dnum)/float(nsigma)
print *, mype, '------------------------------------'
print *, mype, d(1:5)
print *, mype, d(nsigma-5:nsigma)
!print *, mype, coeff,dnum, nsigma
!     call set_deepsigma(coeff,nsigma, d,dnum)
!     sigma = 1-d(1:nsigma)/4400.0_wp
  case(1)

     allocate(sigma(nsigma))

     IF(pow_sig > 1 .AND. MOD(nsigma,2) == 0)THEN
        WRITE(*,*) 'nsigma should be an odd number,stop ....'
        STOP
     END IF

     IF(pow_sig == 1)THEN
        DO K=1,nsigma
           sigma(K) = -((K-1)/DFLOAT(nsigma-1))**pow_sig + 1.0_WP
        END DO
     ELSE
        DO K=1,(nsigma+1)/2
           sigma(K) = -((K-1)/DFLOAT((nsigma+1)/2-1))**pow_sig*0.5_WP +1.0_WP
        END DO

        DO K=(nsigma+1)/2+1,nsigma
           sigma(K) = ((nsigma-K)/DFLOAT((nsigma+1)/2-1))**pow_sig*0.5_WP
        END DO
     END IF

     !==============================================================================|

  case(2)

     allocate(sigma(nsigma))
     ! Sigma
     DO K=1,nsigma
        X1=TANH((lev_bot+lev_sur)*(nsigma-K)/(nsigma-1) - lev_bot)
        X2=TANH(lev_bot)
        X3=X2+TANH(lev_sur)
        sigma(K)=(X1+X2)/X3
     END DO

     !==============================================================================|

  case(3)

     allocate(Zl(nsigma,nod2D))

     !  IF((sum(KSw)+sum(KBw))>D_ref_min)THEN
     !     WRITE(*,*) '(sum(KSw)+sum(KBw)) should be less then or equal to D_ref_min....'
     !    STOP
     !  END IF

     D_ref_min=sum(KSw)+sum(KBw)

     DO I=1,nod2D
        IF(depth(I) < D_ref_min)THEN
           DO K=1,nsigma
              Zl(K,I)=-((K-1)/DFLOAT(nsigma-1)) + 1.0_WP
           END DO
        ELSE
           IF(depth(I) > D_ref_max)THEN
              lev_bot=1.0_WP+size(KBw)/nsigma
              lev_sur=1.0_WP+size(KSw)/nsigma

              DO K=1,nsigma
                 X1=TANH((lev_bot+lev_sur)*(nsigma-K)/(nsigma-1) - lev_bot)
                 X2=TANH(lev_bot)
                 X3=X2+TANH(lev_sur)
                 Zl(K,I)=(X1+X2)/X3
              END DO

           ELSE
              lb=size(KBw)
              ls=size(KSw)
              allocate(ZKSw(ls),ZKBw(lb))
              CS=-sum(KSw)/depth(I)
              CB=sum(KBw)/depth(I)-1.0_WP
              DR=(CB-CS)/(nsigma-ls-lb-1)

              DO K=1,ls
                 ZKSw(K)=KSw(K)/depth(I)
              END DO
              DO K=1,lb
                 ZKBw(K)=KBw(K)/depth(I)
              END DO

              Zl(1,I)=1.0_WP

              DO K=2,ls+1
                 Zl(K,I)=Zl(K-1,I)-ZKSw(K-1)
              END DO

              DO K=ls+2,nsigma-lb
                 Zl(K,I)=Zl(K-1,I)+DR
              END DO

              KK=0
              DO K=nsigma-lb+1,nsigma
                 KK=KK+1
                 Zl(K,I)=Zl(K-1,I)-ZKBw(KK)
              END DO
              Zl(nsigma,I)=0.0_WP !can be very small, but non-zero value, after the loop above
              deallocate(ZKSw,ZKBw)
           END IF
        END IF
     END DO
  end select

  if (mype==0) print *,'... SET_SIGMA done'

  return
end subroutine SET_SIGMA
!==============================================================================|
