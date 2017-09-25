Ccccc This is the new Kalman IG filter which tries to eliminate errors from data sets adopting 
Ccccc Information Geometry techniques.

      Program Kalman

      implicit none

      integer dim,i,j,k,m,n,Status,hours, 
     & hour_obs, history_index, index

      real E,L1,L2,a1,S1,S2,S3,b1,Covar,Var,diff,yV,xV,
     +temp,A,B,C,Tobs,Tmod1,Tmod2,Uobs,Umod, 
     +D,W,V,H,P,PR,S,T,x,y,z,ss,x_matrix,Bias, 
     +Q1,P1,KG,ID,pre,Mean,Obs,Model,KAL,ErrorModel,ErrorKAL,Sum,Q,
     +Sumsquare,MeanError,MeanAbsError,StDevError,StDevAbsError,RMSE,
     +Probab,Sumabs,aa, history_index_r, index_r

      Parameter (dim=2)  
      Parameter (hours=24)  
      Parameter (history_index=7)  
      Parameter (index=history_index-1)  
      character*12 date_obs, date_mod1, date_mod2
     
      dimension L1(0:index),L2(0:index),E(dim,0:index),a1(0:index),
     +b1(0:index),diff(dim,dim),yV(0:index),xV(dim,0:index),
     +temp(dim,dim),A(dim,dim),
     +C(dim,dim),H(dim,dim),PR(dim,dim),x_matrix(dim,0:history_index), 
     +KG(dim,dim),S(dim,dim),T(dim,dim),x(dim,dim),Q(dim,dim), 
     +P1(dim,dim),ID(dim,dim),B(dim,dim),Obs(9000), 
     +Model(9000),KAL(9000),ErrorModel(9000),ErrorKAL(9000),
     +W(dim,dim), P(dim,dim),
     +Tobs(hours), Tmod1(hours), Tmod2(hours),
     +date_mod1(hours), date_mod2(hours), date_obs(hours)

      character*100 none 

      do m=1,dim
      do n=1,dim
         if (m.eq.n) then
          ID(m,n)=1
         else 
          ID(m,n)=0
         endif 
      enddo
      enddo

      history_index_r = history_index*1.0
      index_r = index*1.0

      write (*,*) 'history_index=', history_index
      write (*,*) 'index=', index
      write (*,*) 'history_index_r=', history_index_r
      write (*,*) 'index_r=', index_r

         open (unit=10, file='./obs.txt',status='old')    ! obs
         open (unit=15, file='./model1.txt',status='old')  ! model's yest fcst
         open (unit=18, file='./model2.txt',status='old')  ! model's new fcst

        do i=1,hours 
          read (10,*) Tobs(i)
        enddo
 1000 format(a12,1x,f2.0) 

        do i=1,hours
          read (18,*) Tmod2(i)
        enddo

        do i=1,hours 
          read (15,*) Tmod1(i) 
        enddo
 1500 format(a12,1x,f7.2) 

 
      DO i=1,hours   ! Loop over number of  hours 
       write (*,*) 'HOUR=', i 

      pre=Tmod1(i)

      open (unit=20, file='Kalman.txt', status='old')

105      read (20,*,ioStat=Status) none
      if (Status.lt.0) then 
        go to 106
      else 
        go to 105
      endif

106   continue

         ss = Tmod1(i)
         z  = Tobs(i)

110     y=z-ss

        do n=1,dim
         H(1,n)=ss**(n-1)
        enddo

         open (unit=16, file='yV.txt', status='old')
         
          do j=0,index 
           read(16,*) yV(j)
          enddo

         open (unit=17, file='x_matrix.txt', status='old')
          do j=0,history_index 
             read(17,*) (x_matrix(m,j), m=1,dim)
          enddo

          do j=1,history_index
          do m=1,dim 
            xV(m,j) = x_matrix(m,j) - x_matrix(m,j-1)
          enddo 
          enddo 

!!!!!!!	Define Variance !!!!!!!!!!!!!!!!!!!!!!!
	open (unit=69, file='IG_val', status='old')

	read(69,*) V
	V = V
	
	close (69)
       open (unit=19, file='P_matrix.txt', status='old')


          do n=1,dim 
             read(19,*) (P(m,n), m=1,dim)
          enddo

         open (unit=21, file='x.txt', status='old')
         read(21,*) (x(m,1), m=1,dim)


         do m=1,dim
         do n=1,dim
          call Line(xV,dim,m,L1,index)
          call Line(xV,dim,n,L2,index)
        call Covariance(L1,L2,S1,S2,S3,W(m,n),dim,history_index,index)
         enddo
         enddo


Ccccccccc KALMAN GAIN cccccccccccccccccccccccccccccccccccccccccccccc
Ccccccc    here we calculate the value of the Kalman gain     ccccccc

120     do m=1,dim
        do n=1,dim     
         P1(m,n) = P(m,n) + W(m,n)
        enddo
        enddo

        call sf(H,P1,Q,dim)
        Q1=Q(1,1)
        call mT(H,T,dim)
        call mprod(P1,T,PR,dim)

        do m=1,dim
        do n=1,dim
         KG(m,n)=PR(m,n)/(Q1+V)
        enddo
        enddo


Cccccccccc CALCULATION OF THE NEW VALUE ccccccccccccccccccccccc
Cccc    here we calculate the new value of the system equation, 
Ccccc   i.e the output of the program 

        call mprod(H,x,PR,dim)
        D=y-PR(1,1)

        call mprod(KG,D,PR,dim)

        
        do m=1,dim
        do n=1,dim
         S(m,n)=x(m,n)+PR(m,n)
         temp(m,n)=x(m,n) 
         x(m,n)=S(m,n)
         diff(m,n)=x(m,n)-temp(m,n)
        enddo
        enddo

         if (dim.eq.1) then 
          pre=x(1,1)
         else 
          pre=x(1,1)+x(2,1)*Tmod1(i)
         endif

         if (dim.ge.3) then  
          do m=3,dim
           pre=pre+x(m,1)*(Tmod1(i)**(m-1))
          enddo
         endif

        Bias = pre 
        pre = pre + Tmod2(i)

        if ((abs(pre-Tmod2(i)).gt.20.)) then
         pre = Tmod2(i)
        endif
        if (pre.lt.0.) then
         pre = 0
        endif



       write (20,*) Tobs(i), Tmod1(i), Tmod2(i), pre



3000  FORMAT(f8.2,2x,f8.2)

Ccccccccccccccc UPDATE OF COEFFICIENTS ccccccccccccccccccccccccc

       call mprod(KG,H,PR,dim)
       call mdif(ID,PR,D,dim)
       call mprod(D,P1,PR,dim)
      
        do m=1,dim
        do n=1,dim
          P(m,n)=PR(m,n)
        enddo
        enddo

Cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


150     close (10)
        close (15)
        close (16)
        close (17)
        close (18)
        close (19)
        close (21)

C_IG         call system ("rm -f ./yV.txt")
         call system ("rm -f ./x_matrix.txt")
         call system ("rm -f ./P_matrix.txt")
         call system ("rm -f ./x.txt")

CCC update of yV file CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

C_IG         open (unit=16, file='yV.txt', status='new')
C_IG          do j=1,index 
C_IG          write (16,*) yV(j)
C_IG          enddo 
C_IG  
C_IG          yV(index) = Tobs(i)-Tmod1(i)
C_IG 
C_IG          do m=1,dim
C_IG           aa = Tmod1(i)
C_IG           yV(index) = yV(index)-x(m,1)*(aa**(m-1))
C_IG          enddo           
C_IG
C_IG          write (16,*) yV(index) 
C_IG         close (16)

CcccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccC

Ccccccccccc update of x_matrix file  ccccccccccccccccccccccccccccccccccC
         open (unit=17, file='x_matrix.txt', status='new')
          do j=1,history_index
           write (17,*) (x_matrix(m,j), m=1,dim) 
          enddo 
           write (17,*) (x(m,1), m=1,dim) 
         close (17)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

Ccccccccccc update of P_matrix file  ccccccccccccccccccccccccccccccccccC
        open (unit=19, file='P_matrix.txt', status='new')
         do n=1,dim
            write(19,*) (P(m,n), m=1,dim)
         enddo
        close (19)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

Ccccccccccc update of x  file   ccccccccccccccccccccccccccccccccccC
         open (unit=21, file='x.txt', status='new')
         write (21,*) (x(m,1), m=1,dim)
         close (21)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

100    continue      
       END DO   ! end do loop over days 

        stop
        end 

Cccccccccccccc TRANSPOSED MATRIX cccccccccccccccccccccccccccccccccccccc
Cccccccccc   T = the transpose of A  cccccccccccccccccccccccccccccccccc

        subroutine mT(A,T,dim)
        integer dim
        real A, T  
        dimension A(dim,dim), T(dim,dim) 

        do m=1,dim
        do n=1,dim
          T(m,n)=A(n,m)
        enddo
        enddo
      
        return
        end 

Cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
Cccccccccccccc PRODUCT OF MATRICES ccccccccccccccccccccccccccccccccccccc
Ccccccccccccc   PR = A * B         ccccccccccccccccccccccccccccccccccccc
    
        subroutine mprod(A,B,PR,dim)
        integer dim
        real A,B,PR
        DIMENSION A(dim,dim), B(dim,dim), PR(dim,dim)
    
        do m=1,dim
        do n=1,dim
         PR(m,n)=0
        enddo
        enddo
    
        do m=1,dim
        do n=1,dim
        do k=1,dim
         PR(m,n)=PR(m,n)+A(m,k)*B(k,n)
        enddo
        enddo
        enddo
    
        return 
        end 

Ccccccccccc TRIPLE PRODUCT  ccccccccccccccccccccccccccccccccccccccccccccc
Cccccccc    Q = A * B * TR(A) cccccccccccccccccccccccccccccccc
    
        subroutine sf(A,B,Q,dim)
        integer dim
        real A, B, C, Q, PR, T
      dimension A(dim,dim),B(dim,dim),C(dim,dim),Q(dim,dim),PR(dim,dim),
     +T(dim,dim)
      
       call mprod(A,B,PR,dim)
      
       do m=1,dim
       do n=1,dim
          C(m,n)=PR(m,n)
       enddo
       enddo 
     
      call mT(A,T,dim)
      call mprod(C,T,PR,dim)

       do m=1,dim
       do n=1,dim
          Q(m,n)=PR(m,n)
       enddo
       enddo 
    
       return 
       end 

Ccccccccccccccccc DIFFERENCE OF MATRICES cccccccccccccccccccccccccccccccccccccc     
Cccccccccccccccc       D = A - B          ccccccccccccccccccccccccccccccccccc
    
      subroutine mdif(A,B,D,dim)
      integer dim 
      real A, B, D 
      dimension A(dim,dim), B(dim,dim), D(dim,dim)

      do m=1,dim
      do n=1,dim
         D(m,n)=A(m,n)-B(m,n)
      enddo
      enddo 
    
      return
      end

Ccccccccccccccc  VARIANCE  ccccccccccccccccccccccccccccccccccccccccccccc

        subroutine Variance(a1,S1,S2,Mean,Var,dim,history_index,index)
        integer dim,history_index,index
        real a1,S1,S2,Mean,Var,history_index_r,index_r 
        dimension a1(0:index)
        integer i
  
        S1=0
        S2=0
        Mean=0
        Var=0

        history_index_r = history_index*1.0
        index_r = index*1.0
     
        do i=0,index
          S1=S1+a1(i)
          S2=S2+(a1(i)**2)
        enddo

        Mean=S1/history_index_r
   
        Var=(1.0/index_r)*S2-(history_index_r/index_r)*(Mean**2)
   
        return
        end

Ccccccccccccccc  COVARIANCE  ccccccccccccccccccccccccccccccccccccccccccccc
 
      subroutine Covariance(a1,b1,S1,S2,S3,Covar,dim,history_index,
     +index)
        integer dim,history_index,index
        real a1,b1,S1,S2,S3,Covar,history_index_r,index_r 
        dimension a1(0:index), b1(0:index)
        integer i
   
        S1=0
        S2=0
        S3=0
        Covar=0
     
        history_index_r = history_index*1.0
        index_r = index*1.0

        do i=0,index
         S1=S1+a1(i)
         S2=S2+b1(i)
         S3=S3+a1(i)*b1(i)
        enddo

      Covar=(1.0/history_index_r)*S3 - 
     + (S1/history_index_r)*(S2/history_index_r)
  
        return
        end

Cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
Ccccc   This subroutine returns the i-line of the matrix E  cccccccccccc

        subroutine Line(E,dim,i,L1,index)
        integer i,n,dim,index
        real E,L1
        dimension E(dim,0:index), L1(0:index)
  
        do n=0,index
         L1(n)=E(i,n)
        enddo

        return
        end 
