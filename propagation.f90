subroutine propagation_1D
	
use data_grid
use data_au

 implicit none   
 include "/usr/include/fftw3.f"	 

 integer I, J, K, FNAME1, FNAME2
 integer L, M, N
 integer*8 planF, planB

 double precision:: time		! Die Zeit
 double precision:: R			! Kernabstand
 double precision:: cpm			! Abschneideparameter d. cut-off Funktion
 double precision:: verw(N_states), terw(N_states), perw(N_states) ! Erwartungswerte
 double precision:: norm(N_states), evR(N_states)	! Population, Ortserwartungswert
 double precision:: in_norm(N_states), out_norm(N_states), outnorm_sum ! Hilfsgrößen um auslaufende Norm zu berechnen

 double precision, allocatable, dimension(:):: cof	! Abschneidefunktion
 double precision, allocatable, dimension(:,:):: pdens	! Impulsraumdichten

 complex*16:: tout(N_states,N_states)			! Diagonalisierte Kopplungsmatrix
 complex*16, allocatable, dimension(:,:):: psi_ges	! Gesamtwellenfunktion in el. Zustaenden
 complex*16, allocatable, dimension(:,:):: psi_out_imp      ! auslaufende Gesamtwellenfunktion
 complex*16, allocatable, dimension(:):: psi, kprop	! Hilfsfunktion, kinetischer Propagator
 
 open(77,file='cut_off.out',status='unknown')
 open(78,file='psi0.out',status='unknown')
 open(79,file='psi0_imp.out',status='unknown')	 	
! open(80,file='R_erw.out',status='unknown')
! open(81,file='pot_erw.out',status='unknown')
! open(82,file='kin_erw.out',status='unknown')
! open(83,file='imp_erw.out',status='unknown')
! open(84,file='norm.out',status='unknown')
 open(85,file='outnorm.out',status='unknown')
 open(86,file='psi_out_imp.out',status='unknown')

 allocate(psi(NR),kprop(NR),psi_ges(NR,N_states),cof(NR),psi_out_imp(NR,N_states))
 allocate(pdens(Nr,N_states))

! print*
! print*,'Tuning FFTW...'

 call dfftw_plan_dft_1d(planF, NR, psi, psi, FFTW_FORWARD,FFTW_MEASURE)
 call dfftw_plan_dft_1d(planB, NR, psi, psi, FFTW_BACKWARD,FFTW_MEASURE)
  
             
! print*,'Done.'
! print*	  


 FNAME1 = 100
 FNAME2 = 1000
 psi_ges = (0.d0,0.d0) 
 psi_out_imp = (0.d0,0.d0)
 in_norm = 0.d0
 out_norm = 0.d0

 do I = 1, NR
  R = R0 + (i-1)* dR  
   kprop(I) = exp(-im *dt * PR(I)**2/(4.d0*Mass))  ! Propagator, kin. 
   psi_ges(I,initial_state) = exp(- 0.5d0 * (kappa**2/5.545177232d0) * (R - R_init)**2) * exp(im * P_init * (R-R_init))   ! Anfangsverteilung		
 end do 


 cpm = 3.d0  !................ Definition einer cut_off-Funktion
 
 do j = 1, NR
   R = R0 + (j - 1) * dR
    if(R.lt.(Rend - cpm)) then
    cof(j) = 1.d0
    else
       cof(j) = cos(((R - Rend + cpm) / -cpm) * (0.5d0 * pi))
       cof(j) = cof(j)**2
!       cof(j) = sin(((R - Rend + cpm) / -cpm) * (0.5d0 * pi))**4
!       cof(j) = 1-cof(j)
    end if

 end do

     
 
 call integ(psi_ges, norm)	 !... Normierung
 print*
 print*,'norm =', sngl(norm) 
 	 
 psi_ges(:,initial_state) = psi_ges(:,initial_state) / sqrt(norm(initial_state))
 
 call integ(psi_ges, norm)
 print*,'norm =', sngl(norm) 


 do I = 1, NR      
  R = R0 + (I-1) * dR  
  write(78,*) sngl(R ),sngl(abs(psi_ges(I,initial_state))) ! Anfangszustand  
  write(77,*) sngl(R), sngl(cof(i))		! Abschneidefunktion
 end do

    do I = 1,NR
       psi(I) = psi_ges(I,initial_state)   ! Hilfsgroesse, Kopie in 1D Groesse
    end do
    call dfftw_execute(planF)	 ! psi(r) --> psi(p)
    
    psi = psi / sqrt(dble(NR))

    do I = 1, NR      
       write(79,*) sngl(PR(I)),sngl(abs(psi(I))**1) ! Anfangszustand  
    end do

 
!______________________________________________________________________
!                                                                     
!                   Propagation Loop                         
!______________________________________________________________________
                             

 print*
 print*,'1D propagation...'
 print*
	

timeloop: do K = 0, Nt

    time = K * dt

!    if(mod(K,200).eq.0) then
!      print*,'time:', sngl(time)	     
!    end if


      terw = 0.d0
      perw = 0.d0 
      evR = 0.d0
      verw = 0.d0
 
  
    
!----------------- Propagation via Split Operator 3rd Ordnung ---------
!
!   
! ............... kinetische Propagation, 1. Teil ...................


    do J = 1,N_states	
      do I = 1,NR
	 psi(I) = psi_ges(I,J)   ! Hilfsgroesse, Kopie in 1D Groesse
      end do 	
      call dfftw_execute(planF)	 ! psi(r) --> psi(p)

      psi = psi * kprop  	 ! kin. Propagation
      call dfftw_execute(planB)	 ! psi(p) --> psi(r)
      psi = psi / dble(NR)	 ! Normieren, wg. numerischer FT
	 
      do I = 1,NR
	 psi_ges(I,J) = psi(I)   ! Hilfsgroesse, zurueck kopieren
      end do      
    end do
 
 
! ............... potentielle Propagation ..........................
  
   
    do i = 1, NR				! Kopplung zwischen den el. Zustaenden,
       call pulse2(tout,i) 		        ! geht ein in die off-Diagonalen -- hier
       psi_ges(i,1:N_states) = matmul(tout(1:N_states,1:N_states),psi_ges(i,1:N_states))  ! wird diagonalisiert	
    end do


   do j = 1, N_states
    do i = 1, NR  
    psi_ges(i,j) = psi_ges(i,j) * exp(-im * dt * pot(i,j)) ! potentielle Propagation    
    end do
   end do
   
  
  
! ............... kinetische Propagation, 2. Teil ...................  


    do J = 1,N_states
      do I = 1, NR
	 psi(I) = psi_ges(I,J)  	! Hilfsgroesse
      end do 
      call dfftw_execute(planF)       	! psi(r) --> psi(p)
      psi = psi * kprop		      	! kin. Propagation
      psi = psi / sqrt(dble(NR))	! Normieren, wg. numerischer FT

       do i = 1, NR
         terw(j) = terw(j) + abs(psi(i))**2 * (PR(I)**2/(2.d0*Mass)) ! Erwartungswerte, kin. Energie
         perw(j) = perw(j) + abs(psi(i))**2 * PR(I)	! Impulserwartungswerte
         pdens(i,j) = abs(psi(i))**2			! Impulsraumdichte
       end do  

      call dfftw_execute(planB)         ! psi(p) --> psi(r)                  
      psi = psi / sqrt(dble(NR))	! Normieren, wg. numerischer FT

      do I = 1, NR
	 psi_ges(I,J) = psi(I)		! Hilfsgroesse, zurueck kopieren
      end do      
    end do     
 

!-------------------- Erwartungswerte und Output -------------------------------

 
!    do j = 1, N_states
!     do i = 1, NR     
!      R = R0 + (i-1) *dR   
!      evR(j) = evR(j) + abs(psi_ges(i,j))**2 * R  ! Erwartungswerte, Kernabstand
!      verw(j) = verw(j)+ abs(psi_ges(i,j))**2 * pot(i,j) ! Erwartungswerte, pot. Energie
!     end do   
!    end do
!
!    evR = evR * dR
!    verw = verw *dR
!    terw = terw *dpR
!    perw = perw *dpR
!  
! 
!    call integ(psi_ges, norm)			! Norm in den elektr. Zustaenden
!    
!
!    do j = 1, N_states
!	if (norm(j).ge.1.d-10) then
!	evR(j) = evR(j) / norm(j)
!        verw(j) = verw(j) / norm(j)
!        terw(j) = terw(j) / norm(j)
!        perw(j) = perw(j) / norm(j)
!	end if
!    end do
         
!    write(80,"(F15.8,F15.8,F15.8)") time, evR(:)
!    write(81,"(F15.8,F15.8,F15.8)") time, verw(:)
!    write(82,"(F15.8,F15.8,F15.8)") time, terw(:)
!    write(83,"(F15.8,F15.8,F15.8)") time, perw(:)
!    write(84,"(F15.8,F15.8,F15.8)") time, norm(:)	  

!---------------------------------
! perform timestep for psi_out_imp
!---------------------------------
    do J = 1,N_states

       do I = 1,NR
          psi_out_imp(I,J) = psi_out_imp(I,J) * kprop(I) * exp(-im * dt * pot(I,J)) * kprop(I)
       end do

    end do

!-------------------------------------
! add new outgoing part to psi_out_imp
!-------------------------------------

    do J = 1, N_states
      do I = 1,NR
	 psi(I) = psi_ges(I,J) * (1.d0 - cof(I))  ! Hilfsgroesse, Kopie in 1D Groesse
      end do 	
      call dfftw_execute(planF)	 ! psi(r) --> psi(p)

      psi = psi / sqrt(dble(NR)) ! Normieren, wg. numerischer FT
	 
      do I = 1,NR
	 psi_out_imp(I,J) = psi_out_imp(I,J) + psi(I)   ! Hilfsgroesse, zurueck kopieren
      end do      
    end do

!------------------------------------------------------
! print (outgoing) densities in real and momentum space
!------------------------------------------------------

 if(mod(K,50).eq.0) then

!  do I = 1, NR  		! Dichten in den elektronischen Zustaenden
!     R = R0 + (I-1) *dR 	           
!     if ((abs(psi_ges(I,1)**2).gt.0.001d0).AND.(FNAME2-1000.gt.33)) write(FNAME2,*) sngl(R ), sngl(abs(psi_ges(I,1)**2)+pot(int(evR(1)/dr),1))
!     if (abs(psi_ges(I,2)**2).gt.0.001d0) write(FNAME1,*) sngl(R ), sngl(abs(psi_ges(I,2)**2)+pot(int(evR(2)/dr),2))
!!     write(FNAME1,*) sngl(R ), sngl(abs(psi_ges(I,3)**2))
!  end do 
!
!  do I = nr/2+1, Nr 	! Impulsraumdichten in den elektr. Zustaenden
!!     write(FNAME2,*) sngl(pr(i)), sngl(abs(psi_out_imp(i,1))**2)
!     write(FNAME2,*) sngl(pr(i)), sngl(abs(psi_out_imp(i,2))**2)
!!     write(FNAME2,*) sngl(pr(i)), sngl(pdens(i,2))
!!     write(FNAME2,*) sngl(pr(i)), sngl(abs(psi_out_imp(i,3))**2)
!  end do
!  do I = 1, Nr/2
!!     write(FNAME2,*) sngl(pr(i)), sngl(abs(psi_out_imp(i,1))**2)  
!     write(FNAME2,*) sngl(pr(i)), sngl(abs(psi_out_imp(i,2))**2)
!!     write(FNAME2,*) sngl(pr(i)), sngl(pdens(i,2))
!!     write(FNAME2,*) sngl(pr(i)), sngl(abs(psi_out_imp(i,3))**2)
!  end do

  FNAME1 = FNAME1 + 1
  FNAME2 = FNAME2 + 1

 end if   

!-------------------------------------
! cut off the outgoing part of psi_ges
!-------------------------------------

 call integ(psi_ges,norm)

 do j = 1, N_states
    do i = 1, NR
       psi_ges(i,j) = psi_ges(i,j) * cof(i)	! Abschneidefunktion; die 
    end do
 end do

!-----------------------------------------------------------
! calculate and print parts of norm running out of comp. box
!-----------------------------------------------------------

 call integ(psi_ges,in_norm)

 out_norm(:) = out_norm(:) + (norm(:) - in_norm(:))

 write(85,"(F15.8,F15.8,F15.8)") time, out_norm(:)
 
!------------------------------------------------------
! check sum of out_norm and end program if norm_sum ~ 1
!------------------------------------------------------

 outnorm_sum = 0.d0
 do i = 1, N_states
    outnorm_sum = outnorm_sum + out_norm(i)
 end do
 
 if (outnorm_sum.gt.0.995d0) then
    print*, 'out_norm reached 1'
    exit
 end if

end do timeloop

  do I = nr/2+1, Nr 	! Impulsraumdichten in den elektr. Zustaenden
     write(86,"(F15.8,F15.8,F15.8)") pr(i), abs(psi_out_imp(i,:))**2
  end do
  do I = 1, Nr/2
     write(86,"(F15.8,F15.8,F15.8)") pr(i), abs(psi_out_imp(i,:))**2
  end do

 close(77, status='keep')	    
 close(78, status='keep')
 close(79, status='keep')
! close(80, status='keep')
! close(81, status='keep')
! close(82, status='keep')
! close(83, status='keep')
! close(84, status='keep')
 close(85, status='keep')
 close(86, status='keep')
 
                               
 call dfftw_destroy_plan(planF)
 call dfftw_destroy_plan(planB)
 	    
 deallocate(psi, kprop, psi_ges, cof, psi_out_imp, pdens)
 
return
end subroutine


!_________________________________________________________


subroutine integ(psi, norm)
      
use data_grid
 implicit none
 integer I, J
      
 double precision,intent(out):: norm(N_states)
 complex*16,intent(in):: psi(Nr,N_states)
      
 norm = 0.d0
 
 do J = 1, N_states ! Norm in den N_states Zustaenden
  do I = 1, NR		  
    norm(J)= norm(J) + abs(psi(I,J))**2     				  
  end do
 end do
 
   norm = norm * dR
     
     
return 
end subroutine
!______________________________________________________________


subroutine pulse2(tout,r_index)

use data_grid, only:N_states,dt,diab_coupling_matrix
use data_au, only:im

implicit none

 integer:: J,info,r_index
 double precision :: w(N_states,N_states), u(N_states,N_states), d(N_states)
 
 complex*16 :: b(N_states,N_states), z(N_states,N_states)
 complex*16, intent(out):: tout(N_states,N_states)
 
!  ........... Definition der Kopplungsmatrix

 u(:,:) = diab_coupling_matrix(r_index,:,:)

 call jacobi2(u,N_states,N_states,d,w,info)

 b= (0.d0,0.d0)
 
 do J = 1,N_states
    b(J,J) = cmplx(dcos(dt*d(J)),-dsin(dt*d(J)))
    !b(J,J) = exp (-im * dt * d(J))
 end do

 !jacobi
 !z = matmul(u,b)
 !tout = matmul(z,transpose(u))
 
 ! jacobi2
 z = matmul(w,b)
 tout = matmul(z,transpose(w))
  
 return
end subroutine pulse2


!!-------------------------


      subroutine jacobi (mat,dim,ewerte)

         implicit none

         double precision genau
         parameter (genau=1.d-15)

         integer Jmax,mmax
         parameter (Jmax=15,mmax=18)
         integer matdim
         parameter (matdim=3)

         double precision mat(matdim,matdim)
         integer dim
         double precision ewerte(matdim)

         double precision s(matdim,matdim)
         integer ca,cb,p,q
         double precision c1,c2,t1,t2,t3,v1,v2,v3
         double precision tmp,l,n,t,m1,w,m
         logical flag

	s= 0.d0
         !!!!call fillmt(s,dim,dim,0.d0,matdim,matdim)

         do 1 ca=1,dim,1
            s(ca,ca)=1.d0
1           continue

         l=0.d0
         do 2 ca=2,dim,1
            do 12 cb=1,dim,1
               tmp=mat(ca,cb)
               l=l+2.d0*tmp*tmp
12          continue
2        continue

         n=dsqrt(l)
         m=genau*n/dim
         t=n

3        t=t/dim
4           do 6 q=2,dim,1
               do 16 p=1,q-1,1
                  flag=.false.
                  if (dabs(mat(p,q)).gt.t) then
                     flag=.true.
                     v1=mat(p,p)
                     v2=mat(p,q)
                     v3=mat(q,q)
                     m1=(v1-v3)/2.d0
                     if (m1.eq.0.d0) then
                           w=-1.d0
                        else
                           if (m1.gt.0.d0) then
                                 w=-v2/(dsqrt(v2*v2+m1*m1))
                              else
                                 w=v2/(dsqrt(v2*v2+m1*m1))
                              endif
                        endif

                     t1=w/dsqrt(2.d0*(1+dsqrt(1.d0-w/2.d0)))
                     t2=t1*t1
                     c1=dsqrt(1.d0-t2)
                     c2=c1*c1
                     t3=t1*c1

                     do 7 ca=1,dim,1
                        l=mat(ca,p)*c1-mat(ca,q)*t1
                        mat(ca,q)=mat(ca,p)*t1+mat(ca,q)*c1
                        mat(ca,p)=l
                        l=s(ca,p)*c1-s(ca,q)*t1
                        s(ca,q)=s(ca,p)*t1+s(ca,q)*c1
                        s(ca,p)=l
7                       continue
                     do 8 ca=1,dim,1
                        mat(p,ca)=mat(ca,p)
                        mat(q,ca)=mat(ca,q)
8                       continue
                     mat(p,p)=v1*c2+v3*t2-2*v2*t3
                     mat(q,q)=v1*t2+v3*c2+2*v2*t3
                     tmp=(v1-v3)*t3+v2*(c2-t2)
                     mat(p,q)=tmp
                     mat(q,p)=tmp
                     end if
16                continue
6              continue
               if (flag) go to 4
            if (m.lt.t) go to 3
            ewerte=0.d0
         !!!call fillvc(ewerte,dim,0.d0)
         do 9 ca=1,dim,1
            ewerte(ca)=mat(ca,ca)
9           continue
         do 10 ca=1,dim,1
            do 11 cb=1,dim,1
               mat(ca,cb)=s(ca,cb)
11          continue
10       continue

         return
         end


!______________________________________________________

      SUBROUTINE jacobi2(a,n,np,d,v,nrot)
      INTEGER n,np,nrot,NMAX
      DOUBLE PRECISION a(np,np),d(np),v(np,np)
      PARAMETER (NMAX=1000)
      INTEGER i,ip,iq,j
      DOUBLE PRECISION c,g,h,s,sm,t,tau,theta,tresh,b(NMAX),z(NMAX)
      do 12 ip=1,n
        do 11 iq=1,n
          v(ip,iq)=0.d0
11      continue
        v(ip,ip)=1.d0
12    continue
      do 13 ip=1,n
        b(ip)=a(ip,ip)
        d(ip)=b(ip)
        z(ip)=0.d0
13    continue
      nrot=0
      do 24 i=1,50
        sm=0.d0
        do 15 ip=1,n-1
          do 14 iq=ip+1,n
            sm=sm+abs(a(ip,iq))
14        continue
15      continue
        if(sm.eq.0.d0)return
        if(i.lt.4)then
          tresh=0.2d0*sm/n**2
        else
          tresh=0.d0
        endif
        do 22 ip=1,n-1
          do 21 iq=ip+1,n
            g=100.d0*abs(a(ip,iq))
            if((i.gt.4).and.(abs(d(ip))+g.eq.abs(d(ip))).and.(abs(d(iq))+g.eq.abs(d(iq))))then
              a(ip,iq)=0.d0
            else if(abs(a(ip,iq)).gt.tresh)then
              h=d(iq)-d(ip)
              if(abs(h)+g.eq.abs(h))then
                t=a(ip,iq)/h
              else
                theta=0.5d0*h/a(ip,iq)
                t=1.d0/(abs(theta)+sqrt(1.d0+theta**2))
                if(theta.lt.0.d0)t=-t
              endif
              c=1.d0/sqrt(1+t**2)
              s=t*c
              tau=s/(1.d0+c)
              h=t*a(ip,iq)
              z(ip)=z(ip)-h
              z(iq)=z(iq)+h
              d(ip)=d(ip)-h
              d(iq)=d(iq)+h
              a(ip,iq)=0.d0
              do 16 j=1,ip-1
                g=a(j,ip)
                h=a(j,iq)
                a(j,ip)=g-s*(h+g*tau)
                a(j,iq)=h+s*(g-h*tau)
16            continue
              do 17 j=ip+1,iq-1
                g=a(ip,j)
                h=a(j,iq)
                a(ip,j)=g-s*(h+g*tau)
                a(j,iq)=h+s*(g-h*tau)
17            continue
              do 18 j=iq+1,n
                g=a(ip,j)
                h=a(iq,j)
                a(ip,j)=g-s*(h+g*tau)
                a(iq,j)=h+s*(g-h*tau)
18            continue
              do 19 j=1,n
                g=v(j,ip)
                h=v(j,iq)
                v(j,ip)=g-s*(h+g*tau)
                v(j,iq)=h+s*(g-h*tau)
19            continue
              nrot=nrot+1
            endif
21        continue
22      continue
        do 23 ip=1,n
          b(ip)=b(ip)+z(ip)
          d(ip)=b(ip)
          z(ip)=0.d0
23      continue
24    continue
      pause 'too many iterations in jacobi'
      return
      END

