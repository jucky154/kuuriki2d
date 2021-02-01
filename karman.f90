  program karman
  !$ use omp_lib
  implicit none
  integer::mdx=500,mdy=300
  real(8)::re,cfl,omegap,errorp,dx,dy,dt,res,resp,dp
  real(8)::ux,uy,vx,vy
  real(8)::t1,t2
  real(8)::x(-2:500,-2:300),y(-2:500,-2:300)
  real(8)::u(-1:500,-1:300),v(-1:500,-1:300),p(-1:500,-1:300)
  real(8)::rhs(-1:500,-1:300),urhs(-1:500,-1:300),vrhs(-1:500,-1:300)

  integer::maxitp,nlast,i1,i2,j1,j2,itr,i,j,nlp,mx,my,n,itrp
  integer::icent,jcent
  !call omp_set_numthread(4)
  t1= omp_get_wtime()
  !define
  re=70.0
  cfl=0.2

  omegap=1.0
  maxitp=100
  errorp=0.0001

  nlast=8000
  nlp=10

  mx=401
  i1=96
  i2=106

  my=201
  j1=96
  j2=106

  dx=1.0/float((i2-i1))
  dy=1.0/float((j2-j1))
  dt=cfl*dmin1(dx,dy)
  
  

  !make xygrid
  icent=(i1+i2)/2.0
  jcent=(j1+j2)/2.0
  do i=1,mx
     do j=1,my
        x(i,j)=dx*float(i-icent)
        y(i,j)=dy*float(j-jcent)
     end do
  end do
  
  !init condition
  do i=1,mx
     do j=1,my
        u(i,j)=1.0
        v(1,j)=0.0
        p(i,j)=0.0
     end do
  end do

  !bcforp
  i=1
  do j=1,my
     p(i,j)=0.0
  end do

  i=mx
  do j=1,my
     p(i,j)=0.0
  end do

  j=1
  do i=1,mx
     p(i,j)=0.0
  end do

  j=my
  do i=1,mx
     p(i,j)=0.0
  end do

  p(i1,j1)=p(i1-1,j1-1)
  p(i1,j2)=p(i1-1,j2+1)
  p(i2,j1)=p(i2+1,j1-1)
  p(i2,j2)=p(i2+1,j2+1)

  i=i1
  do j=j1+1,j2-1
     p(i,j)=p(i-1,j)
  end do
  
  i=i2
  do j=j1+1,j2-1
     p(i,j)=p(i+1,j)
  end do

  j=j1
  do i=i1+1,i2-1
     p(i,j)=p(i,j-1)
  end do

  j=j2
  do i=i1+1,i2-1
     p(i,j)=p(i,j+1)
  end do

  !bcforv
  i=1
  do j=1,my
     u(i,j)=1.0
     v(i,j)=0.0
     u(i-1,j)=1.0
     v(i-1,j)=0.0
  end do

  i=mx
  do j=1,my
     u(i,j)=2.0*u(i-1,j)-u(i-2,j)
     v(i,j)=2.0*v(i-1,j)-v(i-2,j)
     u(i+1,j)=2.0*u(i,j)-u(i-1,j)
     v(i+1,j)=2.0*v(i,j)-v(i-1,j)
  end do

  j=1
  do i=1,mx
     u(i,j)=2.0*u(i,j+1)-u(i,j+2)
     v(i,j)=2.0*v(i,j+1)-v(i,j+2)
     u(i,j-1)=2.0*u(i,j)-u(i,j+1)
     v(i,j-1)=2.0*v(i,j)-v(i,j+1)
  end do

  j=my
  do i=1,mx
     u(i,j)=2.0*u(i,j-1)-u(i,j-2)
     v(i,j)=2.0*v(i,j-1)-v(i,j-2)
     u(i,j+1)=2.0*u(i,j)-u(i,j-1)
     v(i,j+1)=2.0*v(i,j)-v(i,j-1)
  end do

  do i=i1,i2
     do j=j1,j2
        u(i,j)=0.0
        v(i,j)=0.0
     end do
  end do

  do n=1,nlast
     !poiseq
     !$omp parallel do private(i,j,ux,uy,vx,vy)
     do i=2,mx-1
        do j=2,my-1
           if(i .ge. i1 .and. i .le. i2 .and. j .ge. j1 .and. j .le. j2) go to 1000
           ux=(u(i+1,j)-u(i-1,j))/(2.0*dx)
           uy=(u(i,j+1)-u(i,j-1))/(2.0*dy)
           vx=(v(i+1,j)-v(i-1,j))/(2.0*dx)
           vy=(v(i,j+1)-v(i,j-1))/(2.0*dy)
           rhs(i,j)=(ux+vy)/dt -(ux*ux+2.0*uy*vy+vy*vy)
1000       continue
        end do
     end do
     !$omp end parallel do 

     do itr=1,maxitp
        res=0.0
        do i=2,mx-1
           do j=2,my-1
              if(i .ge. i1 .and. i .le. i2 .and. j .ge. j1 .and. j .le. j2) go to 2000
              dp=(p(i+1,j)+p(i-1,j))/(dx*dx)+(p(i,j+1)+p(i,j-1))/(dy*dy)-rhs(i,j)
              dp=dp/(2.0/(dx*dx)+2.0/(dy*dy))-p(i,j)
              res=res+dp*dp
              p(i,j)=p(i,j)+omegap*dp
2000          continue
           end do
        end do

        !bcforp
        i=1
        do j=1,my
           p(i,j)=0.0
        end do

        i=mx
        do j=1,my
           p(i,j)=0.0
        end do

        j=1
        do i=1,mx
           p(i,j)=0.0
        end do

        j=my
        do i=1,mx
           p(i,j)=0.0
        end do

        p(i1,j1)=p(i1-1,j1-1)
        p(i1,j2)=p(i1-1,j2+1)
        p(i2,j1)=p(i2+1,j1-1)
        p(i2,j2)=p(i2+1,j2+1)

        i=i1
        do j=j1+1,j2-1
           p(i,j)=p(i-1,j)
        end do
  
        i=i2
        do j=j1+1,j2-1
           p(i,j)=p(i+1,j)
        end do

        j=j1
        do i=i1+1,i2-1
           p(i,j)=p(i,j-1)
        end do

        j=j2
        do i=i1+1,i2-1
           p(i,j)=p(i,j+1)
        end do

 !       print*,'p',itr
        res=sqrt(res/float(mx*my))
        if(res .lt. errorp) go to 2999
     end do

2999 continue

     resp=res
     itrp=itr

     !bcforp
     i=1
     do j=1,my
        p(i,j)=0.0
     end do

     i=mx
     do j=1,my
        p(i,j)=0.0
     end do

     j=1
     do i=1,mx
        p(i,j)=0.0
     end do
     
     j=my
     do i=1,mx
        p(i,j)=0.0
     end do

     p(i1,j1)=p(i1-1,j1-1)
     p(i1,j2)=p(i1-1,j2+1)
     p(i2,j1)=p(i2+1,j1-1)
     p(i2,j2)=p(i2+1,j2+1)

     i=i1
     do j=j1+1,j2-1
        p(i,j)=p(i-1,j)
     end do
  
     i=i2
     do j=j1+1,j2-1
        p(i,j)=p(i+1,j)
     end do

     j=j1
     do i=i1+1,i2-1
        p(i,j)=p(i,j-1)
     end do

     j=j2
     do i=i1+1,i2-1
        p(i,j)=p(i,j+1)
     end do

     !veloeq
     do i=2,mx-1
        do j=2,my-1
           if(i .ge. i1 .and. i .le. i2 .and. j .ge. j1 .and. j .le. j2) go to 3000
           urhs(i,j)=-(p(i+1,j)-p(i-1,j))/(2.0*dx)
           vrhs(i,j)=-(p(i,j+1)-p(i,j-1))/(2.0*dy)
3000       continue
        end do
     end do

     !$omp parallel do private(i,j)
     do i=2,mx-1
        do j=2,my-1
           if(i .ge. i1 .and. i .le. i2 .and. j .ge. j1 .and. j .le. j2) go to 4000
           urhs(i,j)=urhs(i,j)+(u(i+1,j)-2.0*u(i,j)+u(i-1,j))/(re*dx*dx)+(u(i,j+1)-2.0*u(i,j)+u(i,j-1))/(re*dy*dy)
           vrhs(i,j)=vrhs(i,j)+(v(i+1,j)-2.0*v(i,j)+v(i-1,j))/(re*dx*dx)+(v(i,j+1)-2.0*v(i,j)+v(i,j-1))/(re*dy*dy)
4000       continue
        end do
     end do
     !$omp end parallel do
     
     do j=j1+1,j2-1
        u(i1+1,j)=2.0*u(i1,j)-u(i1-1,j)
        u(i2-1,j)=2.0*u(i2,j)-u(i2+1,j)
        v(i1+1,j)=2.0*v(i1,j)-v(i1-1,j)
        v(i2-1,j)=2.0*v(i2,j)-v(i2+1,j)
     end do

     !$omp parallel do private(i,j)
     do i=2,mx-1
        do j=2,my-1
           if(i .ge. i1 .and. i .le. i2 .and. j .ge. j1 .and. j .le. j2) go to 4001
           urhs(i,j)=urhs(i,j)-u(i,j)*(-u(i+2,j)+8.0*(u(i+1,j)-u(i-1,j))+u(i-2,j))/(12.0*dx)
           urhs(i,j)=urhs(i,j)-abs(u(i,j))*(u(i+2,j)-4.0*u(i+1,j)+6.0*u(i,j)-4.0*u(i-1,j)+u(i-2,j))/(4.0*dx)
           vrhs(i,j)=vrhs(i,j)-u(i,j)*(-v(i+2,j)+8.0*(v(i+1,j)-v(i-1,j))+v(i-2,j))/(12.0*dx)
           vrhs(i,j)=vrhs(i,j)-abs(u(i,j))*(v(i+2,j)-4.0*v(i+1,j)+6.0*v(i,j)-4.0*v(i-1,j)+v(i-2,j))/(4.0*dx)
4001       continue
        end do
     end do
     !$omp end parallel do
     
     do i=i1+1,i2-1
        u(i,j1+1)=2.0*u(i,j1)-u(i,j1-1)
        u(i,j2-1)=2.0*u(i,j2)-u(i,j2+1)
        v(i,j1+1)=2.0*v(i,j1)-v(i,j1-1)
        v(i,j2-1)=2.0*v(i,j2)-v(i,j2+1)
     end do

     !$omp parallel do private(i,j)
     do i=2,mx-1
        do j=2,my-1
           if(i .ge. i1 .and. i .le. i2 .and. j .ge. j1 .and. j .le. j2) go to 4002
           urhs(i,j)=urhs(i,j)-v(i,j)*(-u(i,j+2)+8.0*(u(i,j+1)-u(i,j-1))+u(i,j-2))/(12.0*dy)
           urhs(i,j)=urhs(i,j)-abs(v(i,j))*(u(i,j+2)-4.0*u(i,j+1)+6.0*u(i,j)-4.0*u(i,j-1)+u(i,j-2))/(4.0*dy)
           vrhs(i,j)=vrhs(i,j)-v(i,j)*(-v(i,j+2)+8.0*(v(i,j+1)-v(i,j-1))+v(i,j-2))/(12.0*dy)
           vrhs(i,j)=vrhs(i,j)-abs(v(i,j))*(v(i,j+2)-4.0*v(i,j+1)+6.0*v(i,j)-4.0*v(i,j-1)+v(i,j-2))/(4.0*dy)
4002       continue
        end do
     end do
     !$omp end parallel do 
     
     do i=2,mx-1
        do j=2,my-1
           if(i .ge. i1 .and. i .le. i2 .and. j .ge. j1 .and. j .le. j2) go to 5000
           u(i,j)=u(i,j)+dt*urhs(i,j)
           v(i,j)=v(i,j)+dt*vrhs(i,j)
5000       continue
        end do
     end do

     !bcforv
     i=1
     do j=1,my
        u(i,j)=1.0
        v(i,j)=0.0
        u(i-1,j)=1.0
        v(i-1,j)=0.0
     end do

     i=mx
     do j=1,my
        u(i,j)=2.0*u(i-1,j)-u(i-2,j)
        v(i,j)=2.0*v(i-1,j)-v(i-2,j)
        u(i+1,j)=2.0*u(i,j)-u(i-1,j)
        v(i+1,j)=2.0*v(i,j)-v(i-1,j)
     end do

     j=1
     do i=1,mx
        u(i,j)=2.0*u(i,j+1)-u(i,j+2)
        v(i,j)=2.0*v(i,j+1)-v(i,j+2)
        u(i,j-1)=2.0*u(i,j)-u(i,j+1)
        v(i,j-1)=2.0*v(i,j)-v(i,j+1)
     end do

     j=my
     do i=1,mx
        u(i,j)=2.0*u(i,j-1)-u(i,j-2)
        v(i,j)=2.0*v(i,j-1)-v(i,j-2)
        u(i,j+1)=2.0*u(i,j)-u(i,j-1)
        v(i,j+1)=2.0*v(i,j)-v(i,j-1)
     end do

     do i=i1,i2
        do j=j1,j2
           u(i,j)=0.0
           v(i,j)=0.0
        end do
     end do
  ! print*,'n',n
  end do
  
  do i=1,mx
     do j=1,my
        print*,x(i,j),y(i,j),p(i,j)
     end do
     print*,''
  end do
  t2= omp_get_wtime()
  print*,t2-t1

end program karman


              
  
  
  
