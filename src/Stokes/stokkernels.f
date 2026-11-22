c
c     This file contains the Stokes direct kernel evaluators
c
c**********************************************************************
c      
c     We take the following conventions for the Stokes kernels
c
c     For a source y and target x, let r_i = x_i-y_i
c     and let r = sqrt(r_1^2 + r_2^2 + r_3^2)
c
c     The Stokeslet, G_{ij}, and its associated pressure tensor, P_j,
c     are
c
c     G_{ij}(x,y) = 1/(8\pi)*(r_i r_j)/r^3 + 1/(8\pi)*delta_{ij}/r
c     P_j(x,y) = 1/(4\pi) * r_j/r^3
c
c     The (Type I) stresslet, T_{ijk}, and its associated pressure
c     tensor, PI_{jk}, are
c     
c     T_{ijk}(x,y) = 3/(4\pi) r_i r_j r_k/ r^5
c     PI_{jk} = -1/(2\pi) delta_{jk} + 3/(2\pi) r_j r_k/r^5      
c
c     The rotlet, R_{ijk}, and its associated pressure tensor, Q_{jk}, are
c
c     R_{ijk}(x,y) = -\delta_{ik} r_j/(4\pi r^3) + \delta_{ij} r_k/(4\pi r^3)
c     Q_{jk} = 0;
c
c     The doublet, D_{ijk}, and its associated pressure tensor, L_{jk}, are
c
c     D_{ijk}(x,y) = -\delta_{jk} r_i/(4\pi r^3) - \delta_{ik} r_j/(4\\pi r^3) + \delta_{ij} r_k/(4\pi r^3) + 3 r_i r_j r_k/ (4\pi r^5)
c     L_{jk} = -1/(2\pi) \delta_{jk} + 3 r_j r_k/(2\pi r^5)


      subroutine st3ddirectstokg(nd,sources,stoklet,ns,targ,nt,
     1     ifppreg,pot,pre,grad,thresh)
cf2py  intent(in) nd
cf2py  intent(in) sources
cf2py  intent(in) stoklet
cf2py  intent(in) ns
cf2py  intent(in) targ
cf2py  intent(in) nt
cf2py  intent(in) ifppreg
cf2py  intent(in) thresh
cf2py  intent(out) pot,pre,grad
c
c     This subroutine evaluates the potential and gradient due
c     to a collection of Stokeslet and stresslet sources and adds
c     to existing quantities (see definitions at top of file).
c
c       pot(x) = pot(x) + sum_m G_{ij}(x,y^{(m)}) sigma^{(m)}_j
c
c       pre(x) = sum_m P_j(x,y^m) sigma^{(m)}_j
c
c       grad(x) = Grad[sum_m G_{ij}(x,y^m) sigma^{(m)}_j
c      
c     where sigma^{(m)} is the Stokeslet charge, mu^{(m)} is the
c     stresslet charge, and nu^{(m)} is the stresslet orientation
c     (note that each of these is a 3 vector per source point y^{(m)}).
c     For x a source point, the self-interaction in the sum is omitted. 
c
c
c-----------------------------------------------------------------------
      
c     INPUT:
c     
c     nd in: integer *8
c        number of densities
c     
c     nsource in: integer *8  
c        number of sources
c
c     source  in: double precision (3,nsource)
c        source(k,j) is the kth component of the jth
c        source location
c
c     stoklet in: double precision (nd,3,nsource) 
c        Stokeslet charge strengths (sigma vectors above)
c
c     ntarg   in: integer *8  
c        number of targs 
c
c     targ    in: double precision (3,ntarg)
c        targ(k,j) is the kth component of the jth
c        targ location
c     
c     thresh in: double precision
c        threshold for updating potential,
c        potential at target won't be updated if
c        |t - s| <= thresh, where t is the target
c        location and, and s is the source location 
c
c     ifppreg in: integer *8
c        flag for evaluating potential, pressure, and gradient
c           ifppreg = 0, nothing is computed
c           ifppreg = 1, only potential
c           ifppreg = 2, potential and pressure
c           ifppreg = 3, potential, pressure, and gradient
c
c-----------------------------------------------------------------------
c
c   OUTPUT:
c
c     pot out: double precision(nd,3,ntarg) 
c        velocity at the targets
c      
c     pre out: double precision(nd,ntarg)
c        pressure at the targets
c      
c     grad out: double precision(nd,3,3,ntarg) 
c        gradient of velocity at the targets
c        gradtarg(l,i,j,k) is the ith component of the
c        gradient of the jth component of the velocity
c        for the lth density at the kth target
c
c     TODO: implement other stresslet options
c------------------------------------------------------------------
c$    use omp_lib, only: omp_in_parallel
      implicit none

      integer *8 nd, ns, nt, ifppreg
      real *8 sources(3,ns),targ(3,nt)
      real *8 stoklet(nd,3,ns)
      real *8 pot(nd,3,nt),pre(nd,nt),grad(nd,3,3,nt)
      real *8 thresh
      
c     local

      real *8 zdiff(3), d1, d2, d3, tempx1, tempx2, tempx3
      real *8 pl, pv, dmu(3), dnu(3), temp, r, r2, r3, r5
      real *8 dmunu      
      real *8 threshsq
      real *8, parameter :: inv4pi = 7.957747154594766788444188168626d-2

      integer *8 i, j, idim, l

      threshsq = thresh**2

      if (ifppreg .lt. 1) return

c     stokeslet contribution
      
c     Each target owns its outputs; retain source summation order.
c     Avoid nested teams when called inside an FMM parallel loop.
c$omp parallel do default(shared) schedule(static)
c$omp& if (.not.omp_in_parallel() .and. nt*ns .ge. 4096)
c$omp& private(i,j,idim,zdiff,r2,r,r3,pl,d1,d2,d3)
      do i = 1,nt
         do j = 1,ns
            zdiff(1) = targ(1,i)-sources(1,j)
            zdiff(2) = targ(2,i)-sources(2,j)
            zdiff(3) = targ(3,i)-sources(3,j)

            r2 = zdiff(1)**2 + zdiff(2)**2 + zdiff(3)**2
            if (r2 .lt. threshsq) goto 10

            r = sqrt(r2)
            r3 = r*r2
            
            do idim = 1,nd

               pot(idim,1,i) = pot(idim,1,i) +
     1                         stoklet(idim,1,j)/(2*r)*inv4pi
               pot(idim,2,i) = pot(idim,2,i) +
     1                         stoklet(idim,2,j)/(2*r)*inv4pi
               pot(idim,3,i) = pot(idim,3,i) +
     1                         stoklet(idim,3,j)/(2*r)*inv4pi
               
               pl = (zdiff(1)*stoklet(idim,1,j) +
     1              zdiff(2)*stoklet(idim,2,j) +
     2              zdiff(3)*stoklet(idim,3,j))/(r3*2)
               
               pl = pl * inv4pi
               pot(idim,1,i) = pot(idim,1,i) + zdiff(1)*pl
               pot(idim,2,i) = pot(idim,2,i) + zdiff(2)*pl
               pot(idim,3,i) = pot(idim,3,i) + zdiff(3)*pl

               if (ifppreg .ge. 2) pre(idim,i) = pre(idim,i) + pl*2

               if (ifppreg .lt. 3) cycle

               grad(idim,1,1,i) = grad(idim,1,1,i) + pl
               grad(idim,2,2,i) = grad(idim,2,2,i) + pl
               grad(idim,3,3,i) = grad(idim,3,3,i) + pl               

               d1 = stoklet(idim,1,j)/(r3*2)*inv4pi -
     1              zdiff(1)*pl*3.0d0/r2
               d2 = stoklet(idim,2,j)/(r3*2)*inv4pi -
     1              zdiff(2)*pl*3.0d0/r2
               d3 = stoklet(idim,3,j)/(r3*2)*inv4pi -
     1              zdiff(3)*pl*3.0d0/r2
               
               grad(idim,1,1,i) = grad(idim,1,1,i) + d1*zdiff(1)
               grad(idim,2,1,i) = grad(idim,2,1,i) + d2*zdiff(1)
               grad(idim,3,1,i) = grad(idim,3,1,i) + d3*zdiff(1)
               grad(idim,1,2,i) = grad(idim,1,2,i) + d1*zdiff(2)
               grad(idim,2,2,i) = grad(idim,2,2,i) + d2*zdiff(2)
               grad(idim,3,2,i) = grad(idim,3,2,i) + d3*zdiff(2)
               grad(idim,1,3,i) = grad(idim,1,3,i) + d1*zdiff(3)
               grad(idim,2,3,i) = grad(idim,2,3,i) + d2*zdiff(3)
               grad(idim,3,3,i) = grad(idim,3,3,i) + d3*zdiff(3)

               d1 = -stoklet(idim,1,j)/(2*r3)*inv4pi
               d2 = -stoklet(idim,2,j)/(2*r3)*inv4pi
               d3 = -stoklet(idim,3,j)/(2*r3)*inv4pi               
               
               grad(idim,1,1,i) = grad(idim,1,1,i) + zdiff(1)*d1
               grad(idim,2,1,i) = grad(idim,2,1,i) + zdiff(2)*d1
               grad(idim,3,1,i) = grad(idim,3,1,i) + zdiff(3)*d1
               grad(idim,1,2,i) = grad(idim,1,2,i) + zdiff(1)*d2
               grad(idim,2,2,i) = grad(idim,2,2,i) + zdiff(2)*d2
               grad(idim,3,2,i) = grad(idim,3,2,i) + zdiff(3)*d2
               grad(idim,1,3,i) = grad(idim,1,3,i) + zdiff(1)*d3
               grad(idim,2,3,i) = grad(idim,2,3,i) + zdiff(2)*d3
               grad(idim,3,3,i) = grad(idim,3,3,i) + zdiff(3)*d3
            enddo
 10         continue
         enddo
      enddo
c$omp end parallel do

      
      return
      end

      subroutine st3ddirectstokstrsg(nd,sources,stoklet,istress,
     1     strslet,strsvec,ns,targ,nt,ifppreg,pot,pre,grad,thresh)
cf2py  intent(in) nd
cf2py  intent(in) sources
cf2py  intent(in) stoklet,istress
cf2py  intent(in) strslet,strsvec
cf2py  intent(in) ns
cf2py  intent(in) targ
cf2py  intent(in) nt
cf2py  intent(in) ifppreg
cf2py  intent(in) thresh
cf2py  intent(out) pot,pre,grad
c
c     This subroutine evaluates the potential and gradient due
c     to a collection of Stokeslet and stresslet sources and adds
c     to existing quantities (see definitions at top of file).
c
c       pot(x) = pot(x) + sum_m G_{ij}(x,y^{(m)}) sigma^{(m)}_j
c                + sum_m T_{ijk}(x,y^{(m)}) mu^{(m)}_j nu^{(m)}_k
c
c       pre(x) = sum_m P_j(x,y^m) sigma^{(m)}_j
c          + sum_m T_{ijk}(x,y^{(m)}) PI_{jk} mu^{(m)}_j nu^{(m)}_k
c
c       grad(x) = Grad[sum_m G_{ij}(x,y^m) sigma^{(m)}_j
c                + sum_m T_{ijk}(x,y^{(m)}) mu^{(m)}_j nu^{(m)}_k]
c      
c     where sigma^{(m)} is the Stokeslet charge, mu^{(m)} is the
c     stresslet charge, and nu^{(m)} is the stresslet orientation
c     (note that each of these is a 3 vector per source point y^{(m)}).
c     For x a source point, the self-interaction in the sum is omitted. 
c
c
c-----------------------------------------------------------------------
      
c     INPUT:
c     
c     nd in: integer *8
c        number of densities
c     
c     nsource in: integer *8  
c        number of sources
c
c     source  in: double precision (3,nsource)
c        source(k,j) is the kth component of the jth
c        source location
c
c     stoklet in: double precision (nd,3,nsource) 
c        Stokeslet charge strengths (sigma vectors above)
c
c     istress in: integer *8
c        stresslet computation flag
c           istress = 1   =>  include standard stresslet
c                               (type I)
c     
c           NOT YET IMPLEMENTED
c      
c           ifstress = 2   =>  include symmetric stresslet
c                                   (type II)
c           ifstress = 3   =>  include rotlet
c           ifstress = 4   =>  include Stokes doublet
c                      otherwise do not include
c
c     strslet  in: double precision (nd,3,nsource) 
c        stresslet strengths (mu vectors above)
c
c     strsvec  in: double precision (nd,3,nsource)   
c        stresslet orientations (nu vectors above)
c      
c     ntarg   in: integer *8  
c        number of targs 
c
c     targ    in: double precision (3,ntarg)
c        targ(k,j) is the kth component of the jth
c        targ location
c     
c     thresh in: double precision
c        threshold for updating potential,
c        potential at target won't be updated if
c        |t - s| <= thresh, where t is the target
c        location and, and s is the source location 
c
c     ifppreg in: integer *8
c        flag for evaluating potential, pressure, and gradient
c           ifppreg = 0, nothing is computed
c           ifppreg = 1, only potential
c           ifppreg = 2, potential and pressure
c           ifppreg = 3, potential, pressure, and gradient
c
c-----------------------------------------------------------------------
c
c   OUTPUT:
c
c     pot out: double precision(nd,3,ntarg) 
c        velocity at the targets
c      
c     pre out: double precision(nd,ntarg)
c        pressure at the targets
c      
c     grad out: double precision(nd,3,3,ntarg) 
c        gradient of velocity at the targets
c        gradtarg(l,i,j,k) is the ith component of the
c        gradient of the jth component of the velocity
c        for the lth density at the kth target
c
c     TODO: implement other stresslet options
c------------------------------------------------------------------
c$    use omp_lib, only: omp_in_parallel
      implicit none

      integer *8 nd, ns, nt, istress, ifppreg
      real *8 sources(3,ns),targ(3,nt),strslet(nd,3,ns)
      real *8 strsvec(nd,3,ns),stoklet(nd,3,ns)
      real *8 pot(nd,3,nt),pre(nd,nt),grad(nd,3,3,nt)
      real *8 thresh
      
c     local

      real *8 zdiff(3), d1, d2, d3, tempx1, tempx2, tempx3
      real *8 pl, pv, dmu(3), dnu(3), temp, r, r2, r3, r5
      real *8 dmunu      
      real *8 threshsq
      real *8, parameter :: inv4pi = 7.957747154594766788444188168626d-2

      integer *8 i, j, idim, l


      call st3ddirectstokg(nd,sources,stoklet,ns,targ,nt,
     1     ifppreg,pot,pre,grad,thresh)

      
      threshsq = thresh**2

      if (ifppreg .lt. 1) return

      if (istress .ne. 1) goto 100

c     type I stresslet
      
c     Each target owns its outputs; retain source summation order.
c     Avoid nested teams when called inside an FMM parallel loop.
c$omp parallel do default(shared) schedule(static)
c$omp& if (.not.omp_in_parallel() .and. nt*ns .ge. 4096)
c$omp& private(i,j,idim,zdiff,r2,r,r3,r5,dmu,dnu,pl,pv,temp,dmunu)
c$omp& private(tempx1,tempx2,tempx3)
      do i = 1,nt
         do j = 1,ns
            zdiff(1) = targ(1,i)-sources(1,j)
            zdiff(2) = targ(2,i)-sources(2,j)
            zdiff(3) = targ(3,i)-sources(3,j)

            r2 = zdiff(1)**2 + zdiff(2)**2 + zdiff(3)**2
            if (r2 .lt. threshsq) goto 20

            r = sqrt(r2)
            r3 = r*r2
            r5 = r3*r2
            
            do idim = 1,nd

               dmu(1) = strslet(idim,1,j)
               dmu(2) = strslet(idim,2,j)
               dmu(3) = strslet(idim,3,j)
               dnu(1) = strsvec(idim,1,j)
               dnu(2) = strsvec(idim,2,j)
               dnu(3) = strsvec(idim,3,j)

               pl = zdiff(1)*dmu(1) + zdiff(2)*dmu(2) + zdiff(3)*dmu(3)
               pv = zdiff(1)*dnu(1) + zdiff(2)*dnu(2) + zdiff(3)*dnu(3)

               temp = 3.0d0*pl*pv/r5
               temp = temp * inv4pi
               
               pot(idim,1,i) = pot(idim,1,i) + zdiff(1)*temp
               pot(idim,2,i) = pot(idim,2,i) + zdiff(2)*temp
               pot(idim,3,i) = pot(idim,3,i) + zdiff(3)*temp

               if (ifppreg .ge. 2) then
                  dmunu = dmu(1)*dnu(1) + dmu(2)*dnu(2) + dmu(3)*dnu(3)
                  dmunu = dmunu * inv4pi
                  pre(idim,i) = pre(idim,i) - 2.0d0*dmunu/r3 +
     1                          2.0d0*temp
               endif

               if (ifppreg .lt. 3) cycle

               tempx1 = 3.0d0*(dmu(1)*pv + dnu(1)*pl -
     1              5.0d0*zdiff(1)*pl*pv/r2)/ r5
               tempx2 = 3.0d0*(dmu(2)*pv + dnu(2)*pl -
     1              5.0d0*zdiff(2)*pl*pv/r2)/ r5
               tempx3 = 3.0d0*(dmu(3)*pv + dnu(3)*pl -
     1              5.0d0*zdiff(3)*pl*pv/r2)/ r5

               tempx1 = tempx1 * inv4pi
               tempx2 = tempx2 * inv4pi
               tempx3 = tempx3 * inv4pi

               grad(idim,1,1,i) = grad(idim,1,1,i) + temp
               grad(idim,1,1,i) = grad(idim,1,1,i) + zdiff(1)*tempx1
               grad(idim,2,1,i) = grad(idim,2,1,i) + zdiff(1)*tempx2
               grad(idim,3,1,i) = grad(idim,3,1,i) + zdiff(1)*tempx3
               
               grad(idim,2,2,i) = grad(idim,2,2,i) + temp
               grad(idim,1,2,i) = grad(idim,1,2,i) + zdiff(2)*tempx1
               grad(idim,2,2,i) = grad(idim,2,2,i) + zdiff(2)*tempx2
               grad(idim,3,2,i) = grad(idim,3,2,i) + zdiff(2)*tempx3
               
               grad(idim,3,3,i) = grad(idim,3,3,i) + temp
               grad(idim,1,3,i) = grad(idim,1,3,i) + zdiff(3)*tempx1
               grad(idim,2,3,i) = grad(idim,2,3,i) + zdiff(3)*tempx2
               grad(idim,3,3,i) = grad(idim,3,3,i) + zdiff(3)*tempx3
            enddo
            
 20         continue
         enddo
      enddo
c$omp end parallel do
      
      return
      
 100  continue
      
      return
      end
      


      subroutine st3ddirectstokstrsrotdoubg(nd,sources,stoklet,istress,
     1     strslet,strsvec,irotlet,rotlet,rotvec,
     2     idoublet,doublet,doubvec,
     3     ns,targ,nt,ifppreg,pot,pre,grad,thresh)
cf2py  intent(in) nd
cf2py  intent(in) sources
cf2py  intent(in) stoklet,istress
cf2py  intent(in) strslet,strsvec
cf2py  intent(in) irotlet,rotlet,rotvec
cf2py  intent(in) idoublet,doublet,doubvec
cf2py  intent(in) ns
cf2py  intent(in) targ
cf2py  intent(in) nt
cf2py  intent(in) ifppreg
cf2py  intent(in) thresh
cf2py  intent(out) pot,pre,grad
c
c     This subroutine evaluates the potential and gradient due
c     to a collection of Stokeslet and stresslet sources and adds
c     to existing quantities (see definitions at top of file).
c
c       pot(x) = pot(x) + sum_m G_{ij}(x,y^{(m)}) sigma^{(m)}_j
c                + sum_m T_{ijk}(x,y^{(m)}) mu^{(m)}_j nu^{(m)}_k
c                + sum_m R_{ijk}(x,y^{(m)}) rlet^{(m)}_j rvec^{(m)}_k
c                + sum_m D_{ijk}(x,y^{(m)}) dlet^{(m)}_j dvec^{(m)}_k
c
c       pre(x) = sum_m P_j(x,y^m) sigma^{(m)}_j
c                + sum_m PI_{jk} mu^{(m)}_j nu^{(m)}_k
c                + sum_m L_{jk}(x,y^{(m)}) dlet^{(m)}_j dvec^{(m)}_k
c
c       grad(x) = Grad[sum_m G_{ij}(x,y^m) sigma^{(m)}_j
c                + sum_m T_{ijk}(x,y^{(m)}) mu^{(m)}_j nu^{(m)}_k]
c                + sum_m R_{ijk}(x,y^{(m)}) rlet^{(m)}_j rvec^{(m)}_k
c                + sum_m D_{ijk}(x,y^{(m)}) dlet^{(m)}_j dvec^{(m)}_k
c
c     where sigma^{(m)} is the Stokeslet charge, mu^{(m)} is the
c     stresslet charge, nu^{(m)} is the stresslet orientation
c     rlet^{(m)} is the rotlet strength, rvec^{(m)} is the rotlet
c     dlet^{(m)} is the doublet strength, and dvec^{(m)} is the doublet
c     (note that each of these is a 3 vector per source point y^{(m)}).
c     For x a source point, the self-interaction in the sum is omitted.
c
c
c-----------------------------------------------------------------------
c     INPUT:
c
c     nd in: integer *8
c        number of densities
c
c     nsource in: integer *8
c        number of sources
c
c     source  in: double precision (3,nsource)
c        source(k,j) is the kth component of the jth
c        source location
c
c     stoklet in: double precision (nd,3,nsource)
c        Stokeslet charge strengths (sigma vectors above)
c
c     istress in: integer *8
c        stresslet computation flag
c           istress = 1   =>  include standard stresslet
c                               (type I)
c
c           NOT YET IMPLEMENTED
c
c           istress = 2   =>  include symmetric stresslet
c                                   (type II)
c           istress = 3   =>  include rotlet
c           istress = 4   =>  include Stokes doublet
c                      otherwise do not include
c
c     strslet  in: double precision (nd,3,nsource)
c        stresslet strengths (mu vectors above)
c
c     strsvec  in: double precision (nd,3,nsource)
c        stresslet orientations (nu vectors above)
c
c     irotlet in: integer *8
c
c     rotlet  in: double precision (nd,3,nsource)
c        rotlet strengths
c
c     rotvec  in: double precision (nd,3,nsource)
c        rotlet orientations
c
c     idoublet in: integer *8
c
c     doublet  in: double precision (nd,3,nsource)
c        doublet strengths
c
c     doubvec  in: double precision (nd,3,nsource)
c        doublet orientations
c
c     ntarg   in: integer *8
c        number of targs
c
c     targ    in: double precision (3,ntarg)
c        targ(k,j) is the kth component of the jth
c        targ location
c
c     thresh in: double precision
c        threshold for updating potential,
c        potential at target won't be updated if
c        |t - s| <= thresh, where t is the target
c        location and, and s is the source location
c
c     ifppreg in: integer *8
c        flag for evaluating potential, pressure, and gradient
c           ifppreg = 0, nothing is computed
c           ifppreg = 1, only potential
c           ifppreg = 2, potential and pressure
c           ifppreg = 3, potential, pressure, and gradient
c
c-----------------------------------------------------------------------
c
c   OUTPUT:
c
c     pot out: double precision(nd,3,ntarg)
c        velocity at the targets
c
c     pre out: double precision(nd,ntarg)
c        pressure at the targets
c
c     grad out: double precision(nd,3,3,ntarg)
c        gradient of velocity at the targets
c        gradtarg(l,i,j,k) is the ith component of the
c        gradient of the jth component of the velocity
c        for the lth density at the kth target
c
c     TODO: implement other stresslet options
c------------------------------------------------------------------
c$    use omp_lib, only: omp_in_parallel
      implicit none

      integer *8 nd, ns, nt, istress, irotlet, idoublet, ifppreg
      real *8 sources(3,ns),targ(3,nt),strslet(nd,3,ns)
      real *8 rotlet(nd,3,ns),rotvec(nd,3,ns)
      real *8 doublet(nd,3,ns),doubvec(nd,3,ns)
      real *8 strsvec(nd,3,ns),stoklet(nd,3,ns)
      real *8 pot(nd,3,nt),pre(nd,nt),grad(nd,3,3,nt)
      real *8 thresh, threshsq
c     local
      real *8 zdiff(3), tempx1, tempx2, tempx3
      real *8 pl, pv, dmu(3), dnu(3), temp, r, r2, r3, r5
      real *8 r3inv
      real *8 rpl, rpv, dpl, dpv, dlv
      real *8 rlet(3),rvec(3),dlet(3),dvec(3)
      real *8 dmunu
      real *8, parameter :: inv4pi = 7.957747154594766788444188168626d-2

      integer *8 i, j, idim, l

      if (istress .ne. 1) then
         call st3ddirectstokg(nd,sources,stoklet,ns,targ,nt,
     1        ifppreg,pot,pre,grad,thresh)
      else
         call st3ddirectstokstrsg(nd,sources,stoklet,istress,
     1        strslet,strsvec,ns,targ,nt,ifppreg,pot,pre,grad,thresh)
      endif

      threshsq = thresh**2

      if (ifppreg .lt. 1) return

      if (irotlet .ne. 1 .and. idoublet .ne. 1) goto 100

c     rotlet and doublet
c     Each target owns its outputs; retain source summation order.
c     Avoid nested teams when called inside an FMM parallel loop.
c$omp parallel do default(shared) schedule(static)
c$omp& if (.not.omp_in_parallel() .and. nt*ns .ge. 4096)
c$omp& private(i,j,idim,zdiff,r2,r,r3,r5,r3inv,rlet,rvec,rpl,rpv)
c$omp& private(dlet,dvec,dpl,dpv,dlv,temp,tempx1,tempx2,tempx3)
      do i = 1,nt
         do j = 1,ns
            zdiff(1) = targ(1,i)-sources(1,j)
            zdiff(2) = targ(2,i)-sources(2,j)
            zdiff(3) = targ(3,i)-sources(3,j)

            r2 = zdiff(1)**2 + zdiff(2)**2 + zdiff(3)**2
            if (r2 .lt. threshsq) goto 20

            r = sqrt(r2)
            r3 = r*r2
            r5 = r3*r2
            r3inv = inv4pi/r3
            do idim = 1,nd

               if (irotlet .eq. 1) then
                  rlet(1) = rotlet(idim,1,j)
                  rlet(2) = rotlet(idim,2,j)
                  rlet(3) = rotlet(idim,3,j)
                  rvec(1) = rotvec(idim,1,j)
                  rvec(2) = rotvec(idim,2,j)
                  rvec(3) = rotvec(idim,3,j)

                  rpl = zdiff(1)*rlet(1) + zdiff(2)*rlet(2) +
     1                  zdiff(3)*rlet(3)
                  rpv = zdiff(1)*rvec(1) + zdiff(2)*rvec(2) +
     1                  zdiff(3)*rvec(3)

                  pot(idim,1,i) = pot(idim,1,i) - rvec(1)*rpl*r3inv +
     1                            rlet(1)*rpv*r3inv
                  pot(idim,2,i) = pot(idim,2,i) - rvec(2)*rpl*r3inv +
     1                            rlet(2)*rpv*r3inv
                  pot(idim,3,i) = pot(idim,3,i) - rvec(3)*rpl*r3inv +
     1                            rlet(3)*rpv*r3inv

                  if (ifppreg .ge. 3) then
                  tempx1 = -(rlet(1) - 3.0d0*zdiff(1)*rpl/r2) * r3inv
                  tempx2 = -(rlet(2) - 3.0d0*zdiff(2)*rpl/r2) * r3inv
                  tempx3 = -(rlet(3) - 3.0d0*zdiff(3)*rpl/r2) * r3inv

                  grad(idim,1,1,i) = grad(idim,1,1,i) + rvec(1)*tempx1
                  grad(idim,2,1,i) = grad(idim,2,1,i) + rvec(1)*tempx2
                  grad(idim,3,1,i) = grad(idim,3,1,i) + rvec(1)*tempx3

                  grad(idim,1,2,i) = grad(idim,1,2,i) + rvec(2)*tempx1
                  grad(idim,2,2,i) = grad(idim,2,2,i) + rvec(2)*tempx2
                  grad(idim,3,2,i) = grad(idim,3,2,i) + rvec(2)*tempx3

                  grad(idim,1,3,i) = grad(idim,1,3,i) + rvec(3)*tempx1
                  grad(idim,2,3,i) = grad(idim,2,3,i) + rvec(3)*tempx2
                  grad(idim,3,3,i) = grad(idim,3,3,i) + rvec(3)*tempx3

                  tempx1 = (rvec(1) - 3.0d0*zdiff(1)*rpv/r2) * r3inv
                  tempx2 = (rvec(2) - 3.0d0*zdiff(2)*rpv/r2) * r3inv
                  tempx3 = (rvec(3) - 3.0d0*zdiff(3)*rpv/r2) * r3inv

                  grad(idim,1,1,i) = grad(idim,1,1,i) + rlet(1)*tempx1
                  grad(idim,2,1,i) = grad(idim,2,1,i) + rlet(1)*tempx2
                  grad(idim,3,1,i) = grad(idim,3,1,i) + rlet(1)*tempx3

                  grad(idim,1,2,i) = grad(idim,1,2,i) + rlet(2)*tempx1
                  grad(idim,2,2,i) = grad(idim,2,2,i) + rlet(2)*tempx2
                  grad(idim,3,2,i) = grad(idim,3,2,i) + rlet(2)*tempx3

                  grad(idim,1,3,i) = grad(idim,1,3,i) + rlet(3)*tempx1
                  grad(idim,2,3,i) = grad(idim,2,3,i) + rlet(3)*tempx2
                  grad(idim,3,3,i) = grad(idim,3,3,i) + rlet(3)*tempx3
                  endif
               endif

               if (idoublet .eq. 1) then
                  dlet(1) = doublet(idim,1,j)
                  dlet(2) = doublet(idim,2,j)
                  dlet(3) = doublet(idim,3,j)
                  dvec(1) = doubvec(idim,1,j)
                  dvec(2) = doubvec(idim,2,j)
                  dvec(3) = doubvec(idim,3,j)

                  dpl = zdiff(1)*dlet(1) + zdiff(2)*dlet(2) +
     1                  zdiff(3)*dlet(3)
                  dpv = zdiff(1)*dvec(1) + zdiff(2)*dvec(2) +
     1                  zdiff(3)*dvec(3)
                  dlv = dlet(1)*dvec(1) + dlet(2)*dvec(2) +
     1                  dlet(3)*dvec(3)

                  temp = 3.0d0*dpl*dpv*inv4pi/r5

                  pot(idim,1,i) = pot(idim,1,i) - zdiff(1)*dlv*r3inv -
     1                            dvec(1)*dpl*r3inv + dlet(1)*dpv*r3inv+
     2                            zdiff(1)*temp
                  pot(idim,2,i) = pot(idim,2,i) - zdiff(2)*dlv*r3inv -
     1                            dvec(2)*dpl*r3inv + dlet(2)*dpv*r3inv+
     2                            zdiff(2)*temp
                  pot(idim,3,i) = pot(idim,3,i) - zdiff(3)*dlv*r3inv -
     1                            dvec(3)*dpl*r3inv + dlet(3)*dpv*r3inv+
     2                            zdiff(3)*temp

                  if (ifppreg .ge. 2) then
                     pre(idim,i) = pre(idim,i) - 2.0d0*dlv*r3inv
     1                                         + 2.0d0*temp
                  endif

                  if (ifppreg .ge. 3) then
                  tempx1 = 3.0d0*(dlet(1)*dpv + dvec(1)*dpl -
     1                 5.0d0*zdiff(1)*dpl*dpv/r2)*inv4pi/r5
                  tempx2 = 3.0d0*(dlet(2)*dpv + dvec(2)*dpl -
     1                 5.0d0*zdiff(2)*dpl*dpv/r2)*inv4pi/r5
                  tempx3 = 3.0d0*(dlet(3)*dpv + dvec(3)*dpl -
     1                 5.0d0*zdiff(3)*dpl*dpv/r2)*inv4pi/r5

                  grad(idim,1,1,i) = grad(idim,1,1,i) + temp
                  grad(idim,1,1,i) = grad(idim,1,1,i) + zdiff(1)*tempx1
                  grad(idim,2,1,i) = grad(idim,2,1,i) + zdiff(1)*tempx2
                  grad(idim,3,1,i) = grad(idim,3,1,i) + zdiff(1)*tempx3

                  grad(idim,2,2,i) = grad(idim,2,2,i) + temp
                  grad(idim,1,2,i) = grad(idim,1,2,i) + zdiff(2)*tempx1
                  grad(idim,2,2,i) = grad(idim,2,2,i) + zdiff(2)*tempx2
                  grad(idim,3,2,i) = grad(idim,3,2,i) + zdiff(2)*tempx3

                  grad(idim,3,3,i) = grad(idim,3,3,i) + temp
                  grad(idim,1,3,i) = grad(idim,1,3,i) + zdiff(3)*tempx1
                  grad(idim,2,3,i) = grad(idim,2,3,i) + zdiff(3)*tempx2
                  grad(idim,3,3,i) = grad(idim,3,3,i) + zdiff(3)*tempx3

                  tempx1 = -(dlet(1) - 3.0d0*zdiff(1)*dpl/r2) * r3inv
                  tempx2 = -(dlet(2) - 3.0d0*zdiff(2)*dpl/r2) * r3inv
                  tempx3 = -(dlet(3) - 3.0d0*zdiff(3)*dpl/r2) * r3inv

                  grad(idim,1,1,i) = grad(idim,1,1,i) + dvec(1)*tempx1
                  grad(idim,2,1,i) = grad(idim,2,1,i) + dvec(1)*tempx2
                  grad(idim,3,1,i) = grad(idim,3,1,i) + dvec(1)*tempx3

                  grad(idim,1,2,i) = grad(idim,1,2,i) + dvec(2)*tempx1
                  grad(idim,2,2,i) = grad(idim,2,2,i) + dvec(2)*tempx2
                  grad(idim,3,2,i) = grad(idim,3,2,i) + dvec(2)*tempx3

                  grad(idim,1,3,i) = grad(idim,1,3,i) + dvec(3)*tempx1
                  grad(idim,2,3,i) = grad(idim,2,3,i) + dvec(3)*tempx2
                  grad(idim,3,3,i) = grad(idim,3,3,i) + dvec(3)*tempx3

                  tempx1 = (dvec(1) - 3.0d0*zdiff(1)*dpv/r2) * r3inv
                  tempx2 = (dvec(2) - 3.0d0*zdiff(2)*dpv/r2) * r3inv
                  tempx3 = (dvec(3) - 3.0d0*zdiff(3)*dpv/r2) * r3inv

                  grad(idim,1,1,i) = grad(idim,1,1,i) + dlet(1)*tempx1
                  grad(idim,2,1,i) = grad(idim,2,1,i) + dlet(1)*tempx2
                  grad(idim,3,1,i) = grad(idim,3,1,i) + dlet(1)*tempx3

                  grad(idim,1,2,i) = grad(idim,1,2,i) + dlet(2)*tempx1
                  grad(idim,2,2,i) = grad(idim,2,2,i) + dlet(2)*tempx2
                  grad(idim,3,2,i) = grad(idim,3,2,i) + dlet(2)*tempx3

                  grad(idim,1,3,i) = grad(idim,1,3,i) + dlet(3)*tempx1
                  grad(idim,2,3,i) = grad(idim,2,3,i) + dlet(3)*tempx2
                  grad(idim,3,3,i) = grad(idim,3,3,i) + dlet(3)*tempx3

                  tempx1 = 3.0d0*zdiff(1)*dlv*inv4pi/r5
                  tempx2 = 3.0d0*zdiff(2)*dlv*inv4pi/r5
                  tempx3 = 3.0d0*zdiff(3)*dlv*inv4pi/r5

                  grad(idim,1,1,i) = grad(idim,1,1,i) - dlv*r3inv
                  grad(idim,1,1,i) = grad(idim,1,1,i) + zdiff(1)*tempx1
                  grad(idim,2,1,i) = grad(idim,2,1,i) + zdiff(1)*tempx2
                  grad(idim,3,1,i) = grad(idim,3,1,i) + zdiff(1)*tempx3

                  grad(idim,2,2,i) = grad(idim,2,2,i) - dlv*r3inv
                  grad(idim,1,2,i) = grad(idim,1,2,i) + zdiff(2)*tempx1
                  grad(idim,2,2,i) = grad(idim,2,2,i) + zdiff(2)*tempx2
                  grad(idim,3,2,i) = grad(idim,3,2,i) + zdiff(2)*tempx3

                  grad(idim,3,3,i) = grad(idim,3,3,i) - dlv*r3inv
                  grad(idim,1,3,i) = grad(idim,1,3,i) + zdiff(3)*tempx1
                  grad(idim,2,3,i) = grad(idim,2,3,i) + zdiff(3)*tempx2
                  grad(idim,3,3,i) = grad(idim,3,3,i) + zdiff(3)*tempx3
                  endif
               endif

            enddo

 20         continue
         enddo
      enddo
c$omp end parallel do

      return

 100  continue

      return
      end

c**********************************************************************
      subroutine st3ddirectstokg_opt(nd,sources,stoklet,ns,targ,nt,
     1     ifppreg,pot,pre,grad,thresh)
cf2py  intent(in) nd,sources,stoklet,ns,targ,nt,ifppreg,thresh
cf2py  intent(out) pot,pre,grad
c
c     Optimized Stokeslet velocity with the st3ddirectstokg interface.
c     All outputs retain additive semantics. Only ifppreg=1 is optimized;
c     other output flags delegate to the original routine unchanged.
c     Source components are packed once per call. A shared reciprocal
c     distance and three scalar/vector accumulators replace repeated
c     divisions and output updates. Normalization is once per target.
c     nd=1 uses a SIMD source reduction (summation may be reordered).
c     nd>1 shares geometry and vectorizes contiguous densities, keeping
c     source summation order within each density. No fastmath is required.
c     OpenMP targets own their outputs; nested worker teams are avoided.
c
c$    use omp_lib, only: omp_in_parallel
      implicit none
      integer*8 nd,ns,nt,ifppreg,i,j,k
      real*8 sources(3,ns),stoklet(nd,3,ns),targ(3,nt)
      real*8 pot(nd,3,nt),pre(nd,nt),grad(nd,3,3,nt),thresh
      real*8, allocatable :: sx(:),sy(:),sz(:)
      real*8, allocatable :: fx(:,:),fy(:,:),fz(:,:)
      real*8, allocatable :: fx1(:),fy1(:),fz1(:)
      real*8, allocatable :: ax(:),ay(:),az(:)
      real*8 x,y,z,dx,dy,dz,r2,ri,ri3,dot,ux,uy,uz,threshsq
      real*8, parameter :: inv8pi =
     1    3.978873577297383394222094084313d-2

      if (ifppreg .ne. 1) then
         call st3ddirectstokg(nd,sources,stoklet,ns,targ,nt,
     1        ifppreg,pot,pre,grad,thresh)
         return
      endif
      if (ns .le. 0 .or. nt .le. 0 .or. nd .le. 0) return
      threshsq=thresh**2
      allocate(sx(ns),sy(ns),sz(ns))
      sx=sources(1,:)
      sy=sources(2,:)
      sz=sources(3,:)

      if (nd .eq. 1) then
         allocate(fx1(ns),fy1(ns),fz1(ns))
         fx1=stoklet(1,1,:)
         fy1=stoklet(1,2,:)
         fz1=stoklet(1,3,:)
c$omp parallel do default(shared) schedule(static)
c$omp& if (.not.omp_in_parallel() .and. nt*ns .ge. 4096)
c$omp& private(i,j,x,y,z,dx,dy,dz,r2,ri,ri3,dot,ux,uy,uz)
         do i=1,nt
            x=targ(1,i)
            y=targ(2,i)
            z=targ(3,i)
            ux=0.0d0
            uy=0.0d0
            uz=0.0d0
c$omp simd reduction(+:ux,uy,uz)
c$omp& private(dx,dy,dz,r2,ri,ri3,dot)
            do j=1,ns
               dx=x-sx(j)
               dy=y-sy(j)
               dz=z-sz(j)
               r2=dx*dx+dy*dy+dz*dz
               if (r2 .lt. threshsq) cycle
               ri=1.0d0/sqrt(r2)
               ri3=ri*ri*ri
               dot=(dx*fx1(j)+dy*fy1(j)+dz*fz1(j))*ri3
               ux=ux+fx1(j)*ri+dx*dot
               uy=uy+fy1(j)*ri+dy*dot
               uz=uz+fz1(j)*ri+dz*dot
            enddo
            pot(1,1,i)=pot(1,1,i)+inv8pi*ux
            pot(1,2,i)=pot(1,2,i)+inv8pi*uy
            pot(1,3,i)=pot(1,3,i)+inv8pi*uz
         enddo
c$omp end parallel do
         deallocate(fx1,fy1,fz1)
      else
         allocate(fx(nd,ns),fy(nd,ns),fz(nd,ns))
         fx=stoklet(:,1,:)
         fy=stoklet(:,2,:)
         fz=stoklet(:,3,:)
c     Allocatable accumulators are private and allocated once per thread,
c     not per target or pair. Densities remain contiguous in packed data.
c$omp parallel default(shared)
c$omp& if (.not.omp_in_parallel() .and. nt*ns .ge. 4096)
c$omp& private(ax,ay,az,i,j,k,x,y,z,dx,dy,dz,r2,ri,ri3,dot)
         allocate(ax(nd),ay(nd),az(nd))
c$omp do schedule(static)
         do i=1,nt
            x=targ(1,i)
            y=targ(2,i)
            z=targ(3,i)
            ax=0.0d0
            ay=0.0d0
            az=0.0d0
            do j=1,ns
               dx=x-sx(j)
               dy=y-sy(j)
               dz=z-sz(j)
               r2=dx*dx+dy*dy+dz*dz
               if (r2 .lt. threshsq) cycle
               ri=1.0d0/sqrt(r2)
               ri3=ri*ri*ri
c$omp simd private(dot)
               do k=1,nd
                  dot=(dx*fx(k,j)+dy*fy(k,j)+dz*fz(k,j))*ri3
                  ax(k)=ax(k)+fx(k,j)*ri+dx*dot
                  ay(k)=ay(k)+fy(k,j)*ri+dy*dot
                  az(k)=az(k)+fz(k,j)*ri+dz*dot
               enddo
            enddo
            do k=1,nd
               pot(k,1,i)=pot(k,1,i)+inv8pi*ax(k)
               pot(k,2,i)=pot(k,2,i)+inv8pi*ay(k)
               pot(k,3,i)=pot(k,3,i)+inv8pi*az(k)
            enddo
         enddo
c$omp end do
         deallocate(ax,ay,az)
c$omp end parallel
         deallocate(fx,fy,fz)
      endif
      deallocate(sx,sy,sz)
      return
      end
c
c
c
c
      subroutine st3ddirectstrsg_opt(nd,sources,strslet,strsvec,ns,
     1     targ,nt,ifppreg,pot,pre,grad,thresh)
cf2py  intent(in) nd,sources,strslet,strsvec,ns,targ,nt,ifppreg,thresh
cf2py  intent(out) pot,pre,grad
c
c     Optimized type I stresslet velocity, without a Stokeslet term:
c
c       pot(x) = pot(x) + sum_m T_{ijk}(x,y^{(m)}) mu^{(m)}_j nu^{(m)}_k
c
c     with the kernel, sign and scaling of st3ddirectstokstrsg (see top
c     of file); strslet holds mu and strsvec holds nu, both (nd,3,ns).
c     All outputs retain additive semantics. Only ifppreg=1 is
c     optimized; other output flags delegate to st3ddirectstokstrsg
c     with zero Stokeslet strengths, so their results are unchanged.
c     Source components are packed once per call. One reciprocal
c     distance per pair gives 1/r^5; three scalar/vector accumulators
c     replace per-pair output updates; 3/(4 pi) is applied once per
c     target. nd=1 uses a SIMD source reduction (summation may be
c     reordered). nd>1 shares geometry and vectorizes contiguous
c     densities, keeping source summation order within each density.
c     OpenMP targets own their outputs; nested worker teams are avoided.
c
c$    use omp_lib, only: omp_in_parallel
      implicit none
      integer*8 nd,ns,nt,ifppreg,i,j,k,istress
      real*8 sources(3,ns),strslet(nd,3,ns),strsvec(nd,3,ns)
      real*8 targ(3,nt)
      real*8 pot(nd,3,nt),pre(nd,nt),grad(nd,3,3,nt),thresh
      real*8, allocatable :: sx(:),sy(:),sz(:)
      real*8, allocatable :: mx(:,:),my(:,:),mz(:,:)
      real*8, allocatable :: vx(:,:),vy(:,:),vz(:,:)
      real*8, allocatable :: mx1(:),my1(:),mz1(:)
      real*8, allocatable :: vx1(:),vy1(:),vz1(:)
      real*8, allocatable :: ax(:),ay(:),az(:),zstok(:,:,:)
      real*8 x,y,z,dx,dy,dz,r2,ri,ri5,s,ux,uy,uz,threshsq
      real*8, parameter :: c3inv4pi =
     1    2.387324146378430036533256450587713d-1

      if (ifppreg .ne. 1) then
         allocate(zstok(nd,3,ns))
         zstok=0.0d0
         istress=1
         call st3ddirectstokstrsg(nd,sources,zstok,istress,
     1        strslet,strsvec,ns,targ,nt,ifppreg,pot,pre,grad,thresh)
         deallocate(zstok)
         return
      endif
      if (ns .le. 0 .or. nt .le. 0 .or. nd .le. 0) return
      threshsq=thresh**2
      allocate(sx(ns),sy(ns),sz(ns))
      sx=sources(1,:)
      sy=sources(2,:)
      sz=sources(3,:)

      if (nd .eq. 1) then
         allocate(mx1(ns),my1(ns),mz1(ns),vx1(ns),vy1(ns),vz1(ns))
         mx1=strslet(1,1,:)
         my1=strslet(1,2,:)
         mz1=strslet(1,3,:)
         vx1=strsvec(1,1,:)
         vy1=strsvec(1,2,:)
         vz1=strsvec(1,3,:)
c$omp parallel do default(shared) schedule(static)
c$omp& if (.not.omp_in_parallel() .and. nt*ns .ge. 4096)
c$omp& private(i,j,x,y,z,dx,dy,dz,r2,ri,ri5,s,ux,uy,uz)
         do i=1,nt
            x=targ(1,i)
            y=targ(2,i)
            z=targ(3,i)
            ux=0.0d0
            uy=0.0d0
            uz=0.0d0
c$omp simd reduction(+:ux,uy,uz)
c$omp& private(dx,dy,dz,r2,ri,ri5,s)
            do j=1,ns
               dx=x-sx(j)
               dy=y-sy(j)
               dz=z-sz(j)
               r2=dx*dx+dy*dy+dz*dz
               if (r2 .lt. threshsq) cycle
               ri=1.0d0/sqrt(r2)
               ri5=ri*ri
               ri5=ri5*ri5*ri
               s=(dx*mx1(j)+dy*my1(j)+dz*mz1(j))*
     1           (dx*vx1(j)+dy*vy1(j)+dz*vz1(j))*ri5
               ux=ux+dx*s
               uy=uy+dy*s
               uz=uz+dz*s
            enddo
            pot(1,1,i)=pot(1,1,i)+c3inv4pi*ux
            pot(1,2,i)=pot(1,2,i)+c3inv4pi*uy
            pot(1,3,i)=pot(1,3,i)+c3inv4pi*uz
         enddo
c$omp end parallel do
         deallocate(mx1,my1,mz1,vx1,vy1,vz1)
      else
         allocate(mx(nd,ns),my(nd,ns),mz(nd,ns))
         allocate(vx(nd,ns),vy(nd,ns),vz(nd,ns))
         mx=strslet(:,1,:)
         my=strslet(:,2,:)
         mz=strslet(:,3,:)
         vx=strsvec(:,1,:)
         vy=strsvec(:,2,:)
         vz=strsvec(:,3,:)
c     Allocatable accumulators are private and allocated once per thread,
c     not per target or pair. Densities remain contiguous in packed data.
c$omp parallel default(shared)
c$omp& if (.not.omp_in_parallel() .and. nt*ns .ge. 4096)
c$omp& private(ax,ay,az,i,j,k,x,y,z,dx,dy,dz,r2,ri,ri5,s)
         allocate(ax(nd),ay(nd),az(nd))
c$omp do schedule(static)
         do i=1,nt
            x=targ(1,i)
            y=targ(2,i)
            z=targ(3,i)
            ax=0.0d0
            ay=0.0d0
            az=0.0d0
            do j=1,ns
               dx=x-sx(j)
               dy=y-sy(j)
               dz=z-sz(j)
               r2=dx*dx+dy*dy+dz*dz
               if (r2 .lt. threshsq) cycle
               ri=1.0d0/sqrt(r2)
               ri5=ri*ri
               ri5=ri5*ri5*ri
c$omp simd private(s)
               do k=1,nd
                  s=(dx*mx(k,j)+dy*my(k,j)+dz*mz(k,j))*
     1              (dx*vx(k,j)+dy*vy(k,j)+dz*vz(k,j))*ri5
                  ax(k)=ax(k)+dx*s
                  ay(k)=ay(k)+dy*s
                  az(k)=az(k)+dz*s
               enddo
            enddo
            do k=1,nd
               pot(k,1,i)=pot(k,1,i)+c3inv4pi*ax(k)
               pot(k,2,i)=pot(k,2,i)+c3inv4pi*ay(k)
               pot(k,3,i)=pot(k,3,i)+c3inv4pi*az(k)
            enddo
         enddo
c$omp end do
         deallocate(ax,ay,az)
c$omp end parallel
         deallocate(mx,my,mz,vx,vy,vz)
      endif
      deallocate(sx,sy,sz)
      return
      end
