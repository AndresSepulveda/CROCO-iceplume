! $Id: iceplume.F $
!
!======================================================================
!                  ICEPLUME for CROCO
!======================================================================
! Buoyant plume theory (BPT) model of subglacial discharge plumes
! coupled to CROCO through point sources.
!
! Original authorship:
!   - IcePlume for MITgcm: Tom Cowton (Cowton et al., 2015, JGR).
!   - ROMS-ICEPLUME: Chuning Wang (github.com/ChuningWang/roms-iceplume)
!     from which the plume physics below (iceplume_calc ... linint)
!     is ported with minimal changes (PLUME(ng)%x -> pl_x, no ng).
!   - ODEPACK (iceplume_opkd.F90): A.C. Hindmarsh, public domain.
! CROCO adaptation (coupling, MPI reductions, CROCO conventions):
!   see README_ICEPLUME.md. Differences with ROMS-ICEPLUME:
!   * ambient density recomputed with the plume EOS (in-situ UNESCO)
!     instead of the model density (CROCO rho is relative to rho0);
!   * plume temperature converted to potential temperature before
!     injection;
!   * along-front velocity for background melt uses v(i,j),v(i,j+1).
!
! Entry point: iceplume_main, called by the master thread once per
! baroclinic step from prestep3D_thread (step.F).
!======================================================================
!
#include "cppdefs.h"
#ifdef ICEPLUME

      subroutine iceplume_main
      implicit none
# include "param.h"
# include "scalars.h"
# include "sources.h"
# include "iceplume.h"
      integer is, k
      real ramp, fIni
!
! One-time initialisation (geometry of each source)
!
      if (iic.eq.ntstart) call iceplume_init
!
! Ambient profiles at the plume cell (owner node only)
!
      call iceplume_ambient
# ifdef ICEPLUME_DET_AVERAGE
      call iceplume_average
# endif
!
! Plume model (owner node only)
!
# ifdef ICEPLUME_RAMP
      ramp=tanh(dt*sec2day*float(iic-ntstart))
# else
      ramp=1.
# endif
      do is=1,Nsrc
        call iceplume_zero(is)
        if (pl_Isrc(is).gt.0 .and. SG_typ(is).ne.0) then
          call iceplume_owner_test(is, k)
          if (k.eq.1) then
            fIni=abs(SG_bar(is))*ramp
            call iceplume_calc(is, fIni, SG_trc(is,itemp),
     &                         SG_trc(is,isalt), SG_typ(is),
     &                         SG_dep(is), SG_len(is))
            call iceplume_diag_calc(is)
# ifdef ICEPLUME_DEBUG
            call iceplume_debug(is, fIni)
# endif
          endif
        endif
      enddo
!
! Share results with all MPI nodes
!
      call iceplume_reduce
      call iceplume_diag
      return
      end
!
!======================================================================
! Owner test: own=1 if the plume cell of source is is in the interior
! of the local (MPI) subdomain; il,jl are local indices.
!======================================================================
!
      subroutine iceplume_owner_test(is, own)
      implicit none
# include "param.h"
# include "iceplume.h"
      integer is, own, il, jl
      call iceplume_local_ij(pl_Isrc(is), pl_Jsrc(is), il, jl)
      own=0
# ifdef MPI
      if (il.ge.1 .and. il.le.Lmmpi .and.
     &    jl.ge.1 .and. jl.le.Mmmpi) own=1
# else
      if (il.ge.1 .and. il.le.Lm .and.
     &    jl.ge.1 .and. jl.le.Mm) own=1
# endif
      return
      end

      subroutine iceplume_local_ij(ig, jg, il, jl)
      implicit none
# include "param.h"
      integer ig, jg, il, jl
# ifdef MPI
      il=ig-iminmpi+1
      jl=jg-jminmpi+1
# else
      il=ig
      jl=jg
# endif
      return
      end
!
!======================================================================
! Initialisation: locate the wet plume cell and the plume direction.
!   PSOURCE (HMass): the source is a u/v face between land and water;
!     the plume cell is the wet cell, dir=+1 if water is on the
!     positive side of the face, -1 otherwise (as in ROMS).
!   PSOURCE_MASS (VMass): the plume cell is (Isrc,Jsrc), dir=0.
!======================================================================
!
      subroutine iceplume_init
      implicit none
# include "param.h"
# include "scalars.h"
# include "grid.h"
# include "sources.h"
# include "iceplume.h"
# ifdef MPI
      include 'mpif.h'
      real buf(5*Msrc)
      integer ierr
# endif
      integer is, ig, jg, il, jl, own, k
      real rw1, rw0
      do is=1,Msrc
        pl_Isrc(is)=0
        pl_Jsrc(is)=0
        pl_dir(is)=0.
        pl_dx(is)=0.
        pl_dy(is)=0.
# ifdef ICEPLUME_DET_AVERAGE
        do k=1,N
          pl_rhoAm0(is,k)=0.
        enddo
# endif
      enddo
      do is=1,Nsrc
        ig=Isrc(is)
        jg=Jsrc(is)
        call iceplume_local_ij(ig, jg, il, jl)
        pl_dir(is)=0.
# ifdef PSOURCE
        if (il.ge.1 .and. jl.ge.1 .and.
#  ifdef MPI
     &      il.le.Lmmpi+1 .and. jl.le.Mmmpi+1) then
#  else
     &      il.le.Lm+1 .and. jl.le.Mm+1) then
#  endif
          if (Dsrc(is).eq.0) then
            rw1=rmask(il,jl)
            rw0=rmask(il-1,jl)
          else
            rw1=rmask(il,jl)
            rw0=rmask(il,jl-1)
          endif
          if (rw1.lt.0.5 .and. rw0.gt.0.5) then
            if (Dsrc(is).eq.0) then
              ig=ig-1
            else
              jg=jg-1
            endif
            pl_dir(is)=-1.
          elseif (rw1.gt.0.5 .and. rw0.lt.0.5) then
            pl_dir(is)=1.
          endif
        endif
# endif
        pl_Isrc(is)=ig
        pl_Jsrc(is)=jg
        call iceplume_owner_test(is, own)
        if (own.eq.1) then
          call iceplume_local_ij(ig, jg, il, jl)
          if (Dsrc(is).eq.0) then
            pl_dx(is)=om_r(il,jl)
            pl_dy(is)=on_r(il,jl)
          else
            pl_dx(is)=on_r(il,jl)
            pl_dy(is)=om_r(il,jl)
          endif
        else
          pl_dir(is)=0.
        endif
      enddo
# ifdef MPI
!
! Only the owner knows dir, dx, dy: global sum.
!
      if (NNODES.gt.1) then
        do is=1,Msrc
          buf(is)=pl_dir(is)
          buf(Msrc+is)=pl_dx(is)
          buf(2*Msrc+is)=pl_dy(is)
        enddo
        call MPI_Allreduce(MPI_IN_PLACE, buf, 3*Msrc,
     &       MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
        do is=1,Msrc
          pl_dir(is)=buf(is)
          pl_dx(is)=buf(Msrc+is)
          pl_dy(is)=buf(2*Msrc+is)
        enddo
      endif
# endif
      MPI_master_only write(stdout,'(/1x,A)')
     &  'ICEPLUME: source  Iplume Jplume  dir   dx[m]   dy[m] typ'//
     &  '  Q[m3/s]  len[m]  dep[m]'
      do is=1,Nsrc
        MPI_master_only write(stdout,'(10x,3I7,F5.1,2F8.1,I4,3F8.1)')
     &    is, pl_Isrc(is), pl_Jsrc(is), pl_dir(is), pl_dx(is),
     &    pl_dy(is), SG_typ(is), SG_bar(is), SG_len(is), SG_dep(is)
      enddo
      return
      end
!
!======================================================================
! Ambient profiles at the plume cell, time level nstp (owner only).
!======================================================================
!
      subroutine iceplume_ambient
      implicit none
# include "param.h"
# include "scalars.h"
# include "grid.h"
# include "ocean3d.h"
# include "sources.h"
# include "iceplume.h"
      integer is, k, itrc, il, jl, own
      real pr, prRef, rhoRefA
      parameter (prRef=101.d3*1.d-4, rhoRefA=1020.)
      double precision ipl_rho
      external ipl_rho
      do is=1,Nsrc
        do k=0,N
          pl_zW(is,k)=0.
        enddo
        do k=1,N
          pl_zR(is,k)=0.
          do itrc=1,NT
            pl_trcAm(is,k,itrc)=0.
          enddo
        enddo
        call iceplume_owner_test(is, own)
        if (own.eq.1) then
          call iceplume_local_ij(pl_Isrc(is), pl_Jsrc(is), il, jl)
          do k=0,N
            pl_zW(is,k)=z_w(il,jl,k)
          enddo
          do k=1,N
            pl_zR(is,k)=z_r(il,jl,k)
            pl_dz(is,k)=pl_zW(is,k)-pl_zW(is,k-1)
            pl_sAm(is,k)=t(il,jl,k,nstp,isalt)
            pl_tpAm(is,k)=t(il,jl,k,nstp,itemp)
            pr=prRef+(abs(pl_zR(is,k))*rhoRefA*g)*1.d-4
            call sw_temp(pl_sAm(is,k), pl_tpAm(is,k), pr, prRef,
     &                   pl_tAm(is,k))
            do itrc=1,NT
              pl_trcAm(is,k,itrc)=t(il,jl,k,nstp,itrc)
            enddo
            if (Dsrc(is).eq.0) then
              pl_vAm(is,k)=0.5*(v(il,jl,k,nstp)+v(il,jl+1,k,nstp))
            else
              pl_vAm(is,k)=0.5*(u(il,jl,k,nstp)+u(il+1,jl,k,nstp))
            endif
            pl_wAm(is,k)=0.5*(We(il,jl,k-1)+We(il,jl,k))
     &                      *pm(il,jl)*pn(il,jl)
# ifndef ICEPLUME_DET_AVERAGE
            pl_rhoAm(is,k)=ipl_rho(pl_tAm(is,k), pl_sAm(is,k),
     &                             pl_zR(is,k))
# endif
          enddo
          pl_rhoAm(is,N+1)=1.225
          do itrc=1,NT
# ifdef ICEPLUME_TRACER
            pl_trcIni(is,itrc)=SG_trc(is,itrc)
# else
            pl_trcIni(is,itrc)=0.
# endif
          enddo
        endif
      enddo
      return
      end

# ifdef ICEPLUME_DET_AVERAGE
!
!======================================================================
! Ambient density averaged over the window SG_Irange x SG_Jrange
! (global indices, wet points only, area weighted), interpolated to
! the depths of the plume cell. Equivalent of the mp_assemble based
! averaging of ROMS-ICEPLUME v1.1.1, using MPI_Allreduce.
!======================================================================
!
      subroutine iceplume_average
      implicit none
#  include "param.h"
#  include "scalars.h"
#  include "grid.h"
#  include "ocean3d.h"
#  include "sources.h"
#  include "iceplume.h"
#  ifdef MPI
      include 'mpif.h'
      real buf(2*Msrc*N)
      integer ierr, m
#  endif
      integer is, k, i, j, il, jl, own, ilo, ihi, jlo, jhi
      real area, tp, sa, ti, pr, prRef, rhoRefA, avgAlphaA
      parameter (prRef=101.d3*1.d-4, rhoRefA=1020., avgAlphaA=1.0)
      real zc(N), tc(N), sc(N)
      double precision ipl_rho
      external ipl_rho
#  ifdef MPI
!
! Everybody needs the depths of the plume cell.
!
      if (NNODES.gt.1) then
        call MPI_Allreduce(MPI_IN_PLACE, pl_zR, Msrc*N,
     &       MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
      endif
#  endif
      do is=1,Nsrc
        do k=1,N
          pl_rhoSum(is,k)=0.
          pl_avgArea(is,k)=0.
        enddo
        if (SG_typ(is).gt.1) then
          do j=SG_Jrange(is,1),SG_Jrange(is,2)
            do i=SG_Irange(is,1),SG_Irange(is,2)
              call iceplume_local_ij(i, j, il, jl)
#  ifdef MPI
              if (il.ge.1 .and. il.le.Lmmpi .and.
     &            jl.ge.1 .and. jl.le.Mmmpi) then
#  else
              if (il.ge.1 .and. il.le.Lm .and.
     &            jl.ge.1 .and. jl.le.Mm) then
#  endif
               if (rmask(il,jl).gt.0.5) then
                area=om_r(il,jl)*on_r(il,jl)
                do k=1,N
                  zc(k)=z_r(il,jl,k)
                  tc(k)=t(il,jl,k,nstp,itemp)
                  sc(k)=t(il,jl,k,nstp,isalt)
                enddo
                do k=1,N
                  if (pl_zR(is,k).gt.zc(1)) then
                    call linint(N, zc, tc, pl_zR(is,k), tp)
                    call linint(N, zc, sc, pl_zR(is,k), sa)
                    pr=prRef+(abs(pl_zR(is,k))*rhoRefA*g)*1.d-4
                    call sw_temp(sa, tp, pr, prRef, ti)
                    pl_rhoSum(is,k)=pl_rhoSum(is,k)
     &                   +area*ipl_rho(ti, sa, pl_zR(is,k))
                    pl_avgArea(is,k)=pl_avgArea(is,k)+area
                  endif
                enddo
               endif
              endif
            enddo
          enddo
        endif
      enddo
#  ifdef MPI
      if (NNODES.gt.1) then
        m=Msrc*N
        do k=1,N
          do is=1,Msrc
            buf(is+(k-1)*Msrc)=pl_rhoSum(is,k)
            buf(m+is+(k-1)*Msrc)=pl_avgArea(is,k)
          enddo
        enddo
        call MPI_Allreduce(MPI_IN_PLACE, buf, 2*m,
     &       MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
        do k=1,N
          do is=1,Msrc
            pl_rhoSum(is,k)=buf(is+(k-1)*Msrc)
            pl_avgArea(is,k)=buf(m+is+(k-1)*Msrc)
          enddo
        enddo
      endif
#  endif
      do is=1,Nsrc
        call iceplume_owner_test(is, own)
        if (own.eq.1 .and. SG_typ(is).gt.1) then
          do k=1,N
            pl_rhoAm0(is,k)=pl_rhoAm(is,k)
            if (pl_avgArea(is,k).gt.0.) then
              pl_rhoAm(is,k)=pl_rhoSum(is,k)/pl_avgArea(is,k)
            else
              pl_rhoAm(is,k)=ipl_rho(pl_tAm(is,k), pl_sAm(is,k),
     &                               pl_zR(is,k))
            endif
            if (avgAlphaA.gt.0. .and. avgAlphaA.lt.1. .and.
     &          iic.gt.ntstart) then
              pl_rhoAm(is,k)=pl_rhoAm0(is,k)*(1.-avgAlphaA)
     &                      +pl_rhoAm(is,k)*avgAlphaA
            endif
          enddo
        elseif (own.eq.1) then
          do k=1,N
            pl_rhoAm(is,k)=ipl_rho(pl_tAm(is,k), pl_sAm(is,k),
     &                             pl_zR(is,k))
          enddo
        endif
      enddo
      return
      end
# endif /* ICEPLUME_DET_AVERAGE */
!
!======================================================================
! Reset the fluxes of source is (all nodes).
!======================================================================
!
      subroutine iceplume_zero(is)
      implicit none
# include "param.h"
# include "iceplume.h"
      integer is, k, itrc
      do k=1,N
        pl_ent(is,k)=0.
        pl_det(is,k)=0.
        pl_mB(is,k)=0.
        do itrc=1,NT
          pl_trcB(is,k,itrc)=0.
        enddo
      enddo
      do itrc=1,NT
        pl_trc(is,itrc)=0.
      enddo
      pl_trs(is)=0.
      pl_zOS(is)=0.
      pl_zNB(is)=0.
      pl_meltP(is)=0.
      pl_meltB(is)=0.
      return
      end
!
!======================================================================
! Diagnostics of source is (owner): overshoot depth (top of the
! rising plume), neutral buoyancy depth (level of maximum
! detrainment), plume melt and background melt [m3/s].
!======================================================================
!
      subroutine iceplume_diag_calc(is)
      implicit none
# include "param.h"
# include "iceplume.h"
      integer is, k, kmax
      real dmax
      kmax=0
      do k=1,N
        if (pl_f(is,k).gt.0.) pl_zOS(is)=pl_zW(is,k)
      enddo
      dmax=0.
      do k=1,N
        if (pl_det(is,k).gt.dmax) then
          dmax=pl_det(is,k)
          kmax=k
        endif
      enddo
      if (kmax.gt.0) pl_zNB(is)=pl_zR(is,kmax)
      pl_meltP(is)=0.
      do k=0,N
        pl_meltP(is)=max(pl_meltP(is), pl_mInt(is,k))
      enddo
      pl_meltB(is)=0.
      do k=1,N
        pl_meltB(is)=pl_meltB(is)+pl_mB(is,k)
      enddo
      return
      end
!
!======================================================================
! Global sum of the plume fluxes: only the owner has non-zero values,
! so after the sum every node (and every tile) sees the same fluxes.
! Needed because a source face may lie on a subdomain boundary.
!======================================================================
!
      subroutine iceplume_reduce
      implicit none
# include "param.h"
# include "scalars.h"
# include "iceplume.h"
# ifdef MPI
      include 'mpif.h'
      integer nbuf
      parameter (nbuf=Msrc*(3*N+NT+2*N*NT+5))
      real buf(nbuf)
      integer is, k, itrc, m, ierr
      if (NNODES.le.1) return
      m=0
      do is=1,Msrc
        do k=1,N
          buf(m+1)=pl_ent(is,k)
          buf(m+2)=pl_det(is,k)
          buf(m+3)=pl_mB(is,k)
          m=m+3
          do itrc=1,NT
            buf(m+1)=pl_trcAm(is,k,itrc)
            buf(m+2)=pl_trcB(is,k,itrc)
            m=m+2
          enddo
        enddo
        do itrc=1,NT
          m=m+1
          buf(m)=pl_trc(is,itrc)
        enddo
        buf(m+1)=pl_trs(is)
        buf(m+2)=pl_zOS(is)
        buf(m+3)=pl_zNB(is)
        buf(m+4)=pl_meltP(is)
        buf(m+5)=pl_meltB(is)
        m=m+5
      enddo
      call MPI_Allreduce(MPI_IN_PLACE, buf, m,
     &     MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
      m=0
      do is=1,Msrc
        do k=1,N
          pl_ent(is,k)=buf(m+1)
          pl_det(is,k)=buf(m+2)
          pl_mB(is,k)=buf(m+3)
          m=m+3
          do itrc=1,NT
            pl_trcAm(is,k,itrc)=buf(m+1)
            pl_trcB(is,k,itrc)=buf(m+2)
            m=m+2
          enddo
        enddo
        do itrc=1,NT
          m=m+1
          pl_trc(is,itrc)=buf(m)
        enddo
        pl_trs(is)=buf(m+1)
        pl_zOS(is)=buf(m+2)
        pl_zNB(is)=buf(m+3)
        pl_meltP(is)=buf(m+4)
        pl_meltB(is)=buf(m+5)
        m=m+5
      enddo
# endif
      return
      end
!
!======================================================================
! ASCII diagnostics (node 0): file iceplume_diag.txt, every NINFO
! steps: time [days], source, Qsg, net transport, overshoot and
! neutral buoyancy depths, plume melt, background melt [m3/s].
!======================================================================
!
      subroutine iceplume_diag
      implicit none
# include "param.h"
# include "scalars.h"
# include "sources.h"
# include "iceplume.h"
      integer is, iu
      parameter (iu=91)
      if (mynode.ne.0) return
      if (iic.eq.ntstart) then
        open(iu, file='iceplume_diag.txt', status='replace')
        write(iu,'(A)') '# time[days] src Qsg[m3/s] trs[m3/s]'//
     &   ' zOS[m] zNB[m] meltPlume[m3/s] meltBkg[m3/s]'
        close(iu)
      endif
      if (mod(iic-ntstart,ninfo).ne.0) return
      open(iu, file='iceplume_diag.txt', position='append')
      do is=1,Nsrc
        write(iu,'(F12.5,I5,2F12.4,2F9.2,2ES14.5)')
     &    time*sec2day, is, SG_bar(is), pl_trs(is), pl_zOS(is),
     &    pl_zNB(is), pl_meltP(is), pl_meltB(is)
      enddo
      close(iu)
      return
      end

# ifdef ICEPLUME_DEBUG
!
!======================================================================
! Debug: budget of the fluxes delivered by the plume (owner node).
! Net tracer flux F = sum_k (det*Cp + ent*Ca + mB*CB) [tracer m3/s].
!======================================================================
!
      subroutine iceplume_debug(is, fIni)
      implicit none
#  include "param.h"
#  include "scalars.h"
#  include "iceplume.h"
      integer is, k, itrc
      real fIni, ipl_sdet, ipl_sent, ipl_smb, flx(NT)
      ipl_sdet=0.
      ipl_sent=0.
      ipl_smb=0.
      do itrc=1,NT
        flx(itrc)=0.
      enddo
      do k=1,N
        ipl_sdet=ipl_sdet+pl_det(is,k)
        ipl_sent=ipl_sent+pl_ent(is,k)
        ipl_smb=ipl_smb+pl_mB(is,k)
        do itrc=1,NT
          flx(itrc)=flx(itrc)+pl_det(is,k)*pl_trc(is,itrc)
     &           +pl_ent(is,k)*pl_trcAm(is,k,itrc)
     &           +pl_mB(is,k)*pl_trcB(is,k,itrc)
        enddo
      enddo
      write(stdout,'(A,I6,5ES12.4,4ES12.4)') ' IPLDBG ', iic, fIni,
     &   ipl_sdet, ipl_sent, ipl_smb, pl_trs(is),
     &   (flx(itrc),itrc=1,min(NT,4))
      return
      end
# endif
!
!======================================================================
!  Plume physics, ported from ROMS-ICEPLUME (C. Wang)
!======================================================================
!
!
!=======================================================================
!                                                                      !
! This subroutine call the iceplume core functions to calcualte plume  !
!   status/meltrates/etc.                                              ! 
!                                                                      !
!=======================================================================
!
      SUBROUTINE iceplume_calc(I,
     &                         fIni, tIni, sIni0,
     &                         sgTyp, sgDep, sgLen)
!
      implicit none
#include "param.h"
#include "iceplume.h"
#include "iceplume_const.h"
      double precision ipl_rho
      external ipl_rho
!
!=======================================================================
!                                                                      !
! Input variables:                                                     !
!                                                                      !
!=======================================================================
!                                                                      !
! ng      - grid identifier                                            !
! I       - point source index identifier                              !
! fIni    - subglacial runoff volume flux [m^3 s^-1]                   !
! tIni    - subglacial runoff temperature [degC]                       !
! sIni    - subglacial runoff salinity [PSU]                           !
! sgTyp   - subglacial runoff type identifier                          !
! sgDep   - subglacial runoff inject depth [m]                         !
! sgLen   - subglacial runoff discharge length [m]                     !
!                                                                      !
!=======================================================================
!                                                                      !
! Local variables:                                                     !
!                                                                      !
!=======================================================================
!                                                                      !
! iceDepthK   - vertical grid index of grounding line depth            !
! plumeDepthK - vertical grid index of plume depth                     !
! osDepthK    - vertical grid index of overshoot depth                 !
! areaP       - surface area of plume in contact with ice in           !
!               that cell [m^2]                                        !
! areaB       - surface area of plume in contact with ice out          !
!               that cell [m^2]                                        !
! meanVel     - ice tangental velocity [m s^-1]                        !
! det         - total detrainment volume [m^3 s^-1]                    !
! mB          - background meltrate [m s^-1]                           !
!                                                                      !
!=======================================================================
!
! In/out variables
!
      integer, intent(in)  :: I
      real, intent(in) :: fIni, tIni, sIni0
      integer, intent(in)  :: sgTyp
      real, intent(in) :: sgDep, sgLen
!
! Local variables declaration
!
      integer :: iceDepthK, plumeDepthK, osDepthK
      real :: sIni
      real :: areaP, areaB, meanVel, det, mB, sB
!
      integer :: K, itrc
      real :: cff, cff1, cff2
!
! If discharge input salinity is less or equal than zero, set it to
! a small number
!
      IF (sIni0 .LE. 0.0) THEN
        sIni = 0.001
      ELSE
        sIni = sIni0
      ENDIF
!
!=======================================================================
!                                                                      !
! Find iceDepthK                                                       !
!                                                                      !
!=======================================================================
!
      iceDepthK = -1
      IF (sgDep .GE. 0.0) THEN
        iceDepthK = 0
      ELSE
        DO K = 0, N
          IF (pl_zW(I, K) .GE. sgDep) THEN
            iceDepthK = K
            EXIT
          ENDIF
        ENDDO
      ENDIF
!
!=======================================================================
!                                                                      !
! Use entrainment plume model to calculate T, S, & flux in plume       !
!                                                                      !
!=======================================================================
!
! Initialize and prepare profiles
!
      plumeDepthK = iceDepthK+1
      det = 0.0
!
      DO k = 0, N
        pl_f(I, K)    = 0.0
        pl_w(I, K)    = 0.0
        pl_t(I, K)    = 0.0
        pl_s(I, K)    = 0.0
        pl_a(I, K)    = 0.0
        pl_mInt(I, K) = 0.0
        pl_rho(I, K)  = 0.0
        pl_lm(I, K)   = 0.0
        pl_lc(I, K)   = 0.0
      ENDDO
!
      DO k = 1, N
        pl_ent(I, K) = 0.0
        pl_det(I, K) = 0.0
!
        pl_detI(I, K)    = 0
        pl_detFrac(I, K) = 0.0
      ENDDO
!
      IF ( (fIni .GT. 0.0) .AND. (sgTyp .NE. 1) ) THEN
!
! Run the plume model
!
        CALL iceplume_entrain(I, iceDepthK,
     &                        fIni, tIni, sIni,
     &                        sgTyp, sgLen)
!
! Calculate volume flux differential to give entrainment/extrainment
! First clear ent/det
!
        DO K = iceDepthK+1, N
          pl_ent(I, K) =
     &        -(pl_f(I, K) - pl_f(I, K-1))
          pl_ent(I, K) = pl_ent(I, K) +
     &        MAX(pl_mInt(I, K) -
     &            pl_mInt(I, K-1), 0.0)
          pl_ent(I, K) = MIN(pl_ent(I, K), 0.0)
        ENDDO
!
! Separate entrainment/detrainment
!
        DO K = iceDepthK+1, N
          IF (pl_f(I, K) .GT. 0.0) THEN
            osDepthK = K
            det = pl_f(I, K)
          ENDIF
        ENDDO
!
! Replace the last value of lc and lm
!
        pl_lm(I, osDepthK) = pl_lm(I, osDepthK-1)
        pl_lc(I, osDepthK) = pl_lc(I, osDepthK-1)
!
! Find detrainment depth index, Calculate total detrainment volume
! and thickness
!
# if defined ICEPLUME_SHOOT_NO && defined ICEPLUME_DET_NO
        plumeDepthK = osDepthK
# else
        DO K = osDepthK, iceDepthK+1, -1
          cff = ipl_rho(pl_t(I, osDepthK),
     &              pl_s(I, osDepthK),
     &              pl_zR(I, K))
          IF (K .EQ. 1) THEN
            cff1 = pl_rhoAm(I, 1)
            cff2 = 0.5*(pl_rhoAm(I, 1)+
     &                  pl_rhoAm(I, 2))
          ELSEIF (K .EQ. N) THEN
            cff1 = 0.5*(pl_rhoAm(I, N-1)+
     &                  pl_rhoAm(I, N))
            cff2 = pl_rhoAm(I, N)
          ELSE
            cff1 = 0.5*(pl_rhoAm(I, K-1)+pl_rhoAm(I, K))
            cff2 = 0.5*(pl_rhoAm(I, K)+pl_rhoAm(I, K+1))
          ENDIF
          IF ((K .EQ. N) .AND. (cff .LT. cff2)) THEN
            plumeDepthK = N
            EXIT
          ENDIF
          IF ((cff .LT. cff1) .AND. (cff .GT. cff2)) THEN
            plumeDepthK = K
            EXIT
          ENDIF
        ENDDO
# endif
!
!=======================================================================
!                                                                      !
! Use a detrainment model to distribute the detrainment plume in       !
! several layers. Details in iceplume_detrain.F.                       !
!                                                                      !
! *** This part is completely rewritten since Feb 20, 2019.            !
!                                                                      !
!=======================================================================
!
# ifdef ICEPLUME_DET_NO
        pl_det(I, plumeDepthK) = det
        pl_detI(I, plumeDepthK) = 1
        pl_detFrac(I, plumeDepthK) = 1.0
# else
        CALL iceplume_detrain(I,
     &                        iceDepthK, plumeDepthK, osDepthK,
     &                        pl_lc(I, osDepthK),
     &                        det)
# endif
      ENDIF
# ifdef ICEPLUME_MELT
!
!=======================================================================
!                                                                      !
! Calculate melt rates.                                                !
!                                                                      !
!=======================================================================
!
! Initiate
!
      DO K = 1, N
        pl_m(I, K)  = 0.0
        pl_mB(I, K) = 0.0
      ENDDO
!
      DO K = iceDepthK+1, N
!
!=======================================================================
!                                                                      !
! If there is a plume in that cell, then need to calculate plume       !
! meltrate distinct to background melt rate. Plume melt rate is        !
! already encorporated in the plume model, and taken into account      !
! in the temperature and salinity of the plume outflow. It is  useful  !
! though to have it available as a diagnostic.                         !
!                                                                      !
!=======================================================================
!
        areaP = pl_a(I, K) - pl_a(I, K-1)
        areaP = MAX(areaP, 0.0)
        areaB = pl_dz(I, K)*pl_dy(I) - areaP
        areaB = MAX(areaB, 0.0)
!
        IF ( areaP .GT. 0.0 ) THEN
          pl_m(I, K) = pl_mInt(I, K) -
     &                          pl_mInt(I, K-1)
        ENDIF
!
!=======================================================================
!                                                                      !
! Calculate the background melt rate (not generated by plumes).        !
! This will then be used to update the temperature and salinity in     !
! the adjacent cells. Velocities are calculated at cell faces -        !
! find averages for cell centres. Does not include velocity            !
! perpendicular to ice-this differs depending on orientation of        !
! ice front.                                                           !
!                                                                      !
!=======================================================================
!
        IF ( areaB .GT. 0.0 )THEN
!
          meanVel = SQRT((pl_vAm(I, K))**2 +
     &                   (pl_wAm(I, K))**2)
          CALL iceplume_meltrate(pl_tAm(I, K),
     &                           pl_sAm(I, K),
     &                           meanVel,
     &                           pl_zR(I, K),
     &                           mB, pl_tB(I, K), sB)
!
! Calcate integrated background melt rate [m3 s-1]
!
          pl_mB(I, K) = mB*areaB
        ENDIF
      ENDDO
# endif
!
!=======================================================================
!                                                                      !
! Calculate active and passive tracers                                 !
!                                                                      !
! Update log                                                           !
! Add active tracers, T and S into this scheme. This is necessary      !
! when distributing detrainment in several layers, since writting      !
! T & S in PLUME directly is in general not a good idea.               !
! 2018/4/20 Chuning Wang                                               !
!                                                                      !
!=======================================================================
!
! Clear local plume tracer variables
!
      DO itrc = 1, NT
        pl_trc(I, itrc)    = 0.0
        pl_trcCum(I, itrc) = 0.0
!
        DO K = iceDepthK+1, N
          pl_trcB(I, K, itrc) = 0.0
        ENDDO
      ENDDO
!
      IF ( (fIni .GT. 0.0) .AND. (sgTyp .NE. 1) ) THEN
!
! Write the active tracer (T & S) concentration to PLUME (ng) % trc
!
        pl_trc(I, itemp) = pl_t(I, osDepthK)
!  CROCO: the model temperature is potential temperature; convert
!  the plume in-situ temperature at the overshoot level (reference
!  pressure = surface). ROMS-ICEPLUME injects in-situ temperature.
        cff1 = 101.d3*1.d-4 + ABS(pl_zW(I, osDepthK))*rhoRef*g*1.d-4
        CALL sw_ptmp(pl_s(I, osDepthK), pl_t(I, osDepthK), cff1,
     &               101.d3*1.d-4, cff2)
        pl_trc(I, itemp) = cff2
        pl_trc(I, isalt) = pl_s(I, osDepthK)
!
# ifdef ICEPLUME_TRACER
!
! If use passive tracers, calculate passive tracer concentrations
! in the runoff
!
        DO itrc = 3, NT
          pl_trcCum(I, itrc) = pl_trcIni(I, itrc)*fIni
        ENDDO
!
! Add up total sum of each tracer in plume
!
        DO K = iceDepthK+1, osDepthK
          IF (pl_ent(I, K) .LT. 0.0) THEN
            DO itrc = 3, NT
              pl_trcCum(I, itrc) =
     &            pl_trcCum(I, itrc) +
     &            (-pl_ent(I, K) *
     &              pl_trcAm(I, K, itrc))
            ENDDO
          ENDIF
        ENDDO
!
! Calculate concentration of tracer in outflow 
!
        DO itrc = 3, NT
          pl_trc(I, itrc) = pl_trcCum(I, itrc)/det
        ENDDO
#  ifdef ICEPLUME_MELT_TRACER
!
! If use melt water tracer, rewrite tracer concentration
! Calculate accumulated trcer amount from base to detrain depth
!
        pl_trcCum(I, NT-1) = 0.0
        pl_trcCum(I, NT  ) = 0.0
        DO K = iceDepthK+1, osDepthK
          pl_trcCum(I, NT-1) =
     & pl_trcCum(I, NT-1) +
     & (-pl_ent(I, K)*pl_trcAm(I, K, NT-1)) +
     & pl_m(I, K)*pl_trcIni(I, NT-1)
           pl_trcCum(I, NT  ) =
     & pl_trcCum(I, NT  ) +
     & (-pl_ent(I, K)*pl_trcAm(I, K, NT  ))
        ENDDO
!
! Calculate tracer concentration
!
        pl_trc(I,NT-1) = pl_trcCum(I,NT-1)/det
        pl_trc(I,NT  ) = pl_trcCum(I,NT  )/det
#  endif
# endif
      ENDIF  ! (fIni .GT. 0.0) .AND. (sgTyp .NE. 1)
# ifdef ICEPLUME_MELT
!
! Rewrite tracer concentration in background melt water 
!
      DO K = iceDepthK+1, N
        pl_trcB(I, K, itemp) =
     &     (cI*(tIce-pl_tB(I, K)) - L)/cW
        pl_trcB(I, K, isalt) = sIce
#  ifdef ICEPLUME_MELT_TRACER
!
! If use melt water tracer, rewrite the last tracer type.
!
        pl_trcB(I, K, NT) = pl_trcIni(I, NT)
#  endif
      ENDDO
# endif
!
! Calculate depth integrated volume transport [m3 s-1]
!
      pl_trs(I) = 0.0
      DO K = 1, N
        pl_trs(I) = pl_trs(I) +
     &     pl_det(I, K) + pl_ent(I, K) +
     &     pl_mB(I, K)
      ENDDO
      END SUBROUTINE iceplume_calc
!
!=======================================================================
!                                                                      !
! This function calculates background melt rate.                       !
!                                                                      !
!=======================================================================
!
      SUBROUTINE iceplume_meltrate(temp, salt, vel, depth, mdot, sB, tB)
!
      implicit none
#include "param.h"
#include "iceplume.h"
#include "iceplume_const.h"
      double precision ipl_rho
      external ipl_rho
!
      real :: temp, salt, vel
      real :: depth, absVelocity
      real :: a, b, c
      real, intent(inout) :: mdot, sB, tB
!
! Routine cannot cope with zero velocity.
! Unlikely to occur anyway with currents, waves, convection etc.
! This is not very physical, but will do for now.
!
      IF ( ABS(vel) .LT. velBkg )  vel = velBkg
!
! Calculate melt rate from 3 equation formualtion (as for plume
! models)
! Equations for Sb, Tb and mdot
!
      a = lambda1*(GamT*cW-GamS*cI)
      b = GamS*cI*(lambda1*salt-lambda2-lambda3*depth +
     &             tIce - (L/cI)) -
     &    GamT*cW*(temp-lambda2-lambda3*depth+lambda1*sIce)
      c = GamS*salt*(cI*(lambda2+lambda3*depth-tIce) + L) +
     &     GamT*sIce*cW*(temp-lambda2-lambda3*depth)
!
      sB   = (1.0/(2.0*a))*(-b-SQRT(b**2 - 4.0*a*c))
      tB   = lambda1*sB+lambda2+lambda3*depth
      mdot = GamS*SQRT(CdBkg)*ABS(vel)*(salt-sB)/sB
!
      END SUBROUTINE iceplume_meltrate
!
!=======================================================================
!                                                                      !
! These are the core functions of iceplume entrainment model.          !
!                                                                      !
!=======================================================================
!
      SUBROUTINE iceplume_entrain(I, iceDepthK,
     &                            fIni, tIni, sIni,
     &                            sgTyp, sgLen)
!
      implicit none
#include "param.h"
#include "iceplume.h"
#include "iceplume_const.h"
      double precision ipl_rho
      external ipl_rho
      external general_entrain_model, jex
!
!=======================================================================
!                                                                      !
! Input variables:                                                     !
!                                                                      !
!=======================================================================
!                                                                      !
! ng          - grid identifier                                        !
! I           - point source index identifier                          !
! iceDepthK   - vertical grid index of grounding line depth            !
! fIni        - subglacial runoff volume flux [m^3 s^-1]               !
! tIni        - subglacial runoff temperature [degC]                   !
! sIni        - subglacial runoff salinity [PSU]                       !
! sgTyp       - subglacial runoff type identifier                      !
! sgLen       - subglacial runoff discharge length [m]                 !
!                                                                      !
!=======================================================================
!                                                                      !
! Local variables:                                                     !
!                                                                      !
!=======================================================================
!                                                                      !
! rhoA/tA/sA  - ambient density/temperature/salinity                   !
! zIn/zOut    - in/out depth of the model steps                        !
!                                                                      !
!=======================================================================
!
! In/out variables
!
      integer, intent(in)  :: I, iceDepthK
      real, intent(in) :: fIni, tIni, sIni
      integer, intent(in)  :: sgTyp
      real, intent(in) :: sgLen
!
! Local variables declaration
!
      real :: rhoA, tA, sA
      real :: zIn, zOut
      integer  :: K
!
! Local variables for ODEPACK
!
      integer :: IOPT, ISTATE, ITASK, ITOL, IWORK(20), LIW, LRW, MF, NEQ
      real, parameter :: RTOL = 1.0E-5
      real, parameter :: ATOL = 1.0E-5
      real :: RWORK(116), Y(10)
!
!=======================================================================
!                                                                      !
! Y is input/output vector for DLSODE                                  !
!   Y(1)  - plume volume flux                                          !
!   Y(2)  - plume velocity                                             !
!   Y(3)  - plume temperature                                          !
!   Y(4)  - plume salinity                                             !
!   Y(5)  - integrated plume/glacier contact area                      !
!   Y(6)  - area integrated melt                                       !
!   Y(7)  - fake variable to pass in ng                                !
!   Y(8)  - fake variable to pass in I                                 !
!   Y(9)  - fake variable to pass in sgTyp                             !
!   Y(10) - fake variable to pass in sgLen                             !
!                                                                      !
!=======================================================================
!                                                                      !
! For ODEPACK solver. See ODEPACK documentation and source code in     !
! Cowton et al. 2015.                                                  !
!                                                                      !
!=======================================================================
!
      NEQ      = 6
      LRW      = 116
      LIW      = 116
!
      ITOL     = 1
      ITASK    = 1
      ISTATE   = 1
      IOPT     = 0
      MF       = 10
      IWORK(7) = 2  ! To limit # of repeat error messages are printed
!
! Initial conditions
!
      Y(1)  = fIni  ! initial plume volume flux
      Y(2)  = wIni  ! initial vertical velocity
      Y(3)  = tIni  ! initial temperature
      Y(4)  = sIni  ! initial salinity
      Y(5)  = 0.0   ! integrated contact area
      Y(6)  = 0.0   ! integrated melt rate
!
      Y(7)  = 0.0          ! fake variable to pass in ng
      Y(8)  = REAL(I)      ! fake variable to pass in I
      Y(9)  = REAL(sgTyp)  ! fake variable to pass in sgTyp
      Y(10) = sgLen        ! fake variable to pass in sgLen
!
! Start at bottom of ice face
!
      zIn = pl_zW(I, iceDepthK)
      zOut = pl_zW(I, iceDepthK+1)
!
! Set initial conditions
!
      pl_f(I, iceDepthK)    = Y(1)
      pl_w(I, iceDepthK)    = Y(2)
      pl_t(I, iceDepthK)    = Y(3)
      pl_s(I, iceDepthK)    = Y(4)
      pl_a(I, iceDepthK)    = Y(5)
      pl_mInt(I, iceDepthK) = Y(6)
      pl_rho(I, iceDepthK)  = ipl_rho(Y(3), Y(4), zIn)
!
      CALL plume_metrics (sgTyp, Y(1)/Y(2), sgLen,
     &                    pl_lm(I, iceDepthK),
     &                    pl_lc(I, iceDepthK))
!
! Move up through water column from lowest layer
!
      DO K = iceDepthK+1, N
!
!=======================================================================
!                                                                      !
! Use DLSODE to solve plume properties.                                !
!                                                                      !
!=======================================================================
!
! Check to make sure plume has not reached neutral buoyancy in a lower
! layer
!
        IF (ISTATE .GT. -1) THEN
          CALL DLSODE (GENERAL_ENTRAIN_MODEL,
     &                 NEQ, Y, zIn, zOut,
     &                 ITOL, RTOL, ATOL, ITASK,
     &                 ISTATE, IOPT, RWORK, LRW, IWORK,
     &                 LIW, JEX, MF)
!
          CALL plume_metrics (sgTyp, Y(1)/Y(2), sgLen,
     &                        pl_lm(I, K),
     &                        pl_lc(I, K))
!
!=======================================================================
!                                                                      !
! Test to see if neutral buoyancy has now been reached. If solver      !
! returns ISTATE = -1, then it has been unable to meet required        !
! tolerances at this level. This generally occurs because plume        !
! has reached neutral buoyancy and run out of momentum, and so is      !
! no longer rising. At this point, we therefore end the call to the    !
! plume model. Our aim is to catch the plume at the point of           !
! neutral buoyancy. We therefore perform a manual comparrison of       !
! ambient and plume density. If plume density >= ambient density we    !
! assign ISTATE = -1, again ending the call to the plume model.        !
!                                                                      !
!=======================================================================
!
! Calculate plume density
!
          pl_rho(I, K) = ipl_rho(Y(3), Y(4), zIn)
# ifdef ICEPLUME_SHOOT_NO
!
! Calculate ambient density
!
          IF (K .EQ. 0) THEN
            tA = pl_tAm(I, 1)
            sA = pl_sAm(I, 1)
          ELSEIF (K .EQ. N) THEN
            tA = pl_tAm(I, N)
            sA = pl_sAm(I, N)
          ELSE
            tA = 0.5*(pl_tAm(I, K) + pl_tAm(I, K+1))
            sA = 0.5*(pl_sAm(I, K) + pl_sAm(I, K+1))
          ENDIF
          rhoA = ipl_rho(tA, sA, zIn)
          IF ((pl_rho(I, K) .GT. rhoA) .OR. (K .EQ. N)) THEN
# else
!
          IF ((Y(2) .LE. 0.0) .OR. (K .EQ. N)) THEN
# endif
            ISTATE = -1
          ENDIF
!
! If the plume has not reached neutral buoyancy, then we assign a depth
! at which to calculate the next value and loop round to call the plume
! model again. Make sure we are not at the surface
!
          IF (ISTATE .GT. -1) THEN
            zIn = zOut
            zOut = pl_zW(I, K+1)
          ENDIF
!
        ELSE  ! (ISTATE .LE. -1)
!
! This section is entered once the plume has reached neutral buoyancy
! once plume has reached neutral buoyancy, no plume values
!
          Y(1) = 0.0
          Y(2) = 0.0
          Y(3) = 0.0
          Y(4) = 0.0
          Y(5) = 0.0
          Y(6) = 0.0
        ENDIF
!
! Save results.
!
        pl_f(I, K) = Y(1)
        pl_w(I, K) = Y(2)
        pl_t(I, K) = Y(3)
        pl_s(I, K) = Y(4)
        pl_a(I, K) = Y(5)
        pl_mInt(I, K) = Y(6)
!
      ENDDO
      RETURN
      END SUBROUTINE iceplume_entrain
!
!=======================================================================
!                                                                      !
! Below is the generalized plume model.                                !
!                                                                      !
!=======================================================================
!
      SUBROUTINE general_entrain_model (NEQ, T, Y, YDOT)
!
      implicit none
#include "param.h"
#include "iceplume.h"
#include "iceplume_const.h"
      double precision ipl_rho
      external ipl_rho
!
      integer, intent(in) :: NEQ
      integer :: I, sgTyp
      real :: sgLen
      real :: T, Y(10), YDOT(6)
!
      real :: tA, sA, tB, sB, rhoP, rhoA, gRed
      real :: mdot, a, b, c
      real :: Lmc, Lcc
!
! Check if velcoity is positive
!
      IF (Y(2) .LE. 0.0) THEN
        YDOT(1) = 0.0
        YDOT(2) = 0.0
        YDOT(3) = 0.0
        YDOT(4) = 0.0
        YDOT(5) = 0.0
        YDOT(6) = 0.0
      ELSE
!
! Interpolate from imposed ambient profiles
!
        I     = NINT(Y(8))
        sgTyp = NINT(Y(9))
        sgLen = Y(10)
!
        CALL linint(N, pl_zR(I, :),
     &              pl_tAm(I, :), T, tA)
        CALL linint(N, pl_zR(I, :),
     &              pl_sAm(I, :), T, sA)
!
! Calculate reduced gravity
!
        rhoP = ipl_rho(Y(3), Y(4), T)
        rhoA = ipl_rho(tA, sA, T)
        gRed  = g*(rhoA-rhoP)/rhoRef
!
! Calculate plume metrics
!
        CALL plume_metrics (sgTyp, Y(1)/Y(2), sgLen, Lmc, Lcc)
!
! Equations for sB, tB and mdot
!
        a = lambda1*(GamT*cW - GamS*cI)
        b = GamS*cI*(lambda1*Y(4) - lambda2 - lambda3*T +
     &               tIce - (L/cI)) -
     &      GamT*cW*(Y(3) - lambda2 - lambda3*T + lambda1*sIce)
        c = GamS*Y(4)*(cI*(lambda2 + lambda3*T - tIce) + L) +
     &      GamT*sIce*cW*(Y(3) - lambda2 - lambda3*T)
!
        sB   = (1.0/(2.0*a))*(-b - SQRT(b**2 - 4.0*a*c))
        tB   = lambda1*sB + lambda2 + lambda3*T
        mdot = GamS*SQRT(Cd)*Y(2)*(Y(4) - sB)/sB
!
! Plume volume flux
!
        YDOT(1) = alpha*Lcc*Y(2) + Lmc*mdot
!
! Plume vertical velocity
!
        YDOT(2) = (1./Y(1))*(-Y(2)*YDOT(1) + gRed*Y(1)/Y(2) -
     &                       Cd*Lmc*Y(2)**2)
!
! Plume temperature
!
        YDOT(3) = (1./Y(1))*(-Y(3)*YDOT(1) + alpha*Lcc*Ta*Y(2) +
     &                       Lmc*mdot*tB -
     &                       SQRT(Cd)*GamT*Lmc*Y(2)*(Y(3)-tB))
!
! Plume salinity
!
        YDOT(4) = (1./Y(1))*(-Y(4)*YDOT(1) + alpha*Lcc*Sa*Y(2) +
     &                       Lmc*mdot*sB -
     &                       SQRT(Cd)*GamS*Lmc*Y(2)*(Y(4)-sB))
!
! Along-plume integrated contact area and melt rate
!
        YDOT(5) = Lmc
        YDOT(6) = Lmc*mdot
      ENDIF
!
      END SUBROUTINE general_entrain_model
!
!=======================================================================
!                                                                      !
! Use this function to calculate plume metrics.                        !
!                                                                      !
!=======================================================================
!
      SUBROUTINE plume_metrics (sgTyp, area, Ls, Lm, Lc)
!
! Calculate Lm (plume contact length with ice) and C (plume contact
! length with water), based on various of plume metrics.
!
! For each method, refer to the manual for details.
!
! If you wish to add new plume metrics, assign a new sgTyp number
! and code it to a new CASE.
!
      implicit none
#include "iceplume_const.h"
      double precision ipl_rho
      external ipl_rho
!
      integer, intent(in)   :: sgTyp
      real, intent(in)  :: area, Ls
      real, intent(out) :: Lm, Lc
!
      real :: gam, a, b, h
!
      SELECT CASE (sgTyp)
        CASE (2)
!
! Halfcone (Cowton)
!
          Lc = SQRT(2.*pi*area)
          Lm = 2./pi*Lc
!
        CASE (3)
!
! Finite line source
!
          Lc = SQRT(Ls**2 + 2.*pi*area)
          Lm = Ls + 2./pi*(Lc - Ls)
!
        CASE (4)
!
! Sheet (Jenkins)
!
          Lc = Ls
          Lm = Ls
!
        CASE (5)
!
! Fullcone
!
          Lc = 2.*SQRT(2.*pi*area)
          Lm = 0.
!
        CASE (6)
!
! Elipse (Aspect ratio 0.5)
!
          gam = 0.5
          b = SQRT(2.*area/(gam*pi))
          a = gam*b
          h = ((a-b)/(a+b))**2
          Lc = 0.5*pi*(a+b)*(1.+3.*h/(10.+SQRT(4.-3.*h)))
          Lm = 2.*b
!
        CASE DEFAULT
!
! Halfcone (Cowton)
!
          Lc = SQRT(2.*pi*area)
          Lm = 2./pi*Lc
!
      END SELECT
!
      END SUBROUTINE plume_metrics
!
!=======================================================================
!                                                                      !
! Dummy routine for ODEPACK. Necessary for Jacobian matrix if          !
! stiff ODEs.                                                          !
!                                                                      !
!=======================================================================
!
      SUBROUTINE jex()
        RETURN
      END SUBROUTINE jex
!
!=======================================================================
!                                                                      !
! Use this function to calculate the velocity of detrainment.          !
!                                                                      !
!=======================================================================
!
      SUBROUTINE iceplume_detrain(I,
     &                            iceDepthK, plumeDepthK, osDepthK,
     &                            lc, det)
!
      implicit none
#include "param.h"
#include "iceplume.h"
#include "iceplume_const.h"
      double precision ipl_rho
      external ipl_rho
!
!=======================================================================
!                                                                      !
! Input/Output variables:                                              !
!                                                                      !
!=======================================================================
!                                                                      !
! ng          - grid identifier                                        !
! I           - point source index identifier                          !
! iceDepthK   - glacier grounding line depth vertical grid index       !
! plumeDepthK - plume detrainment depth vertical grid index            !
! osDepthK    - plume overshoot depth vertical grid index              !
! lc          - plume/ocean contact length [m]                         !
! det         - detrainment volume flux [m^3 s^-1]                     !
!                                                                      !
!=======================================================================
!                                                                      !
! Local variables:                                                     !
!                                                                      !
!=======================================================================
!                                                                      !
! detDz       - detrainment layer total thickness [m]                  !
! minDz       - minimum detrainment layer thickness [m]                !
! rhoP        - plume density [kg m^-3]                                !
! potE        - potential energy of plume/ambient water                !
!                                                                      !
! detVel      - detrainment velocity [m s^-1]                          !
! detWeight   - detrainment weight function                            !
! detSum      - detrainment volume flux before normalize [m^3 s^-1]    !
!                                                                      !
! isSurface   - logical, if the plume discharge at surface             !
! isBottom    - logical, if the plume discharge at bottom              !
! udSwitch    - switch to determine search one layer up(1)/down(-1)    !
!                                                                      !
! maxVel      - maximum velocity based on Richardson number            !
!               criterion [m s^-1]                                     !
! Fr          - Froude number                                          !
! rho1        - upper layer density [kg m^-3]                          !
! rho2        - lower layer density [kg m^-3]                          !
! h1          - upper layer thickness [m]                              !
! h2          - lower layer thickness [m]                              !
! gRed        - reduced gravity [m s^-2]                               !
# ifdef ICEPLUME_DET_NEUTRAL
!                                                                      !
! volR        - detrainment volume of runoff [m^3]                     !
! volM        - detrainment volume of meltwater                        !
! volF        - detrainment volume of total freshwater                 !
! volE        - detrainment volume of entrainment                      !
! volP        - detrainment volume of plume                            !
!                                                                      !
! tF          - temperature of total freshwater [degC]                 !
! tE          - temperature of entrainment                             !
! tP          - temperature of plume                                   !
!                                                                      !
! sF          - salinity of total freshwater [PSU]                     !
! sE          - salinity of entrainment                                !
! sP          - salinity of plume                                      !
!                                                                      !
! rhoF        - density of total freshwater [kg m^-3]                  !
! rhoE        - density of entrainment                                 !
! rhoP        - density of plume                                       !
# endif
!                                                                      !
!=======================================================================
!
! In/out variables
!
      integer, intent(in)  :: I
      integer, intent(in)  :: iceDepthK, plumeDepthK, osDepthK
      real, intent(in) :: lc, det
!
! Local variables declaration
!
      real :: detDz, minDz, detVel, detWeight, potE
      real :: rhoP
      integer  :: udSwitch
      logical  :: isSurface = .FALSE.
      logical  :: isBottom = .FALSE.
!
! For detrainment weight function
!
      real :: detSum
!
! For calculation of internal wave speed
!
      real :: maxVel, Fr, Ri
      real :: rho1, rho2, h1, h2, gRed
      real :: rhoH, h01, h02
# ifdef ICEPLUME_DET_NEUTRAL
!
! For variable detrainment tracer concentration
!
      real :: volR, volM, volF, volE, volP
      real :: tF, tE, tP, sF, sE, sP, rhoF, rhoE
# endif
!
! Other local variables
!
      integer  :: KI, K, K1, K2, itrc, counter
      real :: cff, cff1, cff2, cff3, cff4
# ifdef ICEPLUME_DET_SMOOTH
      real :: potEi, potEp, potEn, potERatio
# endif
# ifdef ICEPLUME_DET_NEUTRAL
!
! Initialize.
!
      DO K = 1, N
        pl_detF(I, K) = 0.0
        pl_detE(I, K) = 0.0
      ENDDO
      DO itrc = 1, NT
        DO K = 1, N
          pl_detTrc(I, K, itrc) = 0.0
        ENDDO
      ENDDO
# endif
!
!=======================================================================
!                                                                      !
! PART I - Find minimum detrainment depth                              !
!                                                                      !
!=======================================================================
!
! Initiate logical variables
!
      IF (plumeDepthK .EQ. N) THEN
        isSurface = .TRUE.
      ELSE
        isSurface = .FALSE.
      ENDIF
      IF (plumeDepthK .EQ. 1) THEN
        isBottom = .TRUE.
      ELSE
        isBottom = .FALSE.
      ENDIF
!
! Compute minimum detrainment thickness.
!
      rhoP = ipl_rho(pl_t(I, osDepthK),
     &           pl_s(I, osDepthK),
     &           pl_zR(I, plumeDepthK))
!
! Use the Parameterization of Ching et al. (1993) to determine the bulk
! velocity and nt in several layers. First, calculate BV frequency.
!
      IF (isSurface) THEN
!
! Strong density jump. Always use momentum conservation.
!
        maxVel = 0.95*pl_w(I, N)
      ELSE
!
! Use piecewise function.
!
        KI = plumeDepthK
        rhoH = 0.0
        h1 = 0.0
        DO K = KI+1, N
          h1 = h1+pl_dz(I, K)
          rhoH = rhoH + pl_dz(I, K)*pl_rhoAm(I, K)
        ENDDO
        rho1 = rhoH/h1
!
        rhoH = 0.0
        h2 = 0.0
        DO K = 1, KI 
          h2 = h2+pl_dz(I, K)
          rhoH = rhoH + pl_dz(I, K)*pl_rhoAm(I, K)
        ENDDO
        rho2 = rhoH/h2
!
        gRed = MAX(-g*(rho1-rho2)/rhoRef, gRedBkg)
        KI = MAXLOC(pl_w(I, :), 1)
        cff1 = pl_f(I, KI)/pl_w(I, KI)/lc
        Ri = gRed*cff1/(pl_w(I, KI)**2)
# ifdef ICEPLUME_DET_SMOOTH
!
! Add a correction term to generate a smooth transition from subsurface
! plume to surface plume. The idea is to nudge the Richardson number
! towards 6.0 (thus the outflow behaves more like a surface plume) when
! the vertical momentum at surface is non-zero.
!
! This is an engineering fix in order to address Becca's concerns when
! running iceplume with a realistic topography.
!
! Chuning Wang 2021-05-17
!
        cff = 0.0
        DO K = 0,N
          IF (pl_w(I, K) .gt. cff) cff = pl_w(I, K)
        ENDDO
        cff1 = pl_w(I, N) / cff
        Ri = Ri + cff1*(6.0 - Ri)
# endif
!
        IF (Ri .LT. 6.0) THEN
          maxVel = (0.7*(Ri**0.17))*pl_w(I, KI)
        ELSE
          maxVel = 0.95*pl_w(I, KI)
        ENDIF
      ENDIF
      maxVel = maxVel*(SQRT(pi/2)*detSigma)
      minDz = det/lc/maxVel
!
! Determine the vertical extent of detrainment plume.
!
      KI = plumeDepthK
      K1 = KI
      K2 = KI
      detDz = 0.0
!
! Add layers around plumeDepthK until it reaches critical thickness.
!
      counter = 0
      pl_detI(I, KI) = 1
      pl_detFrac(I, KI) = 1.0
      detDz = pl_dz(I, KI)
      potE = g*pl_dz(I, KI)*(pl_rhoAm(I, KI)-rhoP)
      IF (potE .GT. 0.0) THEN
        udSwitch = 1
      ELSE
        udSwitch = -1
      ENDIF
      DO WHILE ( (detDz .LT. minDz) .AND. (counter .LT. 100) )
        counter = counter + 1
        IF (udSwitch .EQ. 1) THEN
!
! Search one layer up
!
          DO K = 2, N
            IF ((pl_detI(I, K  ) .EQ. 0) .AND.
     &          (pl_detI(I, K-1) .EQ. 1)) THEN
              KI = K
            ENDIF
          ENDDO
          K1 = KI
        ELSEIF (udSwitch .EQ. -1) THEN
!
! Search one layer down
!
          DO K = 1, N-1
            IF ((pl_detI(I, K+1) .EQ. 1) .AND.
     &          (pl_detI(I, K  ) .EQ. 0)) THEN
              KI = K
            ENDIF
          ENDDO
          K2 = KI
        ENDIF
        IF (KI .EQ. N) isSurface = .TRUE.
        IF (KI .EQ. 1) isBottom = .TRUE.
!
! Update detrainment layer flag and detrainment depth.
!
        pl_detI(I, KI) = 1
        detDz = detDz + pl_dz(I, KI)
        pl_detFrac(I, KI) = 1.0
!
! Determine next layer to minimize potential energy anomaly.
!
        potE = potE +
     &         g*pl_dz(I, KI)*(pl_rhoAm(I, KI)-rhoP)
        IF (potE .GT. 0.0) THEN
          udSwitch = 1
        ELSE
          udSwitch = -1
        ENDIF
!
! IF layer reached surface/bottom, force udSwitch to -1/1
!
        IF     ( (       isSurface ) .AND. ( .NOT. isBottom ) ) THEN
          udSwitch = -1
        ELSEIF ( ( .NOT. isSurface ) .AND. (       isBottom ) ) THEN
          udSwitch = 1
        ELSEIF ( (       isSurface ) .AND. (       isBottom ) ) THEN
          udSwitch = 0
        ENDIF
!
! Other exit loop criteria.
!
        IF ( (minDz .EQ. (minDz+1.0)) .OR. (udSwitch .EQ. 0) ) EXIT
      ENDDO
!
! Only take partial grid from top/bottom layer to satisfy potE=0
!
      IF ( (.NOT. isSurface) .AND. (counter .GT. 0) ) THEN
        IF ( potE .GT. 0.0 ) THEN
          KI = K2
        ELSE
          KI = K1
        ENDIF
        cff = 1 - potE / (g*(pl_rhoAm(I, KI)-rhoP)) /
     &            pl_dz(I, KI)
        IF ( (cff .GE. 0) .AND. (cff .LT. 1) ) THEN
          pl_detFrac(I, KI) = cff
          detDz = detDz - pl_dz(I, KI)*(1 - cff)
        ENDIF
      ENDIF
# ifdef ICEPLUME_DET_SMOOTH
!
! Make correction for plume shape towards surface.
!
! Chuning Wang 2021-05-17
!
      potEp = 0.0
      potEn = 0.0
      DO K = K2,K1
        potEi = g*pl_detFrac(I, KI)*pl_dz(I, K)*
     &          (pl_rhoAm(I, K)-rhoP)
        IF (potEi .gt. 0.0) THEN
          potEp = potEp + potEi
        ELSE
          potEn = potEn - potEi
        ENDIF
      ENDDO
      IF ((potEp .gt. 0.0) .and. (potEn .gt. 0.0)) THEN
        potERatio = potEp/potEn
      ELSE
        potERatio = 1.0
      ENDIF
# endif
!
!=======================================================================
!                                                                      !
! PART II - Distribute detrainmnet in several layers                   !
!                                                                      !
! Update log                                                           !
! Use a Gause function to smooth the distribution.                     !
! 2018/06/08 Chuning Wang                                              !
!                                                                      !
!=======================================================================
!
! Update detrainment volume flux
!
      detVel = det / (lc * detDz)
      DO K = 1,N
        IF (pl_detI(I, K) .EQ. 1) THEN
!
! First, calculate weight using Gaussian function.
!
          h1 = pl_zR(I, K1)-pl_zR(I, plumeDepthK)+0.0001
          h2 = pl_zR(I, K2)-pl_zR(I, plumeDepthK)+0.0001
# ifdef ICEPLUME_DET_SMOOTH
          h1 = h1*potERatio
# endif
          IF ( K .GT. plumeDepthK ) THEN
            cff = (pl_zR(I, K)-pl_zR(I, plumeDepthK))/h1
          ELSE
            cff = (pl_zR(I, K)-pl_zR(I, plumeDepthK))/h2
          ENDIF
          detWeight = EXP(-0.5 * (cff/detSigma)**2)
!
! Distribute the detrainment.
!
          pl_det(I, K) = detWeight * detVel * lc *
     &      pl_detFrac(I, K) * pl_dz(I, K)
        ENDIF
      ENDDO
!
! Normalize
!
      detSum = SUM(pl_det(I, :))
      DO K = 1,N
        pl_det(I, K) = pl_det(I, K)*det/detSum
      ENDDO
# ifdef ICEPLUME_DET_NEUTRAL
!
! Calculate proportion of each water mass.
! There are three main water masses -
! 1. Runoff (R),
! 2. Melt (M),
! 3. Entrainment (E).
! R and M combined forms a new water mass,
! 4. Freshwater (F).
! M and F combined froms the final water mass,
! 5. Plume (P).
!
      volP = pl_f(I, osDepthK)
      volR = pl_f(I, iceDepthK)
      volM = pl_mInt(I, osDepthK)
      volF = volR + volM 
      volE = volP - volF
!
! Temperature. Note that tGade is effective temperature of Gade line.
!
      tP = pl_trc(I, itemp)
      tF = (tGade*volM + pl_t(I, iceDepthK)*volR)/volF
      tE = (tP*volP-tF*volF)/volE
!
! Salinity.
!
      sP = pl_trc(I, isalt)
      sF = (sIce*volM  + pl_s(I, iceDepthK)*volR)/volF
      sE = (sP*volP-sF*volF)/volE
!
! Density.
!
      rhoF = ipl_rho(tF, sF, pl_zR(I, plumeDepthK))
      rhoE = (rhoP*volP-rhoF*volF)/volE
!
! Calculate proportion from mixing line.
!
      DO K = 1, N
        IF (pl_detI(I, K) .EQ. 1) THEN
          cff = (pl_zR(I, K) -
     &           pl_zR(I, plumeDepthK))*0.005
          cff1 = MIN(MAX((pl_rhoAm(I, K)-rhoE)/(rhoF-rhoE),
     &                   0.01), 0.99)
          cff1 = ((pl_rhoAm(I, K)+cff)-rhoE)/(rhoF-rhoE)
          cff1 = MIN(MAX(cff1, 0.01), 0.99)
          cff2 = 1-cff1
          pl_detF(I, K) = pl_det(I, K)*cff1
          pl_detE(I, K) = pl_det(I, K)*cff2
        ENDIF
      ENDDO
!
! Make corrections to make sure each water mass sums to its previous
! volume.
!
      cff1 = volF/SUM(pl_detF(I, :))
      cff2 = volE/SUM(pl_detE(I, :))
      DO K = 1, N
        IF (pl_detI(I, K) .EQ. 1) THEN
          pl_detF(I, K) = pl_detF(I, K)*cff1
          pl_detE(I, K) = pl_detE(I, K)*cff2
        ENDIF
      ENDDO
!
! Rewrite detrainment volume.
!
      DO K = 1, N
        IF (pl_detI(I, K) .EQ. 1) THEN
          pl_det(I, K) = pl_detF(I, K) +
     &                            pl_detE(I, K)
        ENDIF
      ENDDO
!
!=======================================================================
!                                                                      !
! Update tracer concentration in detrainment model.                    !
!                                                                      !
!=======================================================================
!
! Active tracers.
!
      DO K = 1, N
        IF (pl_detI(I, K) .EQ. 1) THEN
          pl_detTrc(I, K, itemp) =
     &        (tF*pl_detF(I, K)+tE*pl_detE(I, K)) /
     &        (pl_detF(I, K)+pl_detE(I, K))
          pl_detTrc(I, K, isalt) =
     &        (sF*pl_detF(I, K)+sE*pl_detE(I, K)) /
     &        (pl_detF(I, K)+pl_detE(I, K))
        ENDIF
      ENDDO
!
! Passive tracers.
!
      DO itrc = 3, NT
        cff1 = pl_trcIni(I, itrc)*volR/volF
#  ifdef ICEPLUME_MELT_TRACER
        IF ( itrc .EQ. NT-1 ) THEN
          cff1 = pl_trcIni(I, itrc)*volM/volF
        ELSEIF ( itrc .EQ. NT ) THEN
          cff1 = 0.0
        ENDIF
#  endif
        cff2 = (pl_trc(I, itrc)*volP - cff1*volF)/volE
        DO K = 1, N
          IF (pl_detI(I, K) .EQ. 1) THEN
            pl_detTrc(I, K, itrc) =
     &        (cff1*pl_detF(I, K)+cff2*pl_detE(I, K)) /
     &        (pl_detF(I, K)+pl_detE(I, K))
          ENDIF
        ENDDO
      ENDDO
# endif
      END SUBROUTINE iceplume_detrain
!
!=======================================================================
!                                                                      !
! Use this function to calculate seawater density.                     !
!                                                                      !
!=======================================================================
!
      double precision FUNCTION ipl_rho(T,S,z)
!
! Equation of state (UNESCO 1983)
!     T = temperature (deg C)
!     S = salinity (PSU)
!     z = depth (m)
!
      double precision :: T,S,z
      double precision :: P
      double precision :: kw, Aw, Bw, k0
      double precision :: bulk_modulus
      double precision :: A, B, C, rho_w,rho_zero
!
      double precision, parameter :: rho_0 = 1027.0
      double precision, parameter :: g = 9.81
!

      P = rho_0*g*abs(z)*1.0e-5
!
! calculate in-situ rho
      kw = 19652.21 + 148.4206*T - 2.327105*T**2 + 1.360477e-2*(T**3) -
     &     5.155288e-5*(T**4)
      Aw = 3.239908 + 1.43713e-3*T + 1.16092e-4*T**2 - 5.77905e-7*T**3
      Bw = 8.50935e-5 - 6.12293e-6*T + 5.2787e-8*(T**2)
      k0 = kw + (54.6746 - 0.603459*T + 1.09987e-2*(T**2) -
     &           6.1670e-5*(T**3))*S +
     &     (7.944e-2 + 1.6483e-2*T - 5.3009e-4*(T**2))*(S**1.5)
      A = Aw + (2.2838e-3 - 1.0981e-5*T - 1.6078e-6*(T**2))*S +
     &    1.91075e-4*(S**1.5)
      B = Bw + (-9.9348e-7 + 2.0816e-8*T + 9.1697e-10*T**2)*S
      bulk_modulus = k0 + A*P + B*P**2
!
      A = 8.24493e-1 - 4.0899e-3*T + 7.6438e-5*T**2 - 8.2467e-7*T**3 +
     &    5.3875e-9*T**4
      B = -5.72466e-3 + 1.0227e-4*T - 1.6546e-6*T**2
      C = 4.8314e-4
      rho_w = 999.842594 + 6.793952e-2*T - 9.095290e-3*T**2+
     &        1.001685e-4*T**3 - 1.120083e-6*T**4 + 6.536336e-9*T**5
      rho_zero = rho_w + A*S + B*(S**1.5) + C*(S**2)
!
      ipl_rho = rho_zero/(1.0 - (P/bulk_modulus))
!
      END FUNCTION ipl_rho
!
!=======================================================================
!                                                                      !
! These subroutines are modified from MITgcm.                          !
!                                                                      !
!=======================================================================
!
      SUBROUTINE sw_temp(S, T, P, PR, rv)
!
! *=============================================================*
! | S/R  SW_TEMP
! | o compute in-situ temperature from potential temperature
! *=============================================================*
!
! REFERENCES:
! Fofonoff, P. and Millard, R.C. Jr
! Unesco 1983. Algorithms for computation of fundamental properties of
! seawater, 1983. _Unesco Tech. Pap. in Mar. Sci._, No. 44, 53 pp.
! Eqn.(31) p.39
! 
! Bryden, H. 1973.
! New Polynomials for thermal expansion, adiabatic temperature
! gradient and potential temperature of sea water.
! DEEP-SEA RES., 1973, Vol20, 401-408.
!
!
!
! === Global variables ===
! 
! !INPUT/OUTPUT PARAMETERS:
! === Routine arguments ===
! S      :: salinity
! T      :: potential temperature
! P      :: pressure
! PR     :: reference pressure
! rv:: return value (in-situ temeparture in degree C)
!
      real :: S, T, P, PR
      real :: rv
!
! !LOCAL VARIABLES:
! === local variables ===
!
      CALL sw_ptmp  (S, T, PR, P, rv)
!
      END SUBROUTINE sw_temp
!
      SUBROUTINE sw_ptmp  (S, T, P, PR, rv)
!
! !DESCRIPTION: \bv
! *=============================================================*
! | S/R  SW_PTMP
! | o compute potential temperature as per UNESCO 1983 report.
! *=============================================================*
! \ev
! started:
!          Armin Koehl akoehl@ucsd.edu
!
! ==================================================================
! SUBROUTINE SW_PTMP
!
! ==================================================================
!
!
!
! === Global variables ===
!
! !INPUT/OUTPUT PARAMETERS:
! === Routine arguments ===
! S  :: salinity    [psu      (PSS-78) ]
! T  :: temperature [degree C (IPTS-68)]
! P  :: pressure    [db]
! PR:: Reference pressure  [db]
! rv:: return value (potential temeparture in degree C)
!
      real :: S, T, P, PR
      real :: rv
!
! !LOCAL VARIABLES
! === local variables ===
!
      real :: del_P, del_th, th, q
      real :: onehalf, two, three
      parameter ( onehalf = 0.5, two = 2.0, three = 3.0 )
      real :: adtg_val
!
! theta1
!
      del_P   = PR-P
      CALL sw_adtg(S, T, P, adtg_val)
      del_th  = del_P*adtg_val
      th      = T+onehalf*del_th
      q       = del_th
!
! theta2
!
      CALL sw_adtg(S, th, P+onehalf*del_P, adtg_val)
      del_th  = del_P*adtg_val
      th      = th + (1-1/SQRT(two))*(del_th-q)
      q       = (two-sqrt(two))*del_th + (-two+three/SQRT(two))*q
!
! theta3
!
      CALL sw_adtg(S, th, P+onehalf*del_P, adtg_val)
      del_th  = del_P*adtg_val
      th      = th + (1+1/SQRT(two))*(del_th-q)
      q       = (two+SQRT(two))*del_th + (-two-three/SQRT(two))*q
!
! theta4
!
      CALL sw_adtg(S, th, P+del_P, adtg_val)
      del_th  = del_P*adtg_val
      rv      = th + (del_th-two*q)/(two*three)
!
      END SUBROUTINE sw_ptmp
!
      SUBROUTINE sw_adtg  (S, T, P, rv)
!
! !DESCRIPTION: \bv
! *=============================================================*
! | S/R  SW_ADTG
! | o compute adiabatic temperature gradient as per UNESCO 1983 routines.
! *=============================================================*
! \ev
!
! started:
!          Armin Koehl akoehl@ucsd.edu
!
! !USES:
!
!
!
! === Global variables ===
!
! !INPUT/OUTPUT PARAMETERS:
! === Routine arguments ===
!
      real :: S, T, P
      real :: rv
!
! !LOCAL VARIABLES:
! === local variables ===
!
      real :: a0, a1, a2, a3, b0, b1, c0, c1, c2, c3
      real :: d0, d1, e0, e1, e2
      real :: sref
!
      sref = 35.0
      a0 =  3.5803e-05
      a1 =  8.5258e-06
      a2 =  -6.836e-08
      a3 =  6.6228e-10
!
      b0 =  1.8932e-06
      b1 = -4.2393e-08
!
      c0 =  1.8741e-08
      c1 = -6.7795e-10
      c2 =   8.733e-12
      c3 = -5.4481e-14
!
      d0 = -1.1351e-10
      d1 =  2.7759e-12
!
      e0 = -4.6206e-13
      e1 =  1.8676e-14
      e2 = -2.1687e-16
!
      rv =  a0 + (a1 + (a2+a3*T)*T)*T +
     &      (b0+b1*T)*(S-sref) +
     &      ( (c0 + (c1 + (c2+c3*T)*T)*T) + (d0+d1*T)*(S-sref) )*P +
     &      (  e0 + (e1+e2*T)*T )*P*P
!
      END SUBROUTINE sw_adtg
!
!=======================================================================
!                                                                      !
! Use this function to do 1D interpolation.                            !
!                                                                      !
!=======================================================================
!
      SUBROUTINE linint(nx, xtab, ytab, x, y)
!
! Given a value of x return a value of y based on interpolation
! within a table of y values (ytab) corresponding to the x values
! contained in the array xtab. The subroutine assumes that the
! values in xtab increase monotonically
!
! John Mahaffy 2/12/95
! Modified by CW 2020
!
      integer, intent(in)           :: nx
      double precision, intent(in)  :: xtab(nx), ytab(nx)
      double precision, intent(in)  :: x
      double precision, intent(out) :: y
!
! local variables
!
      integer :: i
      double precision :: wx
!
      IF (x .LE. xtab(1)) THEN
        y = ytab(1)
      ELSEIF (x .GE. xtab(nx)) THEN
        y = ytab(nx)
      ELSE
        DO i = 2, nx
          IF (x .LE. xtab(i)) THEN
            wx = (x-xtab(i-1))/(xtab(i)-xtab(i-1))
            y  = (1-wx)*ytab(i-1)+wx*ytab(i)
            EXIT
          ENDIF
        ENDDO
      ENDIF
!
      END SUBROUTINE linint
#else
      subroutine iceplume_empty
      end
#endif /* ICEPLUME */
