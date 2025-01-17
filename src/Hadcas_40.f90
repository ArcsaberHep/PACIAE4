!! Hadcas_40.f90 is a part of the PACIAE event generator.
!! Copyright (C) 2024 PACIAE Group.
!! PACIAE is licensed under the GNU GPL v2 or later, see LICENSE for details.
!! Open source: https://github.com/ArcsaberHep/PACIAE4
!! Author: Ben-Hao Sa, September 2000 - January 2025.

!> This is the program to deal with the hadron cascade (hadronic rescattering).

!!                                             By Ben-Hao at CIAE on 20/09/2000
!!                                  Last updated by An-Ke at UiO  on 17/01/2025


        subroutine hadcas( time_had, ijkk )
!!      Deals with the hadronic rescattering.
!
!       It was composed by Ben-Hao Sa on 20/09/2000.
!       Its input messages are in 'sa1_h', which plays
!        the working block as well.
!       Its output messages are in 'sa1_h'.
!       If ijkk=1, gives up current event avoiding infinite loop.
!       If ijkk=2, no collisions will happen, return.
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        PARAMETER (KSZJ=80000,NSIZE=750000)
        common/sa1/kjp21,nonsa1,bp,iii,neve,nout,nosc
        common/sa1_h/nsa,non1,ksa(kszj,5),psa(kszj,5),vsa(kszj,5)
        common/sa8_h/tau(kszj),ishp(kszj)
        common/sa9_h/kfmax,kfaco(100),numb(100),non9,disbe(100,100)
        common/ctllist_h/nctl,noinel(600),nctl0,noel
        common/ctllist_t/ lc(nsize,5), tc(nsize), tw(nsize)
        common/sa24/adj1(40),nnstop,non24,zstop
!       ifram = 0 for fixed target, = 1 for collider
!       cspipi (fm^2): total cross section of pion + pion
!       cspin (fm^2): total cross section of pion + nucleon
!       cskn (fm^2): total cross section of Kaon + nucleon
!       csnn (fm^2): total cross section of n + n
!       cspsn: total cross section of J/psi (psi') + n
!       cspsm: total cross section of J/psi (psi') + meson
!       rcsit: ratio of inelastic to total cross section
!       kfmax: the maximum # of particles with given flavor code
!       kfaco(i): flavor code of i-th particle among kfmax
!       numb(i): order # of last particle of particles with same flavor of
!        kfaco(i) in particle list
!       disbe(i,j): allowable minimum approaching distance between particles
!                   kfaco(i) & kfaco(j)
!       cspipiKK (fm^2): cross section of pion + pion to Kaon + Kaon
!       edipi: largest interaction distance between two pions.
!       epin: largest interaction distance between pion and nucleon.
!       ekn: largest interaction distance between Kaon and nucleon.
!       ecsnn: largest interaction distance between two nucleons.
!       t0: average proper formation time at rest.
!       ddt: time accuracy
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!       0. is hard core distance between two pions
!       0.5 is hard core distance between two nucleons
!       0. is hard core distance between pion and nucleon
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!       tau(i) : formation time of particle i.
!       ishp(i)=1 if i-th particle inside the simulated volume
!              =0 if i-th particle outside the simulated volume
!       isinel(i) = 0 without i-th inelastic process
!                 = 1 with i-th inelastic process
!       lc(i,1) and lc(i,2) are, respectively, line # in particle
!        list of the colliding particles of i-th collision pair
!       lc(i,3) and lc(i,4) are the flavor codes of scattered particles
!        of i-th collision pair
!       lc(i,5) identifies the different inelastic process,
!        lc(i,5)=592, not used
!       tc(i): collision time of i-th colliding pair
!       tw(i): cross section ratio of (i-th inelas.)/tot
!       pio : 3.141592653589793D0
!       nctl: number of collision pairs in the current collision list
!       nctl0: number of collision pairs in last collision list
!       noinel(i): statistics of the occurring of i-th inelastic process
!       noel: statistics of the occurring of elastic process


        time = time_had
        ijkk = 0

!       Initialization.
        call sysini_h
        nctl = 0
        lc = 0
        ! tc = 0D0
        ! tw = 0D0
        numb = 0
        ! Assume the time origin is 0.
        do i=1,nsa,1
            vsa(i,4) = 0D0
            tau(i)   = 0D0
        end do
        noel = 0
        noinel = 0


!-------------------------------------------------------------------------------
!-----------------------   Interaction Volume Checking   -----------------------
        dpmax = adj1(27)
        drmax = adj1(28)
        do i1=1,nsa,1
            ishp(i1) = 0
            pnn1 = psa(i1,1)
            pnn2 = psa(i1,2)
            pnn3 = psa(i1,3)
            pnn4 = psa(i1,4)
            rnn1 = vsa(i1,1)
            rnn2 = vsa(i1,2)
            rnn3 = vsa(i1,3)
            pnnm = pnn1*pnn1 + pnn2*pnn2 + pnn3*pnn3
            if( pnnm <= 1D-28 ) pnnm = 1D-28
            if( pnnm > 1D28 ) cycle
            pnnm = SQRT(pnnm)
            rnnm = rnn1*rnn1 + rnn2*rnn2 + rnn3*rnn3
            if( rnnm <= 1D-28 ) rnnm = 1D-28
            if( rnnm > 1D28  ) cycle
            rnnm = SQRT(rnnm)
            if( pnnm <= dpmax .AND. pnn4 <= dpmax &
                .AND. rnnm <= drmax )then
                ishp(i1) = 1
            end if
        end do
!-----------------------   Interaction Volume Checking   -----------------------
!-------------------------------------------------------------------------------


!       Changes K0S, K0L to K0, K0bar.
        do j=1,nsa,1
            kf = ksa(j,2)
            if(kf == 130 .OR. kf == 310)then
                rrlu = PYR(1)
                ksa(j,2) = 311
                if(rrlu > 0.5D0) ksa(j,2) = -311
            end if
        end do

!       Filters out particles wanted to study and make in order of proton,
!        neutron, ...
!       Initial particle list is comprised of the arrays in common block
!        'sa1_h', 'tau' and 'ishp' in 'sa8_h', and 'numb' in 'sa9_h'
        call filt_h

!       Calculates position of center of mass of the system. distance of a
!        particle from this cms is used to check whether it is freezes out
!        or not
        call copl_h(time)

!       Creates the initial collision list
        call ctlcre_h(time)
!       No collisions will happen.
        if( nctl <= 0 )then
            ijkk = 2
            return
        end if

!       Administrates the hadron rescattering.
        call scat_h( time, ijkk, iii )
        if(ijkk == 1) return
!       If ijkk=1 (infinite loops), give up the event.
        time_had = time

!       Changes K0, K0bar to K0L and K0S.
        do j=1,nsa,1
            kf = ksa(j,2)
            if(kf == 311 .or. kf == -311)then
                rrlu = pyr(1)
                ksa(j,2) = 130
                if(rrlu > 0.5D0) ksa(j,2) = 310
            end if
        end do


        return
        end



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine sysini_h
!!      Gives initial values to the quantities needed in hadron rescattering.
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        COMMON/PYCIDAT1/KFACOT(100),DISDET(100),ISINELT(600)
        COMMON/PYCIDAT2/KFMAXT,NONT2,PARAM(20),WEIGH(600)
        common/sa9_h/kfmax,kfaco(100),numb(100),non9,disbe(100,100)
        common/sa10/csnn_p,cspin_p,cskn_p,cspipi_p,cspsn_p,cspsm_p, &
         rcsit_p,ifram_p,iabsb_p,iabsm_p,i_sigma_AQM_p,ajpsi_p,csspn_p, &
         csspm_p,csen_p
        common/sa10_h/csnn,cspin,cskn,cspipi,cspsn,cspsm,rcsit,ifram, &
         iabsb,iabsm,i_sigma_AQM,csspn,csspm
        common/sa20_h/t0,cspipiKK,dep,ddt,edipi,epin,ecsnn,ekn,ecspsn,ecspsm, &
         rnt,rnp,ecsspn,ecsspm
        common/sa25/i_inel_proc,i_time_shower,para1_1,para1_2
        common/sa24/adj1(40),nnstop,non24,zstop
        common/syspar_h/pio
        common/count_h/isinel(600)
!       For cross sections of hh collisions etc.
        common/para_h1/ para(20)


        pio = 3.141592653589793D0
        PARAM(1) = para1_2
        PARAM(5) = para(5)

!       cross sections are given in fm^2
        csnn     = PARAM(1)*0.1D0
        cspin    = PARAM(2)*0.1D0
        cskn     = PARAM(3)*0.1D0
        cspipi   = PARAM(4)*0.1D0
        cspipiKK = PARAM(5)*0.1D0
        cspsn    = PARAM(13)*0.1D0
        cspsm    = PARAM(14)*0.1D0
        csspn    = PARAM(15)*0.1D0
        csspm    = PARAM(16)*0.1D0
        iabsb    = 1
        iabsm    = 1
!       PARAM(1) : total cross-section of nucleon + nucleon (in mb)
!       PARAM(5) : cross-section of pi + pi -> K + Kbar (in mb)
!       csnn  : total cross section of nucleon + nucleon
!       1 mb = 10 fm^2
!       cspin : total cross section of pion + nucleon
!       cskn  : total cross section of Kaon + nucleon
!       cspipi: total cross section of pion + pion
!       cspsn : total cross section of J/psi (psi') + nucleon
!       cspsm : total cross section of J/psi (psi') + meson
!       iabsb = 0 : without J/psi (psi') + baryon (nucleon)
!             = 1 : with J/psi (psi')    + baryon (nucleon)
!       iabsm = 0 : without J/psi (psi') + meson
!             = 1 : with J/psi (psi')    + meson

!       i_sigma_AQM: = 1, uses the input total cross sections of piN, KN, pi+pi,
!                         Jpsi(psi')+N, Jpsi(psi')+N, Jpsi(psi')+pi/rho,
!                         and AQM cross sections for D+N/pi/rho in hadcas.
!                         Ignores the channels which are not well defined.
!                  : = 2, uses the input cross sections of piN, KN, ... ,
!                         AQM cross sections for D+N/pi/rho and the channels
!                         which are not well defined (treated as elastic).
!                  : = 3, uses AQM cross sections of piN, KN,... and D+N/pi/rho.
!                         Ignores the channels which are not well defined.
!                  : = 4, uses AQM cross sections of piN, KN, ..., D+N/pi/rho,
!                         and the channels which are not well defined (treated
!                         as elastic).

!       Uses the total cross sections from AQM (additive quark model,
!        arXiv:2203.11601)
        i_sigma_AQM = i_sigma_AQM_p
        if( i_sigma_AQM == 3 .OR. i_sigma_AQM == 4 )then
            sigma_NN_tot = para1_2
            cspin  = sigma_AQM( 211, 2212, sigma_NN_tot ) * 0.1D0
            cskn   = sigma_AQM( 321, 2212, sigma_NN_tot ) * 0.1D0
            cspipi = sigma_AQM( 211,  211, sigma_NN_tot ) * 0.1D0
            cspsn  = sigma_AQM( 443, 2212, sigma_NN_tot ) * 0.1D0
            cspsm  = sigma_AQM( 443,  221, sigma_NN_tot ) * 0.1D0
            csspn  = cspsn
            csspm  = cspsm
        end if

!       largest interaction distance of two colliding particles.
        edipi  = sqrt(cspipi / pio)
        epin   = sqrt(cspin  / pio)
        ekn    = sqrt(cskn   / pio)
        ecsnn  = sqrt(csnn   / pio)
        ecspsn = sqrt(cspsn  / pio)
        ecspsm = sqrt(cspsm  / pio)
        ecsspn = sqrt(csspn  / pio)
        ecsspm = sqrt(csspm  / pio)

        rcsit = PARAM(6)
        t0    = PARAM(7)
!       dep   = PARAM(9)
!       ddt   = PARAM(8)
        ddt   = adj1(11)
!       rao   = PARAM(10)
        kfmax  = KFMAXT
        kfaco  = KFACOT
        isinel = ISINELT
        disbe  = 0D0
        do j=1,kfmax,1
            disbe(1,j)  = DISDET(j)
            disbe(2,j)  = DISDET(j)
            disbe(3,j)  = DISDET(j)
            disbe(4,j)  = DISDET(j)
!           disbe(26,j) = DISDET(j)
!           disbe(27,j) = DISDET(j)
!           disbe(28,j) = DISDET(j)
!           disbe(29,j) = DISDET(j)
        end do
        do i=1,99,1
            do j = i+1, 100, 1
                disbe(j,i) = disbe(i,j)
            end do
        end do


        return
        end



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine filt_h
!!      Filters out desired particles and make them in order of
!               1            2          3          4            5             6             7             8         9           10
!    0     proton,     neutron,     pbar-,      nbar,         pi+,          pi-,          pi0,           K-,    Kbar0,      Sigma0,
!    1     Sigma-,      Sigma+, Sigmabar0, Sigmabar+,   Sigmabar-,      Lambda0,   Lambdabar0,           K0,       K+,         Xi-,
!    2     Xibar+,         Xi0,    Xibar0,    Omega-,   Omegabar+,       Delta-,       Delta0,       Delta+,  Delta++,        rho+,
!    3       rho-,        rho0,     J/psi,      psi',   Deltabar+,    Deltabar0,    Deltabar-,   Deltabar--,       D+,          D-,
!    4         D0,       Dbar0,       D*+,       D*-,         D*0,       D*bar0,    Lambda_c+, Lambda_cbar-,     D_s+,        D_s-,
!    5      D*_s+,       D*_s-,       K*+,       K*-,         K*0,       K*bar0,      Upsilon,     Upsilon',   chi_0c,      chi_1c,
!    6     chi_2c,    Sigma_c0,  Sigma_c+, Sigma_c++, Sigma_cbar0,  Sigma_cbar+, Sigma_cbar++,        omega,       B0,       B0bar,
!    7         B+,          B-,      B_s0,   B_sbar0,        B_c+,         B_c-,          B*0,       B*bar0,      B*+,         B*-,
!    8      B*_s0,    B*_sbar0,     B*_c+,     B*_c-,   Lambda_b0, Lambda_bbar0,     Sigma_b0,  Sigma_bbar0, Sigma_b-, Sigma_bbar+,
!    9   Sigma_b+, Sigma_bbar-,        8*0
!       (92 particle)
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        PARAMETER (KSZJ=80000)
        common/sa1_h/nsa,non1,ksa(kszj,5),psa(kszj,5),vsa(kszj,5)
        common/sa9_h/kfmax,kfaco(100),numb(100),non9,disbe(100,100)


        iiii = 0
        jjjj = 0
        do i=1,kfmax,1
            kf = kfaco(i)
            do j = iiii+1, nsa, 1
                call ord_h( jjjj, j, kf )
            end do
            iiii = jjjj
            numb(i) = jjjj
        end do


        return
        end



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine ord_h(ipi,j,kf)
!!      Makes in order for particles with flavor code kf.
!       j : the particle wanted to order
!       kf: flavor code of j-th particle
!       ipi : j-th particle should order after ipi
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        PARAMETER (KSZJ=80000)
        common/sa1_h/nsa,non1,ksa(kszj,5),psa(kszj,5),vsa(kszj,5)
        dimension kk(5),pp(5),vv(5)


        ik = ksa(j,2)
        if(ik == kf)then
            ipi = ipi + 1
            do jj=1,5,1
                kk(jj) = ksa(ipi,jj)
                pp(jj) = psa(ipi,jj)
                vv(jj) = vsa(ipi,jj)
            end do
            do jj=1,5,1
                ksa(ipi,jj) = ksa(j,jj)
                psa(ipi,jj) = psa(j,jj)
                vsa(ipi,jj) = vsa(j,jj)
            end do
            do jj=1,5,1
                ksa(j,jj) = kk(jj)
                psa(j,jj) = pp(jj)
                vsa(j,jj) = vv(jj)
            end do
        end if


        return
        end



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine copl_h(time)
!!      Calculate position of center of mass of the non-freeze-out system
!!      distance of a particle from this cms is used to checke whether
!!       it freezes out or not.
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        PARAMETER (KSZJ=80000)
        common/sa1_h/nsa,non1,ksa(kszj,5),psa(kszj,5),vsa(kszj,5)
        common/sa8_h/tau(kszj),ishp(kszj)
        common/sa19_h/coor(3)


        tt = time
        coor = 0D0
        samass = 0D0
        do ii=1,nsa,1
!           if( tau(ii) > tt ) cycle
            if( ishp(ii) == 0 ) cycle
            KF = ksa(ii,2)
            damass  = PYMASS( KF )
            samass = samass + damass
            do jj=1,3,1
                coor(jj) = coor(jj) + damass * vsa(ii,jj)
            end do
        end do
        ! 0.14 : pion mass.
        coor = coor / MAX( 0.14D0, samass )


        return
        end



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine ctlcre_h(time)
!!      Creates the initial collision list.
!       The time origin is assumed as 0.
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        PARAMETER (KSZJ=80000,NSIZE=750000)
        common/sa1_h/nsa,non1,ksa(kszj,5),psa(kszj,5),vsa(kszj,5)
        common/sa9_h/kfmax,kfaco(100),numb(100),non9,disbe(100,100)
        common/sa20_h/t0,cspipiKK,dep,ddt,edipi,epin,ecsnn,ekn,ecspsn,ecspsm, &
         rnt,rnp,ecsspn,ecsspm
        common/ctllist_h/nctl,noinel(600),nctl0,noel
        common/ctllist_t/ lc(nsize,5), tc(nsize), tw(nsize)


        time_resolution = ddt
        n_max_considered = numb(kfmax)
        nctl = 1
        do i = 1, nsa-1, 1
            if( i > n_max_considered ) cycle
            do j = i+1, nsa, 1
                if( j > n_max_considered ) cycle
                iflag = 0
                call rsfilt_h( i, j, iflag )
                if( iflag == 0 ) cycle
                tc(nctl) = 0D0
                tw(nctl) = 0D0
                call tcolij_h( i, j, time, nctl )
                tci = tc(nctl)
                if( tci <= 1D-7 ) cycle
                do jj = 1, nctl-1, 1
                    tcj = tc( jj )
                    if( ABS(tcj - tci) < time_resolution )then
                        nctl = nctl - 1
                        exit
                    end if
                end do
                nctl = nctl + 1
            end do
        end do
        nctl = nctl - 1


        return
        end



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine scat_h( time, ijkk, iii )
!!      Administrates the hadronic rescattering.
!       iii: the run number
!       if ijkk=1 (infinite loop) give up the event
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        PARAMETER (KSZJ=80000,NSIZE=750000)
        common/sa1_h/nsa,non1,ksa(kszj,5),psa(kszj,5),vsa(kszj,5)
        common/sa8_h/tau(kszj),ishp(kszj)
        common/sa9_h/kfmax,kfaco(100),numb(100),non9,disbe(100,100)
        common/sa10_h/csnn,cspin,cskn,cspipi,cspsn,cspsm,rcsit,ifram, &
         iabsb,iabsm,i_sigma_AQM,csspn,csspm
        common/sa19_h/coor(3)
        common/sa20_h/t0,cspipiKK,dep,ddt,edipi,epin,ecsnn,ekn,ecspsn,ecspsm, &
         rnt,rnp,ecsspn,ecsspm
        common/ctllist_h/nctl,noinel(600),nctl0,noel
        common/ctllist_t/ lc(nsize,5), tc(nsize), tw(nsize)
        dimension pi(4),pj(4),pii(4),pjj(4),b(3)
        integer winel


        ijkk  = 0
        nctl0 = nctl
        iiii  = 1

10      if(iiii > 1) call copl_h(time)

!       Finds out the binary colli. with minimum collsion time.
        call find_h( icp, tcp )
        if(icp == 0) return
!       icp=0 means the collision list is empty.
        if( tcp > 1D10 )then
            do i1 = icp+1, nctl, 1
                do j1=1,5,1
                    lc( i1-1, j1 ) = lc(i1,j1)
                end do
                tc( i1-1 ) = tc(i1)
                tw( i1-1 ) = tw(i1)
            enddo
            nctl = nctl - 1
            goto 10
        end if
        time0 = time
        l     = lc(icp,1)
        l1    = lc(icp,2)
        kfa   = ksa(l,2)
        kfb   = ksa(l1,2)
!       Records this collision time.
        time  = tcp

!       Performs classical Newton motion in Lab. system.
        call his_h( time, istop )
        if(istop == 1) return
!       istop=1 means all particles have get out of considered volume.

        pi(4) = psa(l,4)
        pj(4) = psa(l1,4)
        do i=1,3,1
            pi(i) = psa(l,i)
            pj(i) = psa(l1,i)
            b(i)  = ( pi(i) + pj(i) ) / ( pi(4) + pj(4) )
        end do
!       Boosts to CMS frame of the colliding pair.
        ilo = 0
        call lorntz( ilo, b, pi, pj )
        ss = pi(4) + pj(4)
        ww = rcsit
!       The cross section ratio of (ela.)/tot = 1 - rcsit.
        rrlu = PYR(1)
        winel = 0
        ! inela.
        if( rrlu <= ww )then
            winel = 1
            call coinel( l,l1, ss, b, pi,pj, icp, pii,pjj, winel, ik3,ik4 )
!       If winel=0 the inelastic reaction has not really happened.
!       If winel=1 inelastic collision happened, ik3 and ik4 are
!        flavors of the scattered particles, pii and pjj are four-momentum of
!        scattered particles in Lab frame, the two colliding particles are
!        still with line numbers of l and l1 in the particle list
        end if

        if( winel /= 0 )then
!           Treats the inelastic collision.
            icp5 = lc(icp,5)
            noinel(icp5) = noinel(icp5) + 1
!           Updates particle list after inelastic collision
!            and truncates collision list correspondingly.
            call updpli( l,l1, icp, pii, pjj, time, icp5 )
            l  = lc(icp,1)
            l1 = lc(icp,2)
!           l and l1 are now the line numbers of scattered particles in
!            particle list.
        else
!           Calculates four-momentum of scattered particles after elastic
!            collistion (pi and pj in CMS frame).
            call coelas_h(l,l1,ss,pi,pj)
!           Updates particle list for elastic scattering, pi and pj have been
!            boosted back to the Lab frame.
            call updple_h(l,l1,b,pi,pj)
            noel = noel + 1
        end if

!       Updates the collision list.
        call updatl_h(l,l1,time)

!#TODO(Lei20241012): particle decay in the transport process.
!       if( nctl <= 1 ) goto 300
!       deltt = time - time0
!       deal with particle decay in the transport processes
!       call decay( time, deltt )

        iiii = iiii + 1

        if( iiii > 100*(nctl0) )then
            write(9,*) "infinite loop may have happened in" &
                    // "subroutine scat iii=", iii
            iii = iii - 1
!           If ijkk=1 (infinite loop), give up the event.
            ijkk = 1
            return
        end if

        goto 10


        return
        end



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine find_h(icp,tcp)
!!      Finds out the binary collision with minimum collision time.
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        PARAMETER(NSIZE=750000)
        common/ctllist_h/nctl,noinel(600),nctl0,noel
        common/ctllist_t/ lc(nsize,5), tc(nsize), tw(nsize)


        icp = 0
        tcp = 20000D0
        do i=1,nctl,1
            if( tc(i) <= 1D-7 ) cycle
            if( tcp < tc(i) ) cycle
            icp = i
            tcp = tc(i)
        end do


        return
        end



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine his_h(t1,istop)
!!      Classical Newton motion in Lab. system.
!       istop=1 means all particles have get out of considered volume.
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        PARAMETER (KSZJ=80000,NSIZE=750000)
        common/sa1_h/nsa,non1,ksa(kszj,5),psa(kszj,5),vsa(kszj,5)
        common/sa8_h/tau(kszj),ishp(kszj)
        common/sa19_h/coor(3)
        common/sa20_h/t0,cspipiKK,dep,ddt,edipi,epin,ecsnn,ekn,ecspsn,ecspsm, &
         rnt,rnp,ecsspn,ecsspm
        common/sa24/adj1(40),nnstop,non24,zstop
        common/ctllist_h/nctl,noinel(600),nctl0,noel
        common/ctllist_t/ lc(nsize,5), tc(nsize), tw(nsize)


        istop = 0
        i_n   = 0
        r0    = adj1(28)
        do i=1,nsa,1
!       do move particles which have not been produced
!           if( t1 <= tau(i) ) cycle
            if( ishp(i) == 0 )then
                i_n = i_n + 1
                cycle
            end if
            aa = 0D0
            pp4 = psa(i,4)
            do j=1,3,1
                vp = psa(i,j) / pp4
                vsa(i,j) = vsa(i,j) + vp * ( t1 - vsa(i,4) )
                aa = aa + ( vsa(i,j) - coor(j) )**2
            end do
            aa = SQRT(aa)
            if(aa < r0)then
                vsa(i,4) = t1
            else
!       if freeze-out already, deduct the distance between the last and
!        current collisions
                do j=1,3,1
                    vp = psa(i,j) / pp4
                    vsa(i,j) = vsa(i,j) - vp * ( t1 - vsa(i,4) )
                end do
                ishp(i) = 0
                do il=1,nctl,1
                    if( lc(il,1) == i .OR. lc(il,2) == i ) tc(il) = 0D0
                end do
            end if
        end do
        if( i_n == nsa ) istop = 1


        return
        end



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine coinel( l,l1, ss,b, pi,pj, icp, pii,pjj, winel, ik1,ik2 )
!!      Treats the inelastic collision
!
!       The inelastic processes considered are:

!       strangeness production reactions:
!        1. pion+ + pion- to k+ + k-
!        2. pion+ + pion- to k0 + k0-
!        3. pion+ + pion0 to k+ + k0-
!        4. pion- + pion0 to k- + k0
!        5. pion0 + pion0 to k+ + k-
!        6. pion0 + pion0 to k0 + k0-

!        7. pion+ + p to k+ + sigma+
!        8. pion+ + n to k+ + sigma0
!        9. pion+ + n to k+ + lambda
!        10. pion+ + n to k0 + sigma+
!        11. pion- + p to k+ + sigma-
!        12. pion- + p to k0 + lambda
!        13. pion- + p to k0 + sigma0
!        14. pion- + n to k0 + sigma-
!        15. pion0 + p to k+ + sigma0
!        16. pion0 + p to k+ + lambda
!        17. pion0 + p to k0 + sigma+
!        18. pion0 + n to k+ + sigma-
!        19. pion0 + n to k0 + lambda
!        20. pion0 + n to k0 + sigma0

!        21. pion+ + pba to k0- + lambdaba
!        22. pion+ + pba to k0- + sigma0ba
!        23. pion+ + pba to k- + sigma-ba
!        24. pion+ + nba to k0- + sigma-ba
!        25. pion- + pba to k- + sigma+ba
!        26. pion- + nba to k- + lambdaba
!        27. pion- + nba to k- + sigma0ba
!        28. pion- + nba to k0- + sigma+ba
!        29. pion0 + pba to k- + lambdaba
!        30. pion0 + pba to k- + sigma0ba
!        31. pion0 + pba to k0- + sigma+ba
!        32. pion0 + nba to k- + sigma-ba
!        33. pion0 + nba to k0- + lambdaba
!        34. pion0 + nba to k0- + sigma0ba

!        35. pion+ + sigma- to k+ + cascade-
!        36. pion- + lambda to k0 + cascade-
!        37. pion- + sigma+ to k+ + cascade-
!        38. pion- + sigma0 to k0 + cascade-
!        39. pion0 + lambda to k+ + cascade-
!        40. pion0 + sigma- to k0 + cascade-
!        41. pion0 + sigma0 to k+ + cascade-
!        42. pion+ + lambdaba to k0- + cascade-ba
!        43. pion+ + sigma+ba to k- + cascade-ba
!        44. pion+ + sigma0ba to k0- + cascade-ba
!        45. pion- + sigma-ba to k- + cascade-ba
!        46. pion0 + lambdaba to k- + cascade-ba
!        47. pion0 + sigma-ba to k0- + cascade-ba
!        48. pion0 + sigma0ba to k- + cascade-ba

!        strangeness exchange reactions:
!        49. k- + p to pion0 + lambda
!        50. k- + p to pion0 + sigma0
!        51. k- + p to pion- + sigma+
!        52. k- + p to pion+ + sigma-
!        53. k- + p to k+ + cascade- (strangeness production)
!        54. k- + n to pion- + sigma0
!        55. k- + n to pion- + lambda
!        56. k- + n to pion0 + sigma-
!        57. k- + n to k0 + cascade- (strangeness production)

!        58. k0- + p to pion+ + lambda
!        59. k0- + p to pion+ + sigma0
!        60. k0- + p to pion0 + sigma+
!        61. k0- + n to pion+ + sigma-
!        62. k0- + n to pion- + sigma+
!        63. k0- + n to pion0 + sigma0
!        64. k0- + n to pion0 + lambda
!        65. k0- + n to k+ + cascade- (strangeness production)

!        66. k+ + p- to pion0 + lambda-
!        67. k+ + p- to pion0 + sigma0-
!        68. k+ + p- to pion+ + sigma+ba
!        69. k+ + p- to pion- + sigma-ba
!        70. k+ + p- to k- + cascade-ba (strangeness production)

!        71. k+ + n- to pion+ + sigma0-
!        72. k+ + n- to pion+ + lambda-
!        73. k+ + n- to pion0 + sigma-ba
!        74. k+ + n- to k0- + cascade-ba(strangeness production)

!        75. k0 + p- to pion- + lambda-
!        76. k0 + p- to pion- + sigma0-
!        77. k0 + p- to pion0 + sigma+ba
!        78. k0 + n- to pion- + sigma-ba
!        79. k0 + n- to pion+ + sigma+ba
!        80. k0 + n- to pion0 + sigma0-
!        81. k0 + n- to pion0 + lambda-
!        82. k0 + n- to k- + cascade-ba (strangeness production)

!        83. k- + lambda to pion0 + cascade-
!        84. k- + sigma+ to pion+ + cascade-
!        85. k- + sigma- to pion- + cascade-
!        86. k- + sigma0 to pion0 + cascade-
!        87. k0- + lambda to pion+ + cascade-
!        88. k0- + sigma0 to pion+ + cascade-
!        89. k0- + sigma- to pion0 + cascade-
!        90. k+ + lambda- to pion0 + cascade-ba
!        91. k+ + sigma+ba to pion- + cascade-ba
!        92. k+ + sigma-ba to pion+ + cascade-ba
!        93. k+ + sigma0- to pion0 + cascade-ba
!        94. k0 + lambda- to pion- + cascade-ba
!        95. k0 + sigma0- to pion- + cascade-ba
!        96. k0 + sigma-ba to pion0 + cascade-ba
!        97. pion+ + sigma- to k0 + cascade0
!        98. pion+ + sigma0 to k+ + cascade0
!        99. pion+ + lambda0 to k+ + cascade0
!       100. pion- + sigma+ to k0 + cascade0
!       101. pion0 + sigma+ to k+ + cascade0
!       102. pion0 + sigma0 to k0 + cascade0
!       103. pion0 + lambda to k0 + cascade0
!       104. pion+ + sigma+ba to k0- + cascade0-
!       105. pion- + sigma-ba to k0- + cascade0-
!       106. pion- + sigma0- to k- + cascade0-
!       107. pion- + lambda- to k- + cascade0-
!       108. pion0 + sigma+- to k- + cascade0-
!       109. pion0 + sigma0- to k0- + cascade0-
!       110. pion0 + lambda- to k0- + cascade0-
!       111. k- + sigma+ to pion0 + cascade0
!       112. k- + sigma0 to pion0 + cascade-
!       113. k- + lambda to pion- + cascade0
!       114. k0- + sigma+ to pion+ + cascade0
!       115. k0- + sigma- to pion- + cascade0
!       116. k0- + sigma0 to pion0 + cascade0
!       117. k0- + lambda to pion0 + cascade0
!       118. k+ + sigma+ba to pion0 + cascade0-
!       119. k+ + sigma0- to pion+ + cascade0-
!       120. k+ + lambda- to pion+ + cascade0-
!       121. k+ + cascade-ba to pion+ + Omega-ba
!       122. k0 + sigma-ba to pion+ + cascade0ba
!       123. k0 + sigma0- to pion0 + cascade0-
!       124. k0 + lambda- to pion0 + cascade0ba
!       125. k- + p to k0 + cascade0
!       126. k0- + p to k+ + cascade0
!       127. k0- + n to k0 + cascade0
!       128. k+ + p- to k0- + cascade0ba
!       129. k0 + p- to k- + cascade0-
!       130. k0 + n- to k0- + cascade0-
!       131. pion+ + cascade- to k+ + Omega-
!       132. pion0 + cascade- to k0 + Omega-
!       133. pion- + cascade-ba to k- + Omega-ba
!       134. pion0 + cascade-ba to k0- + Omega-ba
!       135. pion- + cascade0 to k0 + Omega-
!       136. pion0 + cascade0 to k+ + Omega-
!       137. pion+ + cascade0- to k0- + Omega-ba
!       138. pion0 + cascade0- to k- + Omega-ba
!       139. k- + cascade- to pion- + Omega-
!       140. k0- + cascade- to pion0 + Omega-
!       141. k- + cascade0 to pion0 + Omega-
!       142. k0- + cascade0 to pion+ + Omega-
!       143. k+ + cascade-ba to pion+ + Omega-ba
!       144. k0 + cascade-ba to pion0 + Omega-ba
!       145. k+ + cascade0- to pion0 + Omega-ba
!       146. k0 + cascade0- to pion- + Omega-ba

!       147. pion- + p to delta- + pion+
!       148. pion- + p to rho0 + n
!       149. pion- + p to rho- + p
!       150. pion- + p to delta+ + pion-
!       151. pion- + p to delta0 + pion0
!       152. pion- + n to delta- + pion0
!       153. pion- + n to rho- + n
!       154. pion- + n to delta0 + pion-
!       155. pion+ + p to delta++ + pion0
!       156. pion+ + p to delta+ + pion+
!       157. pion+ + p to rho+ + p
!       158. pion+ + n to delta++ + pion-
!       159. pion+ + n to delta0 + pion+
!       160. pion+ + n to delta+ + pion0
!       161. pion+ + n to rho0 + p
!       162. pion+ + n to rho+ + n
!       163. pion0 + p to delta0 + pion+
!       164. pion0 + p to delta++ + pion-
!       165. pion0 + p to rho+ + n
!       166. pion0 + p to rho0 + p
!       167. pion0 + p to delta+ + pion0
!       168. pion0 + n to delta+ + pion-
!       169. pion0 + n to delta- + pion+
!       170. pion0 + n to delta0 + pion0
!       171. pion0 + n to rho0 + n
!       172. pion0 + n to rho- + p
!       173. p + p to delta+ + p
!       174. p + p to delta++ + n
!       175. p + n to delta+ + n
!       176. p + n to delta0 + p
!       177. n + n to delta0 + n
!       178. n + n to delta- + p

!       179. J/psi + n to lamdac + Dba
!       180. J/psi + n to sigmac + Dba
!       181. J/psi + n to sigmac0 + D0ba
!       182. J/psi + p to lamdac + D0ba
!       183. J/psi + p to sigmac + D0ba
!       184. J/psi + p to sigmac++ + Dba
!       185. J/psi + pion+ to D + D*0ba
!       186. J/psi + pion0 to D0 + D*0ba
!       187. J/psi + pion0 to D + D*ba
!       188. J/psi + pion- to D0 + D*ba
!       189. J/psi + rho+ to D + D0ba
!       190. J/psi + rho0 to D0 + D0ba
!       191. J/psi + rho0 to D + Dba
!       192. J/psi + rho- to D0 + Dba
!       193. psi' + n to lamdac + Dba
!       194. psi' + n to sigmac + Dba
!       195. psi' + n to sigmac0 + D0ba
!       196. psi' + p to lamdac + D0ba
!       197. psi' + p to sigmac + D0ba
!       198. psi' + p to sigmac++ + Dba
!       199. psi' + pion+ to D + D*0ba
!       200. psi' + pion0 to D0 + D*0ba

!       reverse reactions, after '201'
!       201. pion+ + pion- to k+ + k- (r)
!       202. pion0 + pion0 to k+ + k-(r)
!       203. pion+ + pion0 to k+ + k0-(r)
!       204. pion- + pion0 to k- + k0(r)
!       205. pion+ + pion- to k0 + k0-(r)
!       206. pion0 + pion0 to k0 + k0-(r)
!       207. k+ + sigma+ to pion+ + p
!       208. k+ + sigma- to pion- + p
!       209. k+ + sigma- to pion0 + n
!       210. k+ + sigma0 to pion+ + n
!       211. k+ + sigma0 to pion0 + p
!       212. k+ + lambda0 to pion+ + n
!       213. k+ + lambda0 to pion0 + p
!       214. k+ + sigma+ to pion+ + n
!       215. k+ + sigma+ to pion0 + p
!       216. k0 + sigma- to pion- +n
!       217. k0 + sigma0 to pion- + p
!       218. k0 + sigma0 to pion0 + n
!       219. k0 + lambda0 to pion- + p
!       220. k0 + lambda0 to pion0 + n
!       221. k- + sigma+ba to pion- + pba
!       222. k- + sigma-ba to pion+ + pba
!       223. k- + sigma-ba to pion0 + nba
!       224. k- + sigma0ba to pion- + nba
!       225. k- + sigma0ba to pion0 + pba
!       226. k- + lambda0ba to pion- + nba
!       227. k- + lambda0ba to pion0 + pba
!       228. k0ba + sigma+ba to pion- + nba
!       229. k0ba + sigma+ba to pion0 + pba
!       230. k0- + sigma-ba to pion+ + nba
!       231. k0- + sigma0- to pion+ + pba
!       232. k0- + sigma0- to pion0 + nba
!       233. k0- + lambda0- to pion+ + pba
!       234. k0- + lambda0- to pion0 + nba
!       235. k++ cascade- to pi+ sigma-
!       236. k+ + cascade- to pion- + sigma+
!       237. k+ + cascade- to pion0 + lambda0
!       238. k+ + cascade- to pion0 + sigma0
!       239. k- + cascade-ba to pion+ + sigma+ba
!       240. k- + cascade-ba to pion- + sigma-ba
!       241. k- + cascade-ba to pion0 + lambda0-
!       242. k- + cascade-ba to pion0 + sigma0-
!       243. k0 + cascade- to pion- + lambda0
!       244. k0 + cascade- to pion- + sigma0
!       245. k0 + cascade- to pion0 + sigma-
!       246. k0- + cascade-ba to pion+ + lambda0-
!       247. k0- + cascade-ba to pion+ + sigma0-
!       248. k0- + cascade-ba to pion0 + sigma-ba
!       249. pion+ + sigma- to k- + p
!       250. pion+ + sigma- to k0- + n
!       251. pion+ + sigma0 to k0- + p
!       252. pion+ + lambda0 to k0- + p
!       253. k+ + cascade- to k- + p
!       254. pion- + sigma+ to k- + p
!       255. pion- + sigma+ to k0- + n
!       256. pion- + sigma0 to k- + n
!       257. k0 + cascade- to k- + n
!       258. pion- + lambda0 to k- + n
!       259. pion0 + sigma+ to k0- + p
!       260. pion0 + sigma- to k- + n
!       261. pion0 + sigma0 to k- + p
!       262. pion0 + sigma0 to k0- + n
!       263. pion0 + lambda0 to k- + p
!       264. pion0 + lambda0 to k0- + n
!       265. k+ + cascade- to k0- + n
!       266. pion+ + sigma+ba to k+ + pba
!       267. pion+ + sigma+ba to k0 + nba
!       268. pion+ + sigma0- to k+ + nba
!       269. pion+ + lambda0- to k+ + nba
!       270. k- + cascade-ba to k+ + pba
!       271. pion- + sigma-ba to k+ + pba
!       272. pion- + sigma-ba to k0 + nba
!       273. pion- + sigma0- to k0 + pba
!       274. k0- + cascade-ba to k+ + nba
!       275. pion- + lambda0- to k0 + pba
!       276. pion0 + sigma+- to k0 + pba
!       277. pion0 + sigma-ba to k+ + nba
!       278. pion0 + sigma0- to k+ + pba
!       279. pion0 + sigma0- to k0 + nba
!       280. pion0 + slambda0- to k+ + pba
!       281. pion0 + lambda0- to k0 + nba
!       282. k- + cascade-ba to k0 + nba
!       283. pion+ + cascade- to k- + sigma+
!       284. pion+ + cascade- to k0- + lambda0
!       285. pion+ + cascade- to k0- + sigma0
!       286. pion- + cascade- to k- + sigma-
!       287. pion0 + cascade- to k- + lambda0
!       288. pion0 + cascade- to k- + gigma0
!       289. pion0 + cascade- to k0- + sigma-
!       290. pion+ + cascade-ba to k+ + sigma-ba
!       291. pion- + cascade-ba to k+ + sigma+-
!       292. pion- + cascade-ba to k0 + lambda0-
!       293. pion- + cascade-ba to k0 + sigma0-
!       294. pion0 + cascade-ba to k+ + lambda0-
!       295. pion0 + cascade-ba to k+ + gigma0-
!       296. pion0 + cascade-ba to k0 + sigma-ba
!       297. k0 + cascade0 to pion+ + sigma-
!       298. k+ + cascade0 to pion+ + sigma0
!       299. k+ + cascade0 to pion+ + lambda
!       300. k0 + cascade0 to pion- + sigma+
!       301. k+ + cascade0 to pion0 + sigma+
!       302. k0 + cascade0 to pion0 + sigma0
!       303. k0 + cascade0 to pion0 + lambda
!       304. k-0 + cascade0- to pion+ + sigma-ba
!       305. k0- + cascade0- to pion- + sigma-ba
!       306. k- + cascade0- to pion- + sigma0-
!       307. k- + cascade0- to pion- + lambda-
!       308. k- + cascade0- to pion0 + sigma+ba
!       309. k0- + cascade0- to pion0 + sigma0-
!       310. k0- + cascade0- to pion0 + lambda-
!       311. pion0 + cascade0 to k- + sigma+
!       312. pion- + cascade0 to k- + sigma0
!       313. pion- + cascade0 to k- + lambda
!       314. pion+ + cascade0 to k0- + sigma+
!       315. pion- + cascade0 to k0- + sigma-
!       316. pion0 + cascade0 to k0- + sigma0
!       317. pion0 + cascade0 to k0- + lambda
!       318. pion0 + cascade0- to k+ + sigma+ba
!       319. pion+ + cascade0- to k+ + sigma0-
!       320. pion+ + cascade0- to k+ + lambda-
!       321. pion- + cascade0- to k0 + sigma+ba
!       322. pion+ + cascade0- to k0 + sigma-ba
!       323. pion0 + cascade0- to k0 + sigma0-
!       324. pion0 + cascade0- to k0 + lambda-
!       325. k0 + cascade0 to k- + p
!       326. k+ + cascade0 to k0- + p
!       327. k0 + cascade0 to k0- + n
!       328. k0- + cascade0- to k+ + p-
!       329. k- + cascade0- to k0 + p-
!       330. k0- + cascade0- to k0 + n-
!       331. k+ + Omega- to pion+ + cascade-
!       332. k0 + Omega- to pion0 + cascade-
!       333. k- + Omega-ba to pion- + cascade-ba
!       334. k0- + Omega-ba to pion0 + cascade-ba
!       335. k0 + Omega- to pion- + cascade0
!       336. k+ + Omega- to pion0 + cascade0
!       337. k0- + Omega-ba to pion+ + cascade0-
!       338. k- + Omega-ba to pion0 + cascade0-
!       339. pion- + Omega- to k- + cascade-
!       340. pion0 + Omega- to k0- + cascade-
!       341. pion0 + Omega- to k- + cascade0
!       342. pion+ + Omega- to k0- + cascade0
!       343. pion+ + Omega-ba to k+ + cascade-ba
!       344. pion0 + Omega-ba to k0 + cascade-ba
!       345. pion0 + Omega-ba to k+ + cascade0-
!       346. pion- + Omega-ba to k0 + cascade0-
!       347. pion+ + delta- to pion- + p
!       348. pion+ + delta- to pion0 + n
!       349. pion+ + delta0 to pion+ + n
!       350. pion+ + delta0 to pion0 + p
!       351. pion+ + delta+ to pion+ + p
!       352. pion0 + delta++ to pion+ p
!       353. pion0 + delta+ to pion0 + p
!       354. pion0 + delta+ to pion+ + n
!       355. pion0 + delta0 to pion0 + n
!       356. pion0 + delta0 to pion- + p
!       357. pion0 + delta- to pion- + n
!       358. pion- + delta++ to pion0 + p
!       359. pion- + delta++ to pion+ + n
!       360. pion- + delta+ to pion- + p
!       361. pion- + delta+ to pion0 + n
!       362. pion- + delta0 to pion- + n
!       363. rho0 + n to pion- + p
!       364. rho0 + n to pion0 + n
!       365. rho- + n to pion- + n
!       366. rho+ + n to pion0 + p
!       367. rho+ + n to pion+ + n
!       368. rho0 + p to pion0 + p
!       369. rho0 + p to pion+ + n
!       370. rho- + p to pion0 + n
!       371. rho- + p to pion- + p
!       372. rho+ + p to pion+ + p
!       373. delta++ + n to p + p
!       374. delta+ + n to p + n
!       375. delta+ + p to p + p
!       376. delta0 + p to p + n
!       377. delta0 + n to n + n
!       378. delta- + p to n + n

!       follows are direct reactions again
!       379. psi' + pion0 to D + D*ba
!       380. psi' + pion- to D0 + D*ba
!       381. psi' + rho+ to D + D0ba
!       382. psi' + rho0 to D0 + D0ba
!       383. psi' + rho0 to D + Dba
!       384. psi' + rho- to D0 + Dba

!       D meson induced reactions.
!============================================
!       80 channels of D + pion/rho channels: 385-464
!       32 channels of D + N channels: 465-496
!       Total 112 channels: 385-496
!--------------------------------------------
!       D + pion/rho
!       Single production channels.
!
!       Charge ++, --
!       D+- + pion+-/rho+-
!       385. D+ + pi+  -->  D*+ + rho+
!       386. D- + pi-  -->  D*- + rho-
!       387. D*+ + pi+  -->  D+ + rho+
!       388. D*- + pi-  -->  D- + rho-
!       389. D+ + rho+  -->  D*+ + pi+
!       390. D- + rho-  -->  D*- + pi-
!       391. D*+ + rho+  -->  D+ + pi+
!       392. D*- + rho-  -->  D- + pi-
!
!       Charge 0- KF+-, Charge 0+ KF-+
!       D0 + pion+-/rho+-
!       393. D0 + pi-  -->  D*0 + rho-
!       394. Dbar0 + pi+  -->  Dbar*0 + rho+
!       395. D*0 + pi-  -->  D0 + rho-
!       396. D*bar0 + pi+  -->  Dbar0 + rho+
!       397. D0 + rho-  -->  D*0 + pi-
!       398. Dbar0 + rho+  -->  Dbar*0 + pi+
!       399. D*0 + rho-  -->  D0 + pi-
!       400. D*bar0 + rho+  -->  Dbar0 + pi+
!--------------------------------------------
!       D + pion/rho
!       Dual production channels.
!
!       Charge 00
!       D0 + pi0/rho0
!       401. D0 + pi0  -->  D*0 + rho0
!       402. D0 + pi0  -->  D*+ + rho-
!       403. Dbar0 + pi0  -->  Dbar*0 + rho0
!       404. Dbar0 + pi0  -->  Dbar*- + rho+
!       405. D*0 + pi0  -->  D0 + rho0
!       406. D*0 + pi0  -->  D+ + rho-
!       407. D*bar0 + pi0  -->  Dbar0 + rho0
!       408. D*bar0 + pi0  -->  D- + rho+
!       409. D0 + rho0  -->  D*0 + pi0
!       410. D0 + rho0  -->  D*+ + pi-
!       411. Dbar0 + rho0  -->  Dbar*0 + pi0
!       412. Dbar0 + rho0  -->  D*- + pi+
!       413. D*0 + rho0  -->  D0 + pi0
!       414. D*0 + rho0  -->  D+ + pi-
!       415. D*bar0 + rho0  -->  Dbar0 + pi0
!       416. D*bar0 + rho0  -->  D- + pi+
!
!       Charge +-, -+
!       D+- + pi+-/rho+-
!       417. D+ + pi-  -->  D*+ + rho-
!       418. D+ + pi-  -->  D*0 + rho0
!       419. D- + pi+  -->  D*- + rho+
!       420. D- + pi+  -->  D*bar0 + rho0
!       421. D*+ + pi-  -->  D+ + rho-
!       422. D*+ + pi-  -->  D0 + rho0
!       423. D*- + pi+  -->  D- + rho+
!       424. D*- + pi+  -->  Dbar0 + rho0
!       425. D+ + rho-  -->  D*+ + pi-
!       426. D+ + rho-  -->  D*0 + pi0
!       427. D- + rho+  -->  D*- + pi+
!       428. D- + rho+  -->  D*bar0 + pi0
!       429. D*+ + rho-  -->  D+ + pi-
!       430. D*+ + rho-  -->  D0 + pi0
!       431. D*- + rho+  -->  D- + pi+
!       432. D*- + rho+  -->  Dbar0 + pi0
!
!       Charge +0, -0
!       D+- + pi0/rho0
!       433. D+ + pi0  -->  D*+ + rho0
!       434. D+ + pi0  -->  D*0 + rho+
!       435. D- + pi0  -->  D*- + rho0
!       436. D- + pi0  -->  D*bar0 + rho-
!       437. D*+ + pi0  -->  D+ + rho0
!       438. D*+ + pi0  -->  D0 + rho+
!       439. D*- + pi0  -->  D- + rho0
!       440. D*- + pi0  -->  Dbar0 + rho-
!       441. D+ + rho0  -->  D*+ + pi0
!       442. D+ + rho0  -->  D*0 + pi+
!       443. D- + rho0  -->  D*- + pi0
!       444. D- + rho0  -->  D*bar0 + pi-
!       445. D*+ + rho0  -->  D+ + pi0
!       446. D*+ + rho0  -->  D0 + pi+
!       447. D*- + rho0  -->  D- + pi0
!       448. D*- + rho0  -->  Dbar0 + pi-
!
!       Charge 0+ KF++, Charge 0- KF --
!       D0 + pi+-/rho+-
!       449. D0 + pi+  -->  D*0 + rho+
!       450. D0 + pi+  -->  D*+ + rho0
!       451. Dbar0 + pi-  -->  Dbar*0 + rho-
!       452. Dbar0 + pi-  -->  D*- + rho0
!       453. D*0 + pi+  -->  D0 + rho+
!       454. D*0 + pi+  -->  D+ + rho0
!       455. D*bar0 + pi-  -->  Dbar0 + rho-
!       456. D*bar0 + pi-  -->  D- + rho0
!       457. D0 + rho+  -->  D*0 + pi+
!       458. D0 + rho+  -->  D*+ + pi0
!       459. Dbar0 + rho-  -->  Dbar*0 + pi-
!       460. Dbar0 + rho-  -->  D*- + pi0
!       461. D*0 + rho+  -->  D0 + pi+
!       462. D*0 + rho+  -->  D+ + pi0
!       463. D*bar0 + rho-  -->  Dbar0 + pi-
!       464. D*bar0 + rho-  -->  D- + pi0
!--------------------------------------------
!       D+- + N
!       465. D+  + p+    -->  D*+ + p+
!       466. D+  + p-    -->  D*+ + p-
!       467. D+  + n     -->  D*+ + n
!       468. D+  + nbar  -->  D*+ + nbar
!       469. D-  + p+    -->  D*- + p+
!       470. D-  + p-    -->  D*- + p-
!       471. D-  + n     -->  D*- + n
!       472. D-  + nbar  -->  D*- + nbar
!       473. D*+ + p+    -->  D+  + p+
!       474. D*+ + p-    -->  D+  + p-
!       475. D*+ + n     -->  D+  + n
!       476. D*+ + nbar  -->  D+  + nbar
!       477. D*- + p+    -->  D-  + p+
!       478. D*- + p-    -->  D-  + p-
!       479. D*- + n     -->  D-  + n
!       480. D*- + nbar  -->  D-  + nbar
!
!       D0 + N
!       481. D0     + p+    -->  D*0 + p+
!       482. D0     + p-    -->  D*0 + p-
!       483. D0     + n     -->  D*0 + n
!       484. D0     + nbar  -->  D*0 + nbar
!       485. Dbar0  + p+    -->  D*bar0 + p+
!       486. Dbar0  + p-    -->  D*bar0 + p-
!       487. Dbar0  + n     -->  D*bar0 + n
!       488. Dbar0  + nbar  -->  D*bar0 + nbar
!       489. D*0    + p+    -->  D0 + p+
!       490. D*0    + p-    -->  D0 + p-
!       491. D*0    + n     -->  D0 + n
!       492. D*0    + nbar  -->  D0 + nbar
!       493. D*bar0 + p+    -->  Dbar0 + p+
!       494. D*bar0 + p-    -->  Dbar0 + p-
!       495. D*bar0 + n     -->  Dbar0 + n
!       496. D*bar0 + nbar  -->  Dbar0 + nbar
!============================================
!       D meson induced reactions.

!       497-592: null.

!       593. lambda- + p to K*+ + Omega
!       594. lambda- + n to K*0 + Omega
!       595. sigma0- + p to K*+ + Omega
!       596. sigma0- + n to K*0 + Omega
!       597. p- + p to rho0 + Omega
!       598. p- + n to rho- + Omega
!       599. n- + p to rho+ + Omega
!       600. n- + n to rho0 + Omega

!       the reaction processes above are copy from A. Baldini, etal; ``Total
!        cross sections for reaction of high energy particles"; Spring-Veslay
!        Berlin, 1988.

        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        PARAMETER (KSZJ=80000,NSIZE=750000)
        PARAMETER (pio=3.141592653589793D0)
        common/sa1_h/nsa,non1,ksa(kszj,5),psa(kszj,5),vsa(kszj,5)
!       note the name of the arrays in 'sa1' in this subroutine
        common/sa10_h/csnn,cspin,cskn,cspipi,cspsn,cspsm,rcsit,ifram, &
         iabsb,iabsm,i_sigma_AQM,csspn,csspm
        common/ctllist_t/ lc(nsize,5), tc(nsize), tw(nsize)
        dimension pi(4), pj(4), pii(4), pjj(4), b(3)
        integer winel


        winel = 1
        kl    = ksa(l,2)
        kl1   = ksa(l1,2)
        am01  = PYMASS(kl)
        am02  = PYMASS(kl1)
!       The following statements guarantees the line number in particle
!        list of mesons is the first component of prod().
!       tw(i): the cross section ratio of (i-th inela.)/tot
        if( ABS(kl) < 1000 .OR. kl == 100443 .OR. kl == 100553 &
            .OR. kl == 10441 .OR. kl == 20443 )then
            call prod( l, l1, kl, kl1, ss, icp )
        else
            call prod( l1, l, kl1, kl, ss, icp )
        end if
!       the meson in a colliding pair should be first component of prod(),
!        but not for the case of collision occurs between two baryons (or
!        two mesons)
        ik1 = lc(icp,3)
        ik2 = lc(icp,4)
        w1  = tw(icp) / rcsit
!       1/rcsit : the cross section ratio of tot/(inela.)
        if( PYR(1) > w1 )then
            winel = 0
!       treated as elastic
            return
        end if

        icp5 = lc(icp,5)
        ! nonsense indeed
        if(icp5 >= 1593) return

        pj    = pj
        fi1   = PYANGL( pi(1), pi(2) )
        cta1  = PYANGL( pi(3), SQRT( pi(1)**2 + pi(2)**2 ) )
        cfi1  = COS(fi1)
        sfi1  = SIN(fi1)
        ccta1 = COS(cta1)
        scta1 = SIN(cta1)
        am1   = PYMASS(ik1)
        am2   = PYMASS(ik2)
        pp    = (ss**2 - (am1 + am2)**2)*(ss**2 - (am1 - am2)**2) / (4D0*ss**2)
        if(pp < 0D0) pp = 1D-10
        pp = SQRT(pp)
        pii(4) = ( ss**2 + am1**2 - am2**2 ) / ( 2D0*ss )
!       energy of one particle (between two) after scattering
        fis  = 2D0 * pio * PYR(1)
        cfis = COS(fis)
        sfis = SIN(fis)
!       quasi-elastic angular distribution for an inelastic scattering
        call cosin( am01, am1, ss, pi, pp, cctas )
!       isotropical distribution for an inelastic scattering
!       cctas = 2D0*PYR(1) - 1D0
        sctas = SQRT( 1D0 - cctas*cctas )
        pii(1) = cfi1 * ( ccta1*sctas*cfis + scta1*cctas ) - sfi1*sctas*sfis
        pii(2) = sfi1 * ( ccta1*sctas*cfis + scta1*cctas ) + cfi1*sctas*sfis
        pii(3) = ccta1*cctas - scta1*sctas*cfis
        pii(1) = pp*pii(1)
        pii(2) = pp*pii(2)
        pii(3) = pp*pii(3)
        do i=1,3,1
            pjj(i) = 0D0 - pii(i)
        end do
        pjj(4) = ss - pii(4)
        if( pii(4) < 0D0 .or. pjj(4) < 0D0 )then
            write(9,*) "error may happen here in subroutine" &
                    // "coinel(), energy is negative"
            if( pii(4) < 0D0 ) pii(4) = 1D-18
            if( pjj(4) < 0D0 ) pjj(4) = 1D-18
        end if
        ilo = 1
        call lorntz(ilo,b,pii,pjj)
        if( pii(4) < 0D0 .or. pjj(4) < 0D0 )then
            write(9,*) "error may happen here in subroutine" &
                    // "coinel(), energy is negative"
            if( pii(4) < 0D0 ) pii(4) = 1D-18
            if( pjj(4) < 0D0 ) pjj(4) = 1D-18
        end if


        return
        end



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine cosin(am01,am1,eij,pi,pp,cctas)
!!      Calculates cos(theta_s).
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        dimension pi(4)


        eij_dummy = eij
!       d = 3.65D0 * ( eij - am1 - am2 )
!       if( d < 1D-10 )return
        pt = 0.5D0
        a  = MIN( 10.3D0, 1D0 / ( 1.12D0 * pt ) / ( 1.12D0 * pt ) )
!       d6 = d**6
!       B  = d6 * A / ( 1D0 + d6 )
!       changed from 10.3 to 'A'
        B = A
!       if( B < 1D-20 )then
!           B = 1D-20
!       end if
        pm2  = pi(1)**2 + pi(2)**2 + pi(3)**2
        pm   = SQRT(pm2)
        em   = SQRT( pm2 + am01*am01 )
        em1  = SQRT( pp*pp + am1*am1 )
        tmin = am01**2 + am1**2 - 2D0 * ( em*em1 + pm*pp )
        tmax = am01**2 + am1**2 - 2D0 * ( em*em1 - pm*pp )
        abt = EXP( MAX( -7.0D2, B*tmin ) )
        abu = EXP( MAX( -7.0D2, B*tmax ) )
        cc  = PYR(1)
        tt1 = LOG( cc*abu + (1D0 - cc) * abt )
        tt = tt1 / B
        cctas = ( 0.5D0 * ( tt - am01**2 - am1**2 ) + em*em1 ) / ( pm * pp )
        if( ABS(cctas) > 1D0 ) cctas = SIGN( 1D0, cctas )


        return
        end



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine prod( l,l1, kl,kl1, ss, icp )
!!      Calculates particle production weight and fill up lc(i,3-5),tw(i).
!       tw : the ratio of cross section of (given inela.)/tot
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        INTEGER PYCHGE,PYCOMP
        PARAMETER(NSIZE=750000)
        COMMON/PYDAT2/KCHG(500,4),PMAS(500,4),PARF(2000),VCKM(4,4)
        COMMON/PYCIDAT2/KFMAXT,NONT2,PARAM(20),WEIGH(600)
        common/sa10_h/csnn,cspin,cskn,cspipi,cspsn,cspsm,rcsit,ifram, &
         iabsb,iabsm,i_sigma_AQM,csspn,csspm
        common/sa20_h/t0,cspipiKK,dep,ddt,edipi,epin,ecsnn,ekn,ecspsn,ecspsm, &
         rnt,rnp,ecsspn,ecsspm
        common/iloval/ilo1,ilo2,ilo3,ilo4,ilo5,ilo6,ilo7,ilo8,ilo9
        common/count_h/isinel(600)
        common/ctllist_t/ lc(nsize,5), tc(nsize), tw(nsize)
!       integer fact1, fact2


        fact1 = 0D0
        fact2 = 0D0
        para13 = PARAM(1) * 0.1D0 * PARAM(6) * 0.8D0
!       N-Nbar annihilation cross section is assumed to be equal to
!        0.8*(total inelastic cross section), while those of N-Ybar,
!        Nbar-Y annihilation are further divided by 5 (0.8*0.2*sigma_inel).
        ioo=0
        ilo=1
        ilo1=1
        ilo2=1
        ilo3=1
        ilo4=1
        ilo5=1
        ilo6=1
        ilo7=1
        ilo8=1
        ilo9=1

        nchargei = PYCHGE(kl) + PYCHGE(kl1)
! p-p ------------------------>
        if( kl ==  2212 .and. kl1 == 2212 )then
            call ppdelta(ss,icp,ioo)
            if(ioo == 0) goto 13
            goto 10
! p+n ----------------------->
        else if( ( kl ==  2212 .and. kl1 == 2112) .or. &
                 ( kl ==  2112 .and. kl1 == 2212) )then
            call pndelta(ss,icp,ioo)
            if(ioo == 0) goto 13
            goto 10
! n+n ------->
        else if( kl ==  2112 .and. kl1 == 2112 )then
            call nndelta(ss,icp,ioo)
            if(ioo == 0) goto 13
            goto 10
        end if

!       pion+ + pion-
        if( (kl == 211 .and. kl1 == -211) &
         .or.(kl == -211 .and. kl1 == 211) )then
            if(isinel(1) == 0)then
                fact1 = 0D0
                goto 101
            end if
            fact1 = 1D0
            ik3=-321
            ik4=321
            ic3=1
            call spipi(ik3,ik4,ss,ilo1)
            if(ilo1 == 0) fact1 = 0D0
101         if(isinel(2) == 0)then
                fact2 = 0D0
                goto 102
            endif
            fact2 = 1D0
            ik5=-311
            ik6=311
            ic5=2
            call spipi(ik5,ik6,ss,ilo2)
            if(ilo2 == 0) fact2 = 0D0
102         fact=fact1+fact2
            if( fact1 <= 1D-15 .AND. fact2 <= 1D-15 ) goto 13
!           if(ilo1 == 0 .and. ilo2 == 0) goto 13
            lc(icp,3)=ik3
            lc(icp,4)=ik4
            lc(icp,5)=ic3
            if(pyr(1) > fact1/fact)then
                lc(icp,3)=ik5
                lc(icp,4)=ik6
                lc(icp,5)=ic5
            endif
            tw(icp)=fact*cspipiKK/cspipi
            goto 10
        endif

!       pion+ + pion0
        if((kl == 211 .and. kl1 == 111) &
         .or.(kl == 111 .and. kl1 == 211))then
            if(isinel(3) == 0) goto 13
            ik3=-311
            ik4=321
            call spipi(ik3,ik4,ss,ilo1)
            if(ilo1 == 0) goto 13
            lc(icp,3)=ik3
            lc(icp,4)=ik4
            lc(icp,5)=3
            tw(icp)=cspipiKK/cspipi
            goto 10
        endif

!       pion- + pion0
        if((kl == -211 .and. kl1 == 111) &
         .or.(kl == 111 .and. kl1 == -211))then
        if(isinel(4) == 0) goto 13
        ik3=311
        ik4=-321
        call spipi(ik3,ik4,ss,ilo1)
        if(ilo1 == 0) goto 13
        lc(icp,3)=ik3
        lc(icp,4)=ik4
        lc(icp,5)=4
        tw(icp)=cspipiKK/cspipi
        goto 10
        endif

!       pion0 + pion0
        if(kl == 111 .and. kl1 == 111)then
        if(isinel(5) == 0)then
        fact1 = 0D0
        goto 103
        endif
        fact1 = 1D0
        ik3=-321
        ik4=321
        ic3=5
        call spipi(ik3,ik4,ss,ilo1)
        if(ilo1 == 0) fact1 = 0D0
103     if(isinel(6) == 0)then
        fact2 = 0D0
        goto 104
        endif
        fact2 = 1D0
        ik5=-311
        ik6=311
        ic5=6
        call spipi(ik5,ik6,ss,ilo2)
        if(ilo2 == 0) fact2 = 0D0
104     fact=fact1+fact2
        if( fact1 <= 1D-5 .AND. fact2 <= 1D-5 ) goto 13
!       if(ilo1 == 0 .and. ilo2 == 0) goto 13
        lc(icp,3)=ik3
        lc(icp,4)=ik4
        lc(icp,5)=ic3
        if(pyr(1) > fact1/fact)then
        lc(icp,3)=ik5
        lc(icp,4)=ik6
        lc(icp,5)=ic5
        endif
        tw(icp)=fact*cspipiKK/cspipi
        goto 10
        endif

!       pion+ + p
        if(kl == 211 .and. kl1 == 2212)then
        call pip2(ss,icp,ioo)
        if(ioo == 0) goto 13
        goto 10
        endif

!       pion+ + n
        if(kl == 211 .and. kl1 == 2112)then
        call pin1(ss,icp,ioo)
        if(ioo == 0) goto 13
        goto 10
        endif

!       pion- + p
        if(kl == -211 .and. kl1 == 2212)then
        call pip1(ss,icp,ioo)
        if(ioo == 0) goto 13
        goto 10
        endif

!       pion- + n
        if(kl == -211 .and. kl1 == 2112)then
        call pin3(ss,icp,ioo)
        if(ioo == 0) goto 13
        goto 10
        endif

!       pion0 + p
        if(kl == 111 .and. kl1 == 2212)then
        call pip3(ss,icp,ioo)
        if(ioo == 0) goto 13
        goto 10
        endif

!       pion0 + n
        if(kl == 111 .and. kl1 == 2112)then
        call pin2(ss,icp,ioo)
        if(ioo == 0) goto 13
        goto 10
        endif

!       pion+ + pba
        if(kl == 211 .and. kl1 == -2212)then
        if(isinel(21) == 0)then
        si1=0.
        goto 319
        endif
        the=pmas(pycomp(-311),1)+pmas(pycomp(-3122),1)
        si1=s1724(ss,ilo1,0,the)
!       cross section of pion+ + pba to k0- + lambdaba
!       cross section of pion+ + pba to k0- + lambdaba is assumed to be
!        equal to pion- + p to k0 + lambda = s1724
319     if(isinel(22) == 0)then
        si2=0D0
        goto 320
        endif
        the=pmas(pycomp(-311),1)+pmas(pycomp(-3212),1)
        si2=s1724(ss,ilo2,0,the)
!       cross section of pion+ + pba to k0- + sigma0-
320     if(isinel(23) == 0)then
        si3=0D0
        goto 321
        endif
        the=pmas(pycomp(-321),1)+pmas(pycomp(-3112),1)
        si3=s1724(ss,ilo3,0,the)
!       cross section of pion+ + pba to k- + sigma-ba
321     if(ilo1 == 0.and.ilo2 == 0.and.ilo3 == 0) goto 13
        if(si1 < 1D-6.and.si2 < 1D-6.and.si3 < 1D6) goto 13
        si12=si1+si2
        sit=si12+si3
        s1=si1/sit
        s2=si12/sit
        rlus=pyr(1)
        if(rlus <= s1)then
        ik1=-311
        ik2=-3122
        ic=21
!       pion+ + pba to k0- + lambda-
        goto 322
        endif
        if(rlus > s1 .and. rlus <= s2)then
        ik1=-311
        ik2=-3212
        ic=22
!       pion+ + pba to k0- + sigma0-
        goto 322
        endif
        ik1=-321
        ik2=-3112
        ic=23
!       pion+ + pba to k- + sigma-ba
322     lc(icp,3)=ik1
        lc(icp,4)=ik2
        lc(icp,5)=ic
        tw(icp)=sit/cspin/10D0
        goto 10
        endif

!       pion+ + nba to k0- + sigma-ba
        if(kl == 211 .and. kl1 == -2112)then
        if(isinel(24) == 0) goto 13
        the=pmas(pycomp(-311),1)+pmas(pycomp(-3112),1)
        tw(icp)=s1724(ss,ilo,0,the)/cspin/10D0
!       cross section of pion+ + nba to k0- + sigma-ba is assumed to be
!        equal to pion- + n to k + y (isotropic averaged)
        if(ilo == 0) goto 13
        lc(icp,3)=-311
        lc(icp,4)=-3112
        lc(icp,5)=24
        goto 10
        endif

!       pion- + pba to k- + sigma+ba
        if(kl == -211 .and. kl1 == -2212)then
        if(isinel(25) == 0) goto 13
        the=pmas(pycomp(-321),1)+pmas(pycomp(-3222),1)
        tw(icp)=s1724(ss,ilo,0,the)/cspin/10D0
!       cross section of pion- + pba to k- + sigma+ba is assumed to be
!        equal to pion+ + p to k + y (isotropic averaged)
        if(ilo == 0) goto 13
        lc(icp,3)=-321
        lc(icp,4)=-3222
        lc(icp,5)=25
        goto 10
        endif

!       pion- + nba
        if(kl == -211 .and. kl1 == -2112)then
        if(isinel(26) == 0)then
        si1=0D0
        goto 323
        endif
        the=pmas(pycomp(-321),1)+pmas(pycomp(-3122),1)
        si1=s1724(ss,ilo1,0,the)
!       cross section of pion- + nba to k- + y- is assumed to be
!        equal to pion+ + n to k + y
323     if(isinel(27) == 0)then
        si2=0D0
        goto 324
        endif
        the=pmas(pycomp(-321),1)+pmas(pycomp(-3212),1)
        si2=s1724(ss,ilo2,0,the)
324     if(isinel(28) == 0)then
        si3=0D0
        goto 325
        endif
        the=pmas(pycomp(-311),1)+pmas(pycomp(-3222),1)
        si3=s1724(ss,ilo3,0,the)
325     if(ilo1 == 0 .and. ilo2 == 0 .and. ilo3 == 0) goto 13
        if(si1 < 1D-6.and.si2 < 1D-6.and.si3 < 1D-6) goto 13
        si12=si1+si2
        sit=si12+si3
        s1=si1/sit
        s2=si12/sit
        rlus=pyr(1)
        if(rlus <= s1)then
        ik1=-321
        ik2=-3122
        ic=26
!       pion- + nba to k- + lambdaba
        goto 326
        endif
        if(rlus > s1 .and. rlus <= s2)then
        ik1=-321
        ik2=-3212
        ic=27
!       pion- + nba to k- + sigma0ba
        goto 326
        endif
        ik1=-311
        ik2=-3222
        ic=28
!       pion- + nba to k0- + sigma+ba
326     lc(icp,3)=ik1
        lc(icp,4)=ik2
        lc(icp,5)=ic
        tw(icp)=sit/cspin/10D0
        goto 10
        endif

!       pion0 + pba
        if(kl == 111 .and. kl1 == -2212)then
        if(isinel(29) == 0)then
        si1=0D0
        goto 327
        endif
        the=pmas(pycomp(-321),1)+pmas(pycomp(-3122),1)
        si1=s1724(ss,ilo1,0,the)
327     if(isinel(30) == 0)then
        si2=0D0
        goto 328
        endif
        the=pmas(pycomp(-321),1)+pmas(pycomp(-3212),1)
        si2=s1724(ss,ilo2,0,the)
328     if(isinel(31) == 0)then
        si3=0D0
        goto 329
        endif
        the=pmas(pycomp(-311),1)+pmas(pycomp(-3222),1)
        si3=s1724(ss,ilo3,0,the)
329     if(ilo1 == 0.and.ilo2 == 0.and.ilo3 == 0) goto 13
        if(si1 < 1D-6.and.si2 < 1D-6.and.si3 < 1D-6) goto 13
        si12=si1+si2
        sit=si12+si3
        s1=si1/sit
        s2=si12/sit
        rlus=pyr(1)
        if(rlus <= s1)then
        ik1=-321
        ik2=-3122
        ic=29
!       pion0 + pba to k- + lambda-
        goto 330
        endif
        if(rlus > s1 .and. rlus <= s2)then
        ik1=-321
        ik2=-3212
        ic=30
!       pion0 + pba to k- + sigma0-
        goto 330
        endif
        ik1=-311
        ik2=-3222
        ic=31
!       pion0 + pba to k0- + sigma+ba
330     lc(icp,3)=ik1
        lc(icp,4)=ik2
        lc(icp,5)=ic
        tw(icp)=sit/cspin/10D0
        goto 10
        endif

!       pion0 + nba
        if(kl == 111 .and. kl1 == -2112)then
        if(isinel(32) == 0)then
        si1=0D0
        goto 331
        endif
        the=pmas(pycomp(-321),1)+pmas(pycomp(-3112),1)
        si1=s1724(ss,ilo1,0,the)
331     if(isinel(33) == 0)then
        si2=0D0
        goto 332
        endif
        the=pmas(pycomp(-311),1)+pmas(pycomp(-3122),1)
        si2=s1724(ss,ilo2,0,the)
332     if(isinel(34) == 0)then
        si3=0D0
        goto 333
        endif
        the=pmas(pycomp(-311),1)+pmas(pycomp(-3212),1)
        si3=s1724(ss,ilo3,0,the)
333     if(ilo1 == 0.and.ilo2 == 0.and.ilo3 == 0) goto 13
        if(si1 < 1D-6.and.si2 < 1D-6.and.si3 < 1D-6) goto 13
        si12=si1+si2
        sit=si12+si3
        s1=si1/sit
        s2=si12/sit
        rlus=pyr(1)
        if(rlus <= s1)then
        ik1=-321
        ik2=-3112
        ic=32
!       pion0 + nba to k- + sigma-ba
        goto 334
        endif
        if(rlus > s1 .and. rlus <= s2)then
        ik1=-311
        ik2=-3122
        ic=33
!       pion0 + nba to k0- + lambdaba
        goto 334
        endif
        ik1=-311
        ik2=-3212
        ic=34
!       pion0 + nba to k0- + sigma0ba
334     lc(icp,3)=ik1
        lc(icp,4)=ik2
        lc(icp,5)=ic
        tw(icp)=sit/cspin/10D0
        goto 10
        endif

!       pion+ + sigma-
        if(kl == 211 .and. kl1 == 3112)then
        if(isinel(35) == 0)then
        si1=0D0
        goto 701
        endif
        the=pmas(pycomp(321),1)+pmas(pycomp(3312),1)
        si1=s1724(ss,ilo1,0,the)
!       cross section of pion + y to Kaon + cascade is assumed to be
!        equal to pion + n to Kaon + y,but take the different of
!        threshold energy into account

701     if(isinel(97) == 0)then
        si2=0D0
        goto 660
        endif
        the=pmas(pycomp(311),1)+pmas(pycomp(3322),1)
        si2=s1724(ss,ilo2,0,the)

660     if(isinel(249) == 0)then
        si3=0D0
        goto 661
        endif
        the=pmas(pycomp(-321),1)+pmas(pycomp(2212),1)
        ik1=-321
        ik2=2212
        call srev(kl,kl1,ik1,ik2,ss,ilo3,fac, &
         1.,1.,0.,0.5,0.5,0.5,0.,0.5,1.)
        si3=10D0*s1724(ss,ilo3,0,the)*fac
661     if(isinel(250) == 0)then
        si4=0D0
        goto 702
        endif
        the=pmas(pycomp(-311),1)+pmas(pycomp(2112),1)
        ik1=-311
        ik2=2112
        call srev(kl,kl1,ik1,ik2,ss,ilo4,fac, &
         1.,1.,0.,0.5,0.5,0.5,0.,0.5,1.)
        si4=10D0*s1724(ss,ilo4,0,the)*fac
702     if(ilo1 == 0.and.ilo2 == 0.and.ilo3 == 0.and.ilo4 == 0) goto 13
        if(si1 < 1D-6.and.si2 < 1D-6.and.si3 < 1D-6.and. &
         si4 < 1D-6) goto 13
        si12=si1+si2
        si13=si12+si3
        si14=si13+si4
        s1=si1/si14
        s2=si12/si14
        s3=si13/si14
        rlus=pyr(1)
        if(rlus <= s1)then
        ik1=321
        ik2=3312
        ic=35
!       pion+ + sigma- to k+ + cascade-
        goto 703
        endif
        if(rlus > s1 .and. rlus <= s2)then
        ik1=311
        ik2=3322
        ic=97
!       pion+ + sigma- to k0 + cascade0
        goto 703
        endif
        if(rlus > s2 .and. rlus <= s3)then
        ik1=-321
        ik2=2212
        ic=249
!       pion+ + sigma- to k- + p
        goto 703
        endif
        ik1=-311
        ik2=2112
        ic=250
!       pion+ + sigma- to k0- + n
703     lc(icp,3)=ik1
        lc(icp,4)=ik2
        lc(icp,5)=ic
        tw(icp)=si14/cspin/10D0
        goto 10
        endif

!       pion- + sigma-ba
        if(kl == -211 .and. kl1 == -3112)then
        if(isinel(45) == 0)then
        si1=0D0
        goto 7011
        endif
        the=pmas(pycomp(-321),1)+pmas(pycomp(-3312),1)
        si1=s1724(ss,ilo1,0,the)
!       cross section of pion + y to Kaon + cascade is assumed to be
!        equal to pion + n to Kaon + y,but take the different of
!        threshold energy into account

7011    if(isinel(105) == 0)then
        si2=0D0
        goto 6601
        endif
        the=pmas(pycomp(-311),1)+pmas(pycomp(-3322),1)
        si2=s1724(ss,ilo2,0,the)

6601    if(isinel(271) == 0)then
        si3=0D0
        goto 6611
        endif
        the=pmas(pycomp(321),1)+pmas(pycomp(-2212),1)
        ik1=321
        ik2=-2212
        call srev(kl,kl1,ik1,ik2,ss,ilo3,fac, &
         1.,1.,0.,0.5,0.5,0.5,0.,0.5,1.)
        si3=10D0*s1724(ss,ilo3,0,the)*fac
6611    if(isinel(272) == 0)then
        si4=0D0
        goto 7021
        endif
        the=pmas(pycomp(311),1)+pmas(pycomp(-2112),1)
        ik1=311
        ik2=-2112
        call srev(kl,kl1,ik1,ik2,ss,ilo4,fac, &
         1.,1.,0.,0.5,0.5,0.5,0.,0.5,1.)
        si4=10D0*s1724(ss,ilo4,0,the)*fac
7021    if(ilo1 == 0.and.ilo2 == 0.and.ilo3 == 0 &
         .and.ilo4 == 0) goto 13
        if(si1 < 1D-6.and.si2 < 1D-6.and.si3 < 1D-6 &
         .and.si4 < 1D-6) goto 13
        si12=si1+si2
        si13=si12+si3
        si14=si13+si4
        s1=si1/si14
        s2=si12/si14
        s3=si13/si14
        rlus=pyr(1)
        if(rlus <= s1)then
        ik1=-321
        ik2=-3312
        ic=45
!       pion- + sigma-ba to k- + cascade-ba
        goto 7031
        endif
        if(rlus > s1 .and. rlus <= s2)then
        ik1=-311
        ik2=-3322
        ic=105
!       pion- + sigma-ba to k0- + cascade0-
        goto 7031
        endif
        if(rlus > s2 .and. rlus <= s3)then
        ik1=321
        ik2=-2212
        ic=271
!       pion- + sigma-ba to k+ + pba
        goto 7031
        endif
        ik1=311
        ik2=-2112
        ic=272
!       pion- + sigma-ba to k0 + nba
7031    lc(icp,3)=ik1
        lc(icp,4)=ik2
        lc(icp,5)=ic
        tw(icp)=si14/cspin/10D0
        goto 10
        endif

!       pion- + sigma+
        if(kl == -211 .and. kl1 == 3222)then
        if(isinel(37) == 0)then
        si1=0D0
        goto 704
        endif
        the=pmas(pycomp(321),1)+pmas(pycomp(3312),1)
        si1=s1724(ss,ilo1,0,the)

704     if(isinel(100) == 0)then
        si2=0D0
        goto 663
        endif
        the=pmas(pycomp(311),1)+pmas(pycomp(3322),1)
        si2=s1724(ss,ilo2,0,the)

663     if(isinel(254) == 0)then
        si3=0D0
        goto 664
        endif
        ik1=-321
        ik2=2212
        the=pmas(pycomp(-321),1)+pmas(pycomp(2212),1)
        call srev(kl,kl1,ik1,ik2,ss,ilo3,fac, &
         1.,1.,0.,0.5,0.5,0.5,0.,0.5,1.)
        si3=10D0*s1724(ss,ilo3,0,the)*fac
664    if(isinel(255) == 0)then
        si4=0D0
        goto 705
        endif
        ik1=-311
        ik2=2112
        the=pmas(pycomp(-311),1)+pmas(pycomp(2112),1)
        call srev(kl,kl1,ik1,ik2,ss,ilo4,fac, &
         1.,1.,0.,0.5,0.5,0.5,0.,0.5,1.)
        si4=10D0*s1724(ss,ilo4,0,the)*fac
705     if(ilo1 == 0.and.ilo2 == 0.and.ilo3 == 0.and.ilo4 == 0) goto 13
        if(si1 < 1D-6.and.si2 < 1D-6.and.si3 < 1D-6 &
        .and.si4 < 1D-6) goto 13
        si12=si1+si2
        si13=si12+si3
        si14=si13+si4
        s1=si1/si14
        s2=si12/si14
        s3=si13/si14
        rlus=pyr(1)
        if(rlus <= s1)then
        ik1=321
        ik2=3312
        ic=37
!       pion- + sigma+ to k+ + cascade-
        goto 706
        endif
        if(rlus > s1 .and. rlus <= s2)then
        ik1=311
        ik2=3322
        ic=100
!       pion- + sigma+ to k0 + cascade0
        goto 706
        endif
        if(rlus > s2 .and. rlus <= s3)then
        ik1=-321
        ik2=2212
        ic=254
!       pion- + sigma+ to k- + p
        goto 706
        endif
        ik1=-311
        ik2=2112
        ic=255
!       pion- + sigma+ to k0- + n
706     lc(icp,3)=ik1
        lc(icp,4)=ik2
        lc(icp,5)=ic
        tw(icp)=si14/cspin/10D0
        goto 10
        endif

!       pion+ + sigma+bar
        if(kl == 211 .and. kl1 == -3222)then
        if(isinel(43) == 0)then
        si1=0D0
        goto 7041
        endif
        the=pmas(pycomp(-321),1)+pmas(pycomp(-3312),1)
        si1=s1724(ss,ilo1,0,the)

7041    if(isinel(104) == 0)then
        si2=0D0
        goto 6631
        endif
        the=pmas(pycomp(-311),1)+pmas(pycomp(-3322),1)
        si2=s1724(ss,ilo2,0,the)

6631    if(isinel(266) == 0)then
        si3=0D0
        goto 6641
        endif
        ik1=321
        ik2=-2212
        the=pmas(pycomp(321),1)+pmas(pycomp(-2212),1)
        call srev(kl,kl1,ik1,ik2,ss,ilo3,fac, &
         1.,1.,0.,0.5,0.5,0.5,0.,0.5,1.)
        si3=10D0*s1724(ss,ilo3,0,the)*fac
6641    if(isinel(267) == 0)then
        si4=0D0
        goto 7051
        endif
        ik1=311
        ik2=-2112
        the=pmas(pycomp(311),1)+pmas(pycomp(-2112),1)
        call srev(kl,kl1,ik1,ik2,ss,ilo4,fac, &
         1.,1.,0.,0.5,0.5,0.5,0.,0.5,1.)
        si4=10D0*s1724(ss,ilo4,0,the)*fac
7051    if(ilo1 == 0.and.ilo2 == 0.and.ilo3 == 0.and.ilo4 == 0) goto 13
        if(si1 < 1D-6.and.si2 < 1D-6.and.si3 < 1D-6 &
         .and.si4 < 1D-6) goto 13
        si12=si1+si2
        si13=si12+si3
        si14=si13+si4
        s1=si1/si14
        s2=si12/si14
        s3=si13/si14
        rlus=pyr(1)
        if(rlus <= s1)then
        ik1=-321
        ik2=-3312
        ic=43
!       pion+ + sigma+bar to k- + cascade-bar
        goto 7061
        endif
        if(rlus > s1 .and. rlus <= s2)then
        ik1=-311
        ik2=-3322
        ic=104
!       pion+ + sigma+bar to k0- + cascade0bar
        goto 7061
        endif
        if(rlus > s2 .and. rlus <= s3)then
        ik1=321
        ik2=-2212
        ic=266
!       pion+ + sigma+bar to k+ + pbar
        goto 7061
        endif
        ik1=311
        ik2=-2112
        ic=267
!       pion+ + sigma+bar to k0 + nbar
7061    lc(icp,3)=ik1
        lc(icp,4)=ik2
        lc(icp,5)=ic
        tw(icp)=si14/cspin/10D0
        goto 10
        endif
!       pion- + sigma0
        if(kl == -211 .and. kl1 == 3212)then
        if(isinel(38) == 0)then
        si1=0D0
        goto 1681
        endif
        the=pmas(pycomp(311),1)+pmas(pycomp(3312),1)
        si1=s1724(ss,ilo1,0,the)
1681    if(isinel(256) == 0)then
        si2=0D0
        goto 1682
        endif
        the=pmas(pycomp(-321),1)+pmas(pycomp(2112),1)
        ik3=-321
        ik4=2112
        call srev(kl,kl1,ik3,ik4,ss,ilo2,fac, &
         1.,1.,0.,0.5,0.5,0.5,0.,0.5,1.)
        si2=s1724(ss,ilo2,0,the)*10D0*fac
1682    if(ilo1 == 0.and.ilo2 == 0) goto 13
        if(si1 < 1D-6.and.si2 < 1D-6) goto 13
        sit=si1+si2
        s1=si1/sit
        rlus=pyr(1)
        if(rlus <= s1)then
        ik1=311
        ik2=3312
        ic=38
!       pion- + sigma0 to k0 + cascade-
        goto 683
        endif
        ik1=-321
        ik2=2112
        ic=256
!       pion- + sigma0 to k- + n
683     lc(icp,3)=ik1
        lc(icp,4)=ik2
        lc(icp,5)=ic
        tw(icp)=sit/cspin/10D0
        goto 10
        endif


!       pion+ + sigma0 bar
        if(kl == 211 .and. kl1 == -3212)then
        if(isinel(44) == 0)then
        si1=0D0
        goto 1781
        endif
        the=pmas(pycomp(-311),1)+pmas(pycomp(-3312),1)
        si1=s1724(ss,ilo1,0,the)
1781    if(isinel(268) == 0)then
        si2=0D0
        goto 1782
        endif
        the=pmas(pycomp(321),1)+pmas(pycomp(-2112),1)
        ik3=321
        ik4=-2112
        call srev(kl,kl1,ik3,ik4,ss,ilo2,fac, &
         1.,1.,0.,0.5,0.5,0.5,0.,0.5,1.)
        si2=s1724(ss,ilo2,0,the)*10D0*fac
1782    if(ilo1 == 0.and.ilo2 == 0) goto 13
        if(si1 < 1D-6.and.si2 < 1D-6) goto 13
        sit=si1+si2
        s1=si1/sit
        rlus=pyr(1)
        if(rlus <= s1)then
        ik1=-311
        ik2=-3312
        ic=44
!       pion+ + sigma0bar to k0- + cascade-bar
        goto 6831
        endif
        ik1=321
        ik2=-2112
        ic=268
!       pion+ + sigma0bar to k+ + nbar
6831    lc(icp,3)=ik1
        lc(icp,4)=ik2
        lc(icp,5)=ic
        tw(icp)=sit/cspin/10D0
        goto 10
        endif


!       pion0 + sigma0
        if(kl == 111 .and. kl1 == 3212)then
        if(isinel(41) == 0)then
        si1=0D0
        goto 710
        endif
        the=pmas(pycomp(321),1)+pmas(pycomp(3312),1)
        si1=s1724(ss,ilo1,0,the)

710     if(isinel(102) == 0)then
        si2=0D0
        goto 666
        endif
        the=pmas(pycomp(311),1)+pmas(pycomp(3322),1)
        si2=s1724(ss,ilo2,0,the)

666     if(isinel(261) == 0)then
        si3=0D0
        goto 667
        endif
        the=pmas(pycomp(-321),1)+pmas(pycomp(2212),1)
        ik1=-321
        ik2=2212
        call srev(kl,kl1,ik1,ik2,ss,ilo3,fac, &
         1.,1.,0.,0.5,0.5,0.5,0.,0.5,1.)
        si3=10D0*s1724(ss,ilo3,0,the)*fac
667     if(isinel(262) == 0)then
        si4=0D0
        goto 711
        endif
        the=pmas(pycomp(-311),1)+pmas(pycomp(2112),1)
        ik1=-311
        ik2=2112
        call srev(kl,kl1,ik1,ik2,ss,ilo4,fac, &
         1.,1.,0.,0.5,0.5,0.5,0.,0.5,1.)
        si4=10D0*s1724(ss,ilo4,0,the)*fac
711     if(ilo1 == 0.and.ilo2 == 0.and.ilo3 == 0.and.ilo4 == 0) goto 13
        if(si1 < 1D-6.and.si2 < 1D-6.and.si3 < 1D-6 &
         .and.si4 < 1D-6) goto 13
        si12=si1+si2
        si13=si12+si3
        si14=si13+si4
        s1=si1/si14
        s2=si12/si14
        s3=si13/si14
        rlus=pyr(1)
        if(rlus <= s1)then
        ik1=321
        ik2=3312
        ic=41
!       pion0 + sigma0 to k+ + cascade-
        goto 712
        endif
        if(rlus > s1 .and. rlus <= s2)then
        ik1=311
        ik2=3322
        ic=102
!       pion0 + sigma0 to k0 + cascade0
        goto 712
        endif
        if(rlus > s2 .and. rlus <= s3)then
        ik1=-321
        ik2=2212
        ic=261
!      pion0 + sigma0 to k- + p
        goto 712
        endif
       ik1=-311
        ik2=2112
        ic=262
!       pion0 + sigma0 to k0- + n
712     lc(icp,3)=ik1
        lc(icp,4)=ik2
        lc(icp,5)=ic
        tw(icp)=si14/cspin/10D0
        goto 10
        endif
!       pion0 + sigma0bar
        if(kl == 111 .and. kl1 == -3212)then
        if(isinel(48) == 0)then
        si1=0D0
        goto 1910
        endif
        the=pmas(pycomp(-321),1)+pmas(pycomp(-3312),1)
        si1=s1724(ss,ilo1,0,the)

1910    if(isinel(109) == 0)then
        si2=0D0
        goto 1666
        endif
        the=pmas(pycomp(-311),1)+pmas(pycomp(-3322),1)
        si2=s1724(ss,ilo2,0,the)

1666    if(isinel(278) == 0)then
        si3=0D0
        goto 1667
        endif
        the=pmas(pycomp(321),1)+pmas(pycomp(-2212),1)
        ik1=321
        ik2=-2212
        call srev(kl,kl1,ik1,ik2,ss,ilo3,fac, &
         1.,1.,0.,0.5,0.5,0.5,0.,0.5,1.)
        si3=10D0*s1724(ss,ilo3,0,the)*fac
1667    if(isinel(279) == 0)then
        si4=0D0
        goto 5711
        endif
        the=pmas(pycomp(311),1)+pmas(pycomp(-2112),1)
        ik1=311
        ik2=-2112
        call srev(kl,kl1,ik1,ik2,ss,ilo4,fac, &
         1.,1.,0.,0.5,0.5,0.5,0.,0.5,1.)
        si4=10D0*s1724(ss,ilo4,0,the)*fac
5711    if(ilo1 == 0.and.ilo2 == 0.and.ilo3 == 0.and.ilo4 == 0) goto 13
        if(si1 < 1D-6.and.si2 < 1D-6.and.si3 < 1D-6 &
         .and.si4 < 1D-6) goto 13
        si12=si1+si2
        si13=si12+si3
        si14=si13+si4
        s1=si1/si14
        s2=si12/si14
        s3=si13/si14
        rlus=pyr(1)
        if(rlus <= s1)then
        ik1=-321
        ik2=-3312
        ic=48
!       pion0 + sigma0bar to k- + cascade-bar
        goto 7125
        endif
        if(rlus > s1 .and. rlus <= s2)then
        ik1=-311
        ik2=-3322
        ic=109
!       pion0 + sigma0bar to k0- + cascade0bar
        goto 7125
        endif
        if(rlus > s2 .and. rlus <= s3)then
        ik1=321
        ik2=-2212
        ic=278
!      pion0 + sigma0bar to k+ + pbar
        goto 7125
        endif
        ik1=311
        ik2=-2112
        ic=279
!       pion0 + sigma0bar to k0 + nbar
7125    lc(icp,3)=ik1
        lc(icp,4)=ik2
        lc(icp,5)=ic
        tw(icp)=si14/cspin/10D0
        goto 10
        endif

!       pion+ + sigma0
        if(kl == 211 .and. kl1 == 3212)then
        if(isinel(98) == 0)then
        si1=0D0
        goto 2681
        endif
        the=pmas(pycomp(321),1)+pmas(pycomp(3322),1)
        si1=s1724(ss,ilo1,0,the)
2681    if(isinel(251) == 0)then
        si2=0D0
        goto 2682
        endif
        the=pmas(pycomp(-311),1)+pmas(pycomp(2212),1)
        ik3=-311
        ik4=2212
        call srev(kl,kl1,ik3,ik4,ss,ilo2,fac, &
         1.,1.,0.,0.5,0.5,0.5,0.,0.5,1.)
        si2=s1724(ss,ilo2,0,the)*10D0*fac
2682    if(ilo1 == 0.and.ilo2 == 0) goto 13
        if(si1 < 1D-6.and.si2 < 1D-6) goto 13
        sit=si1+si2
        s1=si1/sit
        rlus=pyr(1)
        if(rlus <= s1)then
        ik1=321
        ik2=3322
        ic=98
!       pion+ + sigma0 to k+ + cacade0
        goto 1683
        endif
        ik1=-311
        ik2=2212
        ic=251
!       pion+ + sigma0 to k0- + p
1683    lc(icp,3)=ik1
        lc(icp,4)=ik2
        lc(icp,5)=ic
        tw(icp)=sit/cspin/10D0
        goto 10
        endif
!       pion- + sigma0-
        if(kl == -211 .and. kl1 == -3212)then
        if(isinel(106) == 0)then
        si1=0D0
        goto 2781
        endif
        the=pmas(pycomp(-321),1)+pmas(pycomp(-3322),1)
        si1=s1724(ss,ilo1,0,the)
2781    if(isinel(273) == 0)then
        si2=0D0
        goto 2782
        endif
        the=pmas(pycomp(311),1)+pmas(pycomp(-2212),1)
        ik3=311
        ik4=-2212
        call srev(kl,kl1,ik3,ik4,ss,ilo2,fac, &
         1.,1.,0.,0.5,0.5,0.5,0.,0.5,1.)
        si2=s1724(ss,ilo2,0,the)*10D0*fac
2782    if(ilo1 == 0.and.ilo2 == 0) goto 13
        if(si1 < 1D-6.and.si2 < 1D-6) goto 13
        sit=si1+si2
        s1=si1/sit
        rlus=pyr(1)
        if(rlus <= s1)then
        ik1=-321
        ik2=-3322
        ic=106
!       pion- + sigma0- to k- + cacade0-
        goto 1783
        endif
        ik1=311
        ik2=-2212
        ic=273
!       pion- + sigma0- to k0 + pba
1783    lc(icp,3)=ik1
        lc(icp,4)=ik2
        lc(icp,5)=ic
        tw(icp)=sit/cspin/10D0
        goto 10
        endif

!       pion0 + sigma+
        if(kl == 111 .and. kl1 == 3222)then
        if(isinel(101) == 0)then
        si1=0D0
        goto 3681
        endif
        the=pmas(pycomp(321),1)+pmas(pycomp(3322),1)
        si1=s1724(ss,ilo1,0,the)
3681    if(isinel(259) == 0)then
        si2=0D0
        goto 3682
        endif
        the=pmas(pycomp(-311),1)+pmas(pycomp(2212),1)
        ik3=-311
        ik4=2212
        call srev(kl,kl1,ik3,ik4,ss,ilo2,fac, &
         1.,1.,0.,0.5,0.5,0.5,0.,0.5,1.)
        si2=s1724(ss,ilo2,0,the)*10D0*fac
3682    if(ilo1 == 0.and.ilo2 == 0) goto 13
        if(si1 < 1D-6.and.si2 < 1D-6) goto 13
        sit=si1+si2
        s1=si1/sit
        rlus=pyr(1)
        if(rlus <= s1)then
        ik1=321
        ik2=3322
        ic=101
!       pion0 + sigma+ to k+ + cascade0
        goto 7683
        endif
        ik1=-311
        ik2=2212
        ic=259
!       pion0 + sigma+ to k0- + p
7683    lc(icp,3)=ik1
        lc(icp,4)=ik2
        lc(icp,5)=ic
        tw(icp)=sit/cspin/10D0
        goto 10
        endif
!       pion0 + sigma+bar
        if(kl == 111 .and. kl1 == -3222)then
        if(isinel(108) == 0)then
        si1=0D0
        goto 3781
        endif
        the=pmas(pycomp(-321),1)+pmas(pycomp(-3322),1)
        si1=s1724(ss,ilo1,0,the)
3781    if(isinel(276) == 0)then
        si2=0D0
        goto 3782
        endif
        the=pmas(pycomp(311),1)+pmas(pycomp(-2212),1)
        ik3=311
        ik4=-2212
        call srev(kl,kl1,ik3,ik4,ss,ilo2,fac, &
         1.,1.,0.,0.5,0.5,0.5,0.,0.5,1.)
        si2=s1724(ss,ilo2,0,the)*10D0*fac
3782    if(ilo1 == 0.and.ilo2 == 0) goto 13
        if(si1 < 1D-6.and.si2 < 1D-6) goto 13
        sit=si1+si2
        s1=si1/sit
        rlus=pyr(1)
        if(rlus <= s1)then
        ik1=-321
        ik2=-3322
        ic=108
!       pion0 + sigma+bar to k- + cascade0bar
        goto 6783
        endif
        ik1=311
        ik2=-2212
        ic=276
!       pion0 + sigma+bar to k0 + pbar
6783    lc(icp,3)=ik1
        lc(icp,4)=ik2
        lc(icp,5)=ic
        tw(icp)=sit/cspin/10D0
        goto 10
        endif

!       pion0 + sigma-
        if(kl == 111 .and. kl1 == 3112)then
        if(isinel(40) == 0)then
        si1=0D0
        goto 4681
        endif
        the=pmas(pycomp(311),1)+pmas(pycomp(3312),1)
        si1=s1724(ss,ilo1,0,the)
4681    if(isinel(260) == 0)then
        si2=0D0
        goto 4682
        endif
        the=pmas(pycomp(-321),1)+pmas(pycomp(2112),1)
        ik3=-321
        ik4=2112
        call srev(kl,kl1,ik3,ik4,ss,ilo2,fac, &
         1.,1.,0.,0.5,0.5,0.5,0.,0.5,1.)
        si2=s1724(ss,ilo2,0,the)*10D0*fac
4682    if(ilo1 == 0.and.ilo2 == 0) goto 13
        if(si1 < 1D-6.and.si2 < 1D-6) goto 13
        sit=si1+si2
        s1=si1/sit
        rlus=pyr(1)
        if(rlus <= s1)then
        ik1=311
        ik2=3312
        ic=40
!       pion0 + sigma- to k0+cascade-
        goto 2683
        endif
        ik1=-321
        ik2=2112
        ic=260
!       pion0 + sigma- to k- + n
2683    lc(icp,3)=ik1
        lc(icp,4)=ik2
        lc(icp,5)=ic
        tw(icp)=sit/cspin/10D0
        goto 10
        endif
!       pion0 + sigma-bar
        if(kl == 111 .and. kl1 == -3112)then
        if(isinel(47) == 0)then
        si1=0D0
        goto 7681
        endif
        the=pmas(pycomp(-311),1)+pmas(pycomp(-3312),1)
        si1=s1724(ss,ilo1,0,the)
7681    if(isinel(277) == 0)then
        si2=0D0
        goto 7682
        endif
        the=pmas(pycomp(321),1)+pmas(pycomp(-2112),1)
        ik3=321
        ik4=-2112
        call srev(kl,kl1,ik3,ik4,ss,ilo2,fac, &
         1.,1.,0.,0.5,0.5,0.5,0.,0.5,1.)
        si2=s1724(ss,ilo2,0,the)*10D0*fac
7682    if(ilo1 == 0.and.ilo2 == 0) goto 13
        if(si1 < 1D-6.and.si2 < 1D-6) goto 13
        sit=si1+si2
        s1=si1/sit
        rlus=pyr(1)
        if(rlus <= s1)then
        ik1=-311
        ik2=-3312
        ic=47
!       pion0 + sigma-bar to k0- +cascade-bar
        goto 7688
        endif
        ik1=321
        ik2=-2112
        ic=277
!       pion0 + sigma-bar to k+ + nbar
7688    lc(icp,3)=ik1
        lc(icp,4)=ik2
        lc(icp,5)=ic
        tw(icp)=sit/cspin/10D0
        goto 10
        endif

!       K- + p
        if(kl == -321 .and. kl1 == 2212)then
        if(isinel(49) == 0)then
        si1=0D0
        goto 407
        endif
        the=pmas(pycomp(111),1)+pmas(pycomp(3122),1)
        si1=s1724(ss,ilo1,0,the)
!       cross section of k + n to pion + y is assumed to be equal to
!        ten times of pion + n to k + y
407     if(isinel(50) == 0)then
        si2=0D0
        goto 408
        endif
        the=pmas(pycomp(111),1)+pmas(pycomp(3212),1)
        si2=s1724(ss,ilo2,0,the)
408     if(isinel(51) == 0)then
        si3=0D0
        goto 411
        endif
        the=pmas(pycomp(-211),1)+pmas(pycomp(3222),1)
        si3=s1724(ss,ilo3,0,the)
411     if(isinel(52) == 0)then
        si4=0D0
        goto 431
        endif
        the=pmas(pycomp(211),1)+pmas(pycomp(3112),1)
        si4=s1724(ss,ilo4,0,the)
431     if(isinel(53) == 0)then
        si5=0D0
        goto 412
        endif
        the=pmas(pycomp(321),1)+pmas(pycomp(3312),1)
        si5=s1724(ss,ilo5,0,the)/10D0
!       cross section of Kaon + n to Kaon + cascade is assumed to be
!        equal to pion + n to Kaon + y,but take the different of
!        threshold energy into account

412     if(isinel(125) == 0)then
        si6=0.
        goto 725
        endif
        the=pmas(pycomp(311),1)+pmas(pycomp(3322),1)
        si6=s1724(ss,ilo6,0,the)/10D0

725     if(ilo1 == 0 .and. ilo2 == 0 .and. ilo3 == 0 .and. ilo4 == 0 &
         .and. ilo5 == 0 .and. ilo6 == 0) goto 13
        if(si1 < 1D-6.and.si2 < 1D-6.and.si3 < 1D-6 &
         .and.si4 < 1D-6 .and. si5 <= 1D-6 .and.si6 <= 1D-6) goto 13
        si12=si1+si2
        si13=si12+si3
        si14=si13+si4
        si15=si14+si5
        si16=si15+si6
        s1=si1/si16
        s2=si12/si16
        s3=si13/si16
        s4=si14/si16
        s5=si15/si16
        rlu1=pyr(1)
        if(rlu1 <= s1)then
        ik1=111
        ik2=3122
        ic=49
!       K- + p to pion0 + lambda
        goto 414
        endif
        if(rlu1 > s1 .and. rlu1 <= s2)then
        ik1=111
        ik2=3212
        ic=50
!       K- + p to pion0 + sigma0
        goto 414
        endif
        if(rlu1 > s2 .and. rlu1 <= s3)then
        ik1=-211
        ik2=3222
        ic=51
!       K- + p to pion- + sigma+
        goto 414
        endif
        if(rlu1 > s3 .and. rlu1 <= s4)then
        ik1=211
        ik2=3112
        ic=52
!       K- + p to pion+ + sigma-
        goto 414
        endif
        if(rlu1 > s4 .and. rlu1 <= s5)then
        ik1=321
        ik2=3312
        ic=53
!       K- + p to k+ + cascade-
        goto 414
        endif
        ik1=311
        ik2=3322
        ic=125
!       K- + p to k0 + cascade0
414     lc(icp,3)=ik1
        lc(icp,4)=ik2
        lc(icp,5)=ic
        tw(icp)=si16/cskn
        goto 10
        endif

!       K- + n
        if(kl == -321 .and. kl1 == 2112)then
        if(isinel(54) == 0)then
        si1=0D0
        goto 415
        endif
        the=pmas(pycomp(-211),1)+pmas(pycomp(3212),1)
        si1=s1724(ss,ilo1,0,the)
415     if(isinel(55) == 0)then
        si2=0D0
        goto 416
        endif
        the=pmas(pycomp(-211),1)+pmas(pycomp(3122),1)
        si2=s1724(ss,ilo2,0,the)
416     if(isinel(56) == 0)then
        si3=0D0
        goto 417
        endif
        the=pmas(pycomp(111),1)+pmas(pycomp(3112),1)
        si3=s1724(ss,ilo3,0,the)
417     if(isinel(57) == 0)then
        si4=0D0
        goto 432
        endif
        the=pmas(pycomp(311),1)+pmas(pycomp(3312),1)
        si4=s1724(ss,ilo4,0,the)/10D0

432     if(ilo1 == 0.and.ilo2 == 0.and.ilo3 == 0.and.ilo4 == 0) goto 13
        if(si1 < 1D-6.and.si2 < 1D-6.and.si3 < 1D-6 &
         .and.si4 < 1D-6) goto 13
        si12=si1+si2
        si13=si12+si3
        si14=si13+si4
        s1=si1/si14
        s2=si12/si14
        s3=si13/si14
        rlus=pyr(1)
        if(rlus <= s1)then
        ik1=-211
        ik2=3212
        ic=54
!       K- + n to pion- + sigma0
        goto 418
        endif
        if(rlus > s1 .and. rlus <= s2)then
        ik1=-211
        ik2=3122
        ic=55
!       K- + n to pion- + lambda
        goto 418
        endif
        if(rlus > s2 .and. rlus <= s3)then
        ik1=111
        ik2=3112
        ic=56
!       K- + n to pion0 + sigma-
        goto 418
        endif
        ik1=311
        ik2=3312
        ic=57
!       K- + n to k0 + cascade-
418     lc(icp,3)=ik1
        lc(icp,4)=ik2
        lc(icp,5)=ic
        tw(icp)=si14/cskn
        goto 10
        endif

!       k0- + p
        if(kl == -311 .and. kl1 == 2212)then
        if(isinel(58) == 0)then
        si1=0D0
        goto 419
        endif
        the=pmas(pycomp(211),1)+pmas(pycomp(3122),1)
        si1=s1724(ss,ilo1,0,the)
419     if(isinel(59) == 0)then
        si2=0D0
        goto 420
        endif
        the=pmas(pycomp(211),1)+pmas(pycomp(3212),1)
        si2=s1724(ss,ilo2,0,the)
420     if(isinel(60) == 0)then
        si3=0D0
        goto 421
        endif
        the=pmas(pycomp(111),1)+pmas(pycomp(3222),1)
        si3=s1724(ss,ilo3,0,the)
421     if(isinel(126) == 0)then
        si4=0D0
        goto 726
        endif
        the=pmas(pycomp(321),1)+pmas(pycomp(3322),1)
        si4=s1724(ss,ilo4,0,the)/10D0

726     if(ilo1 == 0.and.ilo2 == 0.and.ilo3 == 0.and.ilo4 == 0) goto 13
        if(si1 < 1D-6.and.si2 < 1D-6.and.si3 < 1D-6.and.si4 < &
         1D-6) goto 13
        si12=si1+si2
        si13=si12+si3
        sit=si13+si4
        s1=si1/sit
        s2=si12/sit
        s3=si13/sit
        rlus=pyr(1)
        if(rlus <= s1)then
        ik1=211
        ik2=3122
        ic=58
!       k0- + p to pion+ + lambda
        goto 422
        endif
        if(rlus > s1 .and. rlus <= s2)then
        ik1=211
        ik2=3212
        ic=59
!       k0- + p to pion+ + sigma0
        goto 422
        endif
        if(rlus > s2 .and. rlus <= s3)then
        ik1=111
        ik2=3222
        ic=60
!       k0- + p to pion0 + sigma+
        goto 422
        endif
        ik1=321
        ik2=3322
        ic=126
!       k0- + p to k+ + cascade0
422     lc(icp,3)=ik1
        lc(icp,4)=ik2
        lc(icp,5)=ic
        tw(icp)=sit/cskn
        goto 10
        endif

!       k0- + n
        if(kl == -311 .and. kl1 == 2112)then
        if(isinel(61) == 0)then
        si1=0D0
        goto 423
        endif
        the=pmas(pycomp(211),1)+pmas(pycomp(3112),1)
        si1=s1724(ss,ilo1,0,the)
423     if(isinel(62) == 0)then
        si2=0D0
        goto 424
        endif
        the=pmas(pycomp(-211),1)+pmas(pycomp(3222),1)
        si2=s1724(ss,ilo2,0,the)
424     if(isinel(63) == 0)then
        si3=0D0
        goto 425
        endif
        the=pmas(pycomp(111),1)+pmas(pycomp(3212),1)
        si3=s1724(ss,ilo3,0,the)
425     if(isinel(64) == 0)then
        si4=0D0
        goto 426
        endif
        the=pmas(pycomp(111),1)+pmas(pycomp(3122),1)
        si4=s1724(ss,ilo4,0,the)
426     if(isinel(65) == 0)then
        si5=0D0
        goto 433
        endif
        the=pmas(pycomp(321),1)+pmas(pycomp(3312),1)
        si5=s1724(ss,ilo5,0,the)/10D0
433     if(isinel(127) == 0)then
        si6=0D0
        goto 729
        endif
        the=pmas(pycomp(311),1)+pmas(pycomp(3322),1)
        si6=s1724(ss,ilo6,0,the)/10D0
729     if(ilo1 == 0 .and. ilo2 == 0 .and. ilo3 == 0 .and. ilo4 == 0 &
         .and. ilo5 == 0 .and. ilo6 == 0) goto 13
        if(si1 < 1D-6.and.si2 < 1D-6.and.si3 < 1D-6 &
         .and.si4 < 1D-6.and.si5 < 1D-6.and.si6 < 1D-6) goto 13
        si12=si1+si2
        si13=si12+si3
        si14=si13+si4
        si15=si14+si5
        si16=si15+si6
        s1=si1/si16
        s2=si12/si16
        s3=si13/si16
        s4=si14/si16
        s5=si15/si16
        rlu1=pyr(1)
        if(rlu1 <= s1)then
        ik1=211
        ik2=3112
        ic=61
!       k0- + n to pion+ + sigma-
        goto 427
        endif
        if(rlu1 > s1 .and. rlu1 <= s2)then
        ik1=-211
        ik2=3222
        ic=62
!       k0- + n to pion- + sigma+
        goto 427
        endif
        if(rlu1 > s2 .and. rlu1 <= s3)then
        ik1=111
        ik2=3212
        ic=63
!       k0- + n to pion0 + sigma0
        goto 427
        endif
        if(rlu1 > s3 .and. rlu1 <= s4)then
        ik1=111
        ik2=3122
        ic=64
!       k0- + n to pion0 + lambda
        goto 427
        endif
        if(rlu1 > s4 .and. rlu1 <= s5)then
        ik1=321
        ik2=3312
        ic=65
!       k0- + n to k+ + cascade-
        goto 427
        endif
        ik1=311
        ik2=3322
        ic=127
!       k0- + n to k0 + cascade0
427     lc(icp,3)=ik1
        lc(icp,4)=ik2
        lc(icp,5)=ic
        tw(icp)=si16/cskn
        goto 10
        endif

!       K+ + p-
        if(kl == 321 .and. kl1 == -2212)then
        if(isinel(66) == 0)then
        si1=0D0
        goto 510
        endif
        the=pmas(pycomp(111),1)+pmas(pycomp(-3122),1)
        si1=s1724(ss,ilo1,0,the)
!       cross section of k + n- to pion + y- is assumed to be equal to
!        ten times of pion + n to k + y
510     if(isinel(67) == 0)then
        si2=0D0
        goto 511
        endif
        the=pmas(pycomp(111),1)+pmas(pycomp(-3212),1)
        si2=s1724(ss,ilo2,0,the)
511     if(isinel(68) == 0)then
        si3=0D0
        goto 512
        endif
        the=pmas(pycomp(211),1)+pmas(pycomp(-3222),1)
        si3=s1724(ss,ilo3,0,the)
512     if(isinel(69) == 0)then
        si4=0D0
        goto 513
        endif
        the=pmas(pycomp(-211),1)+pmas(pycomp(-3112),1)
        si4=s1724(ss,ilo4,0,the)
513     if(isinel(70) == 0)then
        si5=0D0
        goto 514
        endif
        the=pmas(pycomp(-321),1)+pmas(pycomp(-3312),1)
        si5=s1724(ss,ilo5,0,the)/10D0

514     if(isinel(128) == 0)then
        si6=0D0
        goto 730
        endif
        the=pmas(pycomp(-311),1)+pmas(pycomp(-3322),1)

        si6=s1724(ss,ilo6,0,the)/10D0

730     if(ilo1 == 0 .and. ilo2 == 0 .and. ilo3 == 0 .and. ilo4 == 0 &
         .and. ilo5 == 0 .and. ilo6 == 0) goto 13
        if(si1 < 1D-6.and.si2 < 1D-6.and.si3 < 1D-6 &
         .and.si4 < 1D-6 .and. si5 < 1D-6.and.si6 < 1D-6) goto 13
        si12=si1+si2
        si13=si12+si3
        si14=si13+si4
        si15=si14+si5
        si16=si15+si6
        s1=si1/si16
        s2=si12/si16
        s3=si13/si16
        s4=si14/si16
        s5=si15/si16
        rlu1=pyr(1)
        if(rlu1 <= s1)then
        ik1=111
        ik2=-3122
        ic=66
!       K+ + p- to pion0 + lambda-
        goto 515
        endif
        if(rlu1 > s1 .and. rlu1 <= s2)then
        ik1=111
        ik2=-3212
        ic=67
!       K+ + p- to pion0 + sigma0-
        goto 515
        endif
        if(rlu1 > s2 .and. rlu1 <= s3)then
        ik1=211
        ik2=-3222
        ic=68
!       K+ + p- to pion+ + sigma+ba
        goto 515
        endif
        if(rlu1 > s3 .and. rlu1 <= s4)then
        ik1=-211
        ik2=-3112
        ic=69
!       K+ + p- to pion- + sigma-ba
        goto 515
        endif
        if(rlu1 > s4 .and. rlu1 <= s5)then
        ik1=-321
        ik2=-3312
        ic=70
!       K+ + p- to k- + cascade-ba
        goto 515
        endif
        ik1=-311
        ik2=-3322
        ic=128
!       K+ + p- to k0- + cascade0ba
515     lc(icp,3)=ik1
        lc(icp,4)=ik2
        lc(icp,5)=ic
        tw(icp)=si16/cskn
        goto 10
        endif

!       K+ + n-
        if(kl == 321 .and. kl1 == -2112)then
        if(isinel(71) == 0)then
        si1=0D0
        goto 516
        endif
        the=pmas(pycomp(211),1)+pmas(pycomp(-3212),1)
        si1=s1724(ss,ilo1,0,the)
516     if(isinel(72) == 0)then
        si2=0D0
        goto 517
        endif
        the=pmas(pycomp(211),1)+pmas(pycomp(-3122),1)
        si2=s1724(ss,ilo2,0,the)
517     if(isinel(73) == 0)then
        si3=0D0
        goto 518
        endif
        the=pmas(pycomp(111),1)+pmas(pycomp(-3112),1)
        si3=s1724(ss,ilo3,0,the)
518     if(isinel(74) == 0)then
        si4=0D0
        goto 519
        endif
        the=pmas(pycomp(-311),1)+pmas(pycomp(-3312),1)
        si4=s1724(ss,ilo4,0,the)/10D0

519     if(ilo1 == 0.and.ilo2 == 0.and.ilo3 == 0.and.ilo4 == 0) goto 13
        if(si1 < 1D-6.and.si2 < 1D-6.and.si3 < 1D-6 &
         .and.si4 < 1D-6) goto 13
        si12=si1+si2
        si13=si12+si3
        si14=si13+si4
        s1=si1/si14
        s2=si12/si14
        s3=si13/si14
        rlus=pyr(1)
        if(rlus <= s1)then
        ik1=211
        ik2=-3212
        ic=71
!       K+ + n- to pion+ + sigma0-
        goto 520
        endif
        if(rlus > s1 .and. rlus <= s2)then
        ik1=211
        ik2=-3122
        ic=72
!       K+ + n- to pion+ + lambda-
        goto 520
        endif
        if(rlus > s2 .and. rlus <= s3)then
        ik1=111
        ik2=-3112
        ic=73
!       K+ + n- to pion0 + sigma-ba
        goto 520
        endif
        ik1=-311
        ik2=-3312
        ic=74
!       K+ + n- to k0- + cascade-ba
520     lc(icp,3)=ik1
        lc(icp,4)=ik2
        lc(icp,5)=ic
        tw(icp)=si14/cskn
        goto 10
        endif

!       k0 + p-
        if(kl == 311 .and. kl1 == -2212)then
        if(isinel(75) == 0)then
        si1=0D0
        goto 521
        endif
        the=pmas(pycomp(-211),1)+pmas(pycomp(-3122),1)
        si1=s1724(ss,ilo1,0,the)
521     if(isinel(76) == 0)then
        si2=0D0
        goto 522
        endif
        the=pmas(pycomp(-211),1)+pmas(pycomp(-3212),1)
        si2=s1724(ss,ilo2,0,the)
522     if(isinel(77) == 0)then
        si3=0D0
        goto 523
        endif
        the=pmas(pycomp(111),1)+pmas(pycomp(-3222),1)
        si3=s1724(ss,ilo3,0,the)
523     if(isinel(129) == 0)then
        si4=0D0
        goto 731
        endif
        the=pmas(pycomp(-321),1)+pmas(pycomp(-3322),1)
        si4=s1724(ss,ilo4,0,the)/10D0

731     if(ilo1 == 0.and.ilo2 == 0.and.ilo3 == 0.and.ilo4 == 0) goto 13
        if(si1 < 1D-6.and.si2 < 1D-6.and.si3 < 1D-6.and. &
         si4 < 1D-6) goto 13
        si12=si1+si2
        si13=si12+si3
        sit=si13+si4
        s1=si1/sit
        s2=si12/sit
        s3=si13/sit
        rlus=pyr(1)
        if(rlus <= s1)then
        ik1=-211
        ik2=-3122
        ic=75
!       k0 + p- to pion- + lambda-
        goto 524
        endif
        if(rlus > s1 .and. rlus <= s2)then
        ik1=-211
        ik2=-3212
        ic=76
!       k0 + p- to pion- + sigma0-
        goto 524
        endif
        if(rlus > s2 .and. rlus <= s3)then
        ik1=111
        ik2=-3222
        ic=77
!       k0 + p- to pion0 + sigma+ba
        goto 524
        endif
        ik1=-321
        ik2=-3322
        ic=129
!       k0 + p- to k- + cascade0-
524     lc(icp,3)=ik1
        lc(icp,4)=ik2
        lc(icp,5)=ic
        tw(icp)=sit/cskn
        goto 10
        endif

!       k0 + n-
        if(kl == 311 .and. kl1 == -2112)then
        if(isinel(78) == 0)then
        si1=0D0
        goto 525
        endif
        the=pmas(pycomp(-211),1)+pmas(pycomp(-3112),1)
        si1=s1724(ss,ilo1,0,the)
525     if(isinel(79) == 0)then
        si2=0D0
        goto 526
        endif
        the=pmas(pycomp(211),1)+pmas(pycomp(-3222),1)
        si2=s1724(ss,ilo2,0,the)
526     if(isinel(80) == 0)then
        si3=0D0
        goto 527
        endif
        the=pmas(pycomp(111),1)+pmas(pycomp(-3212),1)
        si3=s1724(ss,ilo3,0,the)
527     if(isinel(81) == 0)then
        si4=0D0
        goto 528
        endif
        the=pmas(pycomp(111),1)+pmas(pycomp(-3122),1)
        si4=s1724(ss,ilo4,0,the)
528     if(isinel(82) == 0)then
        si5=0D0
        goto 5299   ! 280224 Lei 529 -> 5299
        endif
        the=pmas(pycomp(-321),1)+pmas(pycomp(-3312),1)
        si5=s1724(ss,ilo5,0,the)/10D0

5299    continue   ! 280224 Lei
        if(isinel(130) == 0)then
        si6=0D0
        goto 529
        endif
        the=pmas(pycomp(-311),1)+pmas(pycomp(-3322),1)
        si6=s1724(ss,ilo6,0,the)/10D0

529     if(ilo1 == 0 .and. ilo2 == 0 .and. ilo3 == 0 .and. ilo4 == 0 &
         .and. ilo5 == 0.and. ilo6 == 0) goto 13
        if(si1 < 1D-6.and.si2 < 1D-6.and.si3 < 1D-6 &
         .and.si4 < 1D-6.and.si5 < 1D-6.and.si6 < 1D-6) goto 13
        si12=si1+si2
        si13=si12+si3
        si14=si13+si4
        si15=si14+si5
        si16=si15+si6
        s1=si1/si16
        s2=si12/si16
        s3=si13/si16
        s4=si14/si16
        s5=si15/si16
        rlu1=pyr(1)
        if(rlu1 <= s1)then
        ik1=-211
        ik2=-3112
        ic=78
!       k0 + n- to pion- + sigma-ba
        goto 530
        endif
        if(rlu1 > s1 .and. rlu1 <= s2)then
        ik1=211
        ik2=-3222
        ic=79
!       k0 + n- to pion+ + sigma+ba
        goto 530
        endif
        if(rlu1 > s2 .and. rlu1 <= s3)then
        ik1=111
        ik2=-3212
        ic=80
!       k0 + n- to pion0 + sigma0-
        goto 530
        endif
        if(rlu1 > s3 .and. rlu1 <= s4)then
        ik1=111
        ik2=-3122
        ic=81
!       k0 + n- to pion0 + lambda-
        goto 530
        endif
        if(rlu1 > s4 .and. rlu1 <= s5)then
        ik1=-321
        ik2=-3312
        ic=82
!       k0 + n- to k- + cascade-ba
        goto 530
        endif
        ik1=-311
        ik2=-3322
        ic=130
!       k0 + n- to k0- + cascade0-
530     lc(icp,3)=ik1
        lc(icp,4)=ik2
        lc(icp,5)=ic
        tw(icp)=si16/cskn
        goto 10
        endif

!       K- + lambda
        if(kl == -321 .and. kl1 == 3122)then
        if(isinel(83) == 0)then
        si1=0D0
        goto 732
        endif
        the=pmas(pycomp(111),1)+pmas(pycomp(3312),1)
        si1=s1724(ss,ilo1,0,the)
!       cross section of k + y to pion + cascade is assumed to be equal to
!        k + n to pion + y
732     if(isinel(113) == 0)then
        si2=0D0
        goto 733
        endif
        the=pmas(pycomp(-211),1)+pmas(pycomp(3322),1)
        si2=s1724(ss,ilo2,0,the)
733     if(ilo1 == 0.and.ilo2 == 0) goto 13
        if(si1 < 1D-6.and.si2 < 1D-6) goto 13
        sit=si1+si2
        s1=si1/sit
        rlus=pyr(1)
        if(rlus <= s1)then
        ik1=111
        ik2=3312
        ic=83
!       K- + lambda to pion0 + cascade-
        goto 734
        endif
        ik1=-211
        ik2=3322
        ic=113
!       K- + lambda to pion- + cascade0
734     lc(icp,3)=ik1
        lc(icp,4)=ik2
        lc(icp,5)=ic
        tw(icp)=sit/cskn
        goto 10
        endif

!       K- + sigma+
        if(kl == -321 .and. kl1 == 3222)then
        if(isinel(84) == 0)then
        si1=0D0
        goto 735
        endif
        the=pmas(pycomp(211),1)+pmas(pycomp(3312),1)
        si1=s1724(ss,ilo1,0,the)
735     if(isinel(111) == 0)then
        si2=0D0
        goto 736
        endif
        the=pmas(pycomp(111),1)+pmas(pycomp(3322),1)
        si2=s1724(ss,ilo2,0,the)
736     if(ilo1 == 0.and.ilo2 == 0) goto 13
        if(si1 < 1D-6.and.si2 < 1D-6) goto 13
        sit=si1+si2
        s1=si1/sit
        rlus=pyr(1)
        if(rlus <= s1)then
        ik1=211
        ik2=3312
        ic=84
!       K- + sigma+ to pion+ + cascade-
        goto 737
        endif
        ik1=111
        ik2=3322
        ic=111
!       K- + sigma+ to pion0 + cascade0
737     lc(icp,3)=ik1
        lc(icp,4)=ik2
        lc(icp,5)=ic
        tw(icp)=sit/cskn
        goto 10
        endif

!       K- + sigma- to pion- + cascade-
        if(kl == -321 .and. kl1 == 3112)then
        if(isinel(85) == 0) goto 13
        the=pmas(pycomp(-211),1)+pmas(pycomp(3312),1)
        tw(icp)=s1724(ss,ilo,0,the)/cskn
        if(ilo == 0) goto 13
        lc(icp,3)=-211
        lc(icp,4)=3312
        lc(icp,5)=85
        goto 10
        endif

!       K- + sigma0
        if(kl == -321 .and. kl1 == 3212)then
        if(isinel(86) == 0)then
        si1=0D0
        goto 738
        endif
        the=pmas(pycomp(111),1)+pmas(pycomp(3312),1)
        si1=s1724(ss,ilo1,0,the)
738     if(isinel(112) == 0)then
        si2=0D0
        goto 739
        endif
        the=pmas(pycomp(-211),1)+pmas(pycomp(3322),1)
        si2=s1724(ss,ilo2,0,the)
739     if(ilo1 == 0.and.ilo2 == 0) goto 13
        if(si1 < 1D-6.and.si2 < 1D-6) goto 13
        sit=si1+si2
        s1=si1/sit
        rlus=pyr(1)
        if(rlus <= s1)then
        ik1=111
        ik2=3312
        ic=86
!       K- + sigma0 to pion0 + cascade-
        goto 740
        endif
        ik1=-211
        ik2=3322
        ic=112
740     lc(icp,3)=ik1
        lc(icp,4)=ik2
        lc(icp,5)=ic
        tw(icp)=sit/cskn
        goto 10
        endif

!       k0- + lambda
        if(kl == -311 .and. kl1 == 3122)then
        if(isinel(87) == 0)then
        si1=0D0
        goto 741
        endif
        the=pmas(pycomp(211),1)+pmas(pycomp(3312),1)
        si1=s1724(ss,ilo1,0,the)
741     if(isinel(117) == 0)then
        si2=0D0
        goto 742
        endif
        the=pmas(pycomp(111),1)+pmas(pycomp(3322),1)
        si2=s1724(ss,ilo2,0,the)
742     if(ilo1 == 0.and.ilo2 == 0) goto 13
        if(si1 < 1D-6.and.si2 < 1D-6) goto 13
        sit=si1+si2
        s1=si1/sit
        rlus=pyr(1)
        if(rlus <= s1)then
        ik1=211
        ik2=3312
        ic=87
!       k0- + lambda to pion+ + cascade-
        goto 743
        endif
        ik1=111
        ik2=3322
        ic=117
!       k0- + lambda to pion0 + cascade0
743     lc(icp,3)=ik1
        lc(icp,4)=ik2
        lc(icp,5)=ic
        tw(icp)=sit/cskn
        goto 10
        endif

!       k0- + sigma0
        if(kl == -311 .and. kl1 == 3212)then
        if(isinel(88) == 0)then
        si1=0D0
        goto 744
        endif
        the=pmas(pycomp(211),1)+pmas(pycomp(3312),1)
        si1=s1724(ss,ilo1,0,the)
744     if(isinel(116) == 0)then
        si2=0D0
        goto 745
        endif
        the=pmas(pycomp(111),1)+pmas(pycomp(3322),1)
        si2=s1724(ss,ilo2,0,the)
745     if(ilo1 == 0.and.ilo2 == 0) goto 13
        if(si1 < 1D-6.and.si2 < 1D-6) goto 13
        sit=si1+si2
        s1=si1/sit
        rlus=pyr(1)
        if(rlus <= s1)then
        ik1=211
        ik2=3312
        ic=88
!       k0- + sigma0 to pion+ + cascade-
        goto 746
        endif
        ik1=111
        ik2=3322
        ic=116
!       k0- + sigma0 to pion0 + cascade0
746     lc(icp,3)=ik1
        lc(icp,4)=ik2
        lc(icp,5)=ic
        tw(icp)=sit/cskn
        goto 10
        endif

!       k0- + sigma-
        if(kl == -311 .and. kl1 == 3112)then
        if(isinel(89) == 0)then
        si1=0.
        goto 747
        endif
        the=pmas(pycomp(111),1)+pmas(pycomp(3312),1)
        si1=s1724(ss,ilo1,0,the)
747     if(isinel(115) == 0)then
        si2=0D0
        goto 748
        endif
        the=pmas(pycomp(-211),1)+pmas(pycomp(3322),1)
        si2=s1724(ss,ilo2,0,the)
748     if(ilo1 == 0.and.ilo2 == 0) goto 13
        if(si1 < 1D-6.and.si2 < 1D-6) goto 13
        sit=si1+si2
        s1=si1/sit
        rlus=pyr(1)
        if(rlus <= s1)then
        ik1=111
        ik2=3312
        ic=89
!       k0- + sigma- to pion0 + cascade-
        goto 749
        endif
        ik1=-211
        ik2=3322
        ic=115
!       k0- + sigma- to pion- + cascade0
749     lc(icp,3)=ik1
        lc(icp,4)=ik2
        lc(icp,5)=ic
        tw(icp)=sit/cskn
        goto 10
        endif

!       K+ + lambda-
        if(kl == 321 .and. kl1 == -3122)then
        if(isinel(90) == 0)then
        si1=0D0
        goto 750
        endif
        the=pmas(pycomp(111),1)+pmas(pycomp(-3312),1)
        si1=s1724(ss,ilo1,0,the)
750     if(isinel(120) == 0)then
        si2=0D0
        goto 751
        endif
        the=pmas(pycomp(211),1)+pmas(pycomp(-3322),1)
        si2=s1724(ss,ilo2,0,the)
751     if(ilo1 == 0.and.ilo2 == 0) goto 13
        if(si1 < 1D-6.and.si2 < 1D-6) goto 13
        sit=si1+si2
        s1=si1/sit
        rlus=pyr(1)
        if(rlus <= s1)then
        ik1=111
        ik2=-3312
        ic=90
!       K+ + lambda- to pion0 + cascade-ba
        goto 752
        endif
        ik1=211
        ik2=-3322
        ic=120
!       K+ + lambda- to pion+ + cascade0-
752     lc(icp,3)=ik1
        lc(icp,4)=ik2
        lc(icp,5)=ic
        tw(icp)=sit/cskn
        goto 10
        endif

!       K+ + sigma+ba
        if(kl == 321 .and. kl1 == -3222)then
        if(isinel(91) == 0)then
        si1=0D0
        goto 753
        endif
        the=pmas(pycomp(-211),1)+pmas(pycomp(-3312),1)
        si1=s1724(ss,ilo1,0,the)
753     if(isinel(118) == 0)then
        si2=0D0
        goto 754
        endif
        the=pmas(pycomp(111),1)+pmas(pycomp(-3322),1)
        si2=s1724(ss,ilo2,0,the)
754     if(ilo1 == 0.and.ilo2 == 0) goto 13
        if(si1 < 1D-6.and.si2 < 1D-6) goto 13
        sit=si1+si2
        s1=si1/sit
        rlus=pyr(1)
        if(rlus <= s1)then
        ik1=-211
        ik2=-3312
        ic=91
!       K+ + sigma+ba to pion- + cascade-ba
        goto 755
        endif
        ik1=111
        ik2=-3322
        ic=118
!       K+ + sigma+ba to pion0 + cascade0-
755     lc(icp,3)=ik1
        lc(icp,4)=ik2
        lc(icp,5)=ic
        tw(icp)=sit/cskn
        goto 10
        endif

!       K+ + sigma-ba to pion+ + cascade-ba
        if(kl == 321 .and. kl1 == -3112)then
        if(isinel(92) == 0) goto 13
        the=pmas(pycomp(211),1)+pmas(pycomp(-3312),1)
        tw(icp)=s1724(ss,ilo,0,the)/cskn
        if(ilo == 0) goto 13
        lc(icp,3)=211
        lc(icp,4)=-3312
        lc(icp,5)=92
        goto 10
        endif

!       K+ + sigma0ba
        if(kl == 321 .and. kl1 == -3212)then
        if(isinel(93) == 0)then
        si1=0D0
        goto 756
        endif
        the=pmas(pycomp(111),1)+pmas(pycomp(-3312),1)
        si1=s1724(ss,ilo1,0,the)
756     if(isinel(119) == 0)then
        si2=0D0
        goto 757
        endif
        the=pmas(pycomp(211),1)+pmas(pycomp(-3322),1)
        si2=s1724(ss,ilo2,0,the)
757     if(ilo1 == 0.and.ilo2 == 0) goto 13
        if(si1 < 1D-6.and.si2 < 1D-6) goto 13
        sit=si1+si2
        s1=si1/sit
        rlus=pyr(1)
        if(rlus <= s1)then
        ik1=111
        ik2=-3312
        ic=93
!       K+ + sigma0- to pion0 + cascade-ba
        goto 758
        endif
        ik1=211
        ik2=-3322
        ic=119
!       K+ + sigma0- to pion+ + cascade0-
758     lc(icp,3)=ik1
        lc(icp,4)=ik2
        lc(icp,5)=ic
        tw(icp)=sit/cskn
        goto 10
        endif

!       k0+ + lambda-
        if(kl == 311 .and. kl1 == -3122)then
        if(isinel(94) == 0)then
        si1=0D0
        goto 759
        endif
        the=pmas(pycomp(-211),1)+pmas(pycomp(-3312),1)
        si1=s1724(ss,ilo1,0,the)
759     if(isinel(124) == 0)then
        si2=0D0
        goto 760
        endif
        the=pmas(pycomp(111),1)+pmas(pycomp(-3322),1)
        si2=s1724(ss,ilo2,0,the)
760     if(ilo1 == 0.and.ilo2 == 0) goto 13
        if(si1 < 1D-6.and.si2 < 1D-6) goto 13
        sit=si1+si2
        s1=si1/sit
        rlus=pyr(1)
        if(rlus <= s1)then
        ik1=-211
        ik2=-3312
        ic=94
!       k0 + lambda- to pion- + cascade-ba
        goto 761
        endif
        ik1=111
        ik2=-3322
        ic=124
!       k0 + lambda- to pion0 + cascade0ba

761     lc(icp,3)=ik1
        lc(icp,4)=ik2
        lc(icp,5)=ic
        tw(icp)=sit/cskn
        goto 10
        endif

!       k0 + sigma0-
        if(kl == 311 .and. kl1 == -3212)then
        if(isinel(95) == 0)then
        si1=0D0
        goto 762
        endif
        the=pmas(pycomp(-211),1)+pmas(pycomp(-3312),1)
        si1=s1724(ss,ilo1,0,the)
762     if(isinel(123) == 0)then
        si2=0D0
        goto 763
        endif
        the=pmas(pycomp(111),1)+pmas(pycomp(-3322),1)
        si2=s1724(ss,ilo2,0,the)
763     if(ilo1 == 0.and.ilo2 == 0) goto 13
        if(si1 < 1D-6.and.si2 < 1D-6) goto 13
        sit=si1+si2
        s1=si1/sit
        rlus=pyr(1)
        if(rlus <= s1)then
        ik1=-211
        ik2=-3312
        ic=95
!       k0 + sigma0- to pion- + cascade-ba
        goto 764
        endif
        ik1=111
        ik2=-3322
        ic=123
!       k0 + sigma0- to pion0 + cascade0-
764     lc(icp,3)=ik1
        lc(icp,4)=ik2
        lc(icp,5)=ic
        tw(icp)=sit/cskn
        goto 10
        endif

!       k0 + sigma-ba
        if(kl == 311 .and. kl1 == -3112)then
        if(isinel(96) == 0)then
        si1=0D0
        goto 765
        endif
        the=pmas(pycomp(111),1)+pmas(pycomp(-3312),1)
        si1=s1724(ss,ilo1,0,the)
765     if(isinel(122) == 0)then
        si2=0D0
        goto 766
        endif
        the=pmas(pycomp(211),1)+pmas(pycomp(-3322),1)
        si2=s1724(ss,ilo2,0,the)
766     if(ilo1 == 0.and.ilo2 == 0) goto 13
        if(si1 < 1D-6.and.si2 < 1D-6) goto 13
        sit=si1+si2
        s1=si1/sit
        rlus=pyr(1)
        if(rlus <= s1)then
        ik1=111
        ik2=-3312
        ic=96
!       k0 + sigma-ba to pion0 + cascade-ba
        goto 767
        endif
        ik1=211
        ik2=-3322
        ic=122
!       k0 + sigma-ba to pion+ + cascade0ba
767     lc(icp,3)=ik1
        lc(icp,4)=ik2
        lc(icp,5)=ic
        tw(icp)=sit/cskn
        goto 10
        endif

!       k0- + sigma+ to pion+ + cascade0
        if(kl == -311 .and. kl1 == 3222)then
        if(isinel(114) == 0) goto 13
        the=pmas(pycomp(211),1)+pmas(pycomp(3322),1)
        tw(icp)=s1724(ss,ilo,0,the)/cskn
        if(ilo == 0) goto 13
        lc(icp,3)=211
        lc(icp,4)=3322
        lc(icp,5)=114
        goto 10
        endif

!       k0 + sigma+bar to pion- + cascade0ba
        if(kl == 311 .and. kl1 == -3222)then
        if(isinel(121) == 0) goto 13
        the=pmas(pycomp(-211),1)+pmas(pycomp(-3322),1)
        tw(icp)=s1724(ss,ilo,0,the)/cskn
        if(ilo == 0) goto 13
        lc(icp,3)=-211
        lc(icp,4)=-3322
        lc(icp,5)=121
        goto 10
        endif

!       K+ + cascade-ba to pion+ + Omega-ba
        if(kl == 321 .and. kl1 == -3312)then
        if(isinel(143) == 0) goto 13
        the=pmas(pycomp(211),1)+pmas(pycomp(-3334),1)
        tw(icp)=s1724(ss,ilo,0,the)/cskn
        if(ilo == 0) goto 13
        lc(icp,3)=211
        lc(icp,4)=-3334
        lc(icp,5)=143
        goto 10
        endif

!       K+ + cascade0- to pion0 + Omega-ba
        if(kl == 321 .and. kl1 == -3322)then
        if(isinel(145) == 0) goto 13
        the=pmas(pycomp(111),1)+pmas(pycomp(-3334),1)
        tw(icp)=s1724(ss,ilo,0,the)/cskn
        if(ilo == 0) goto 13
        lc(icp,3)=111
        lc(icp,4)=-3334
        lc(icp,5)=145
        goto 10
        endif

!       K- + cascade- to pion- + Omega-
        if(kl == -321 .and. kl1 == 3312)then
        if(isinel(139) == 0) goto 13
        the=pmas(pycomp(-211),1)+pmas(pycomp(3334),1)
        tw(icp)=s1724(ss,ilo,0,the)/cskn
        if(ilo == 0) goto 13
        lc(icp,3)=-211
        lc(icp,4)=3334
        lc(icp,5)=139
        goto 10
        endif

!       K- + cascade0 to pion0 + Omega-
        if(kl == -321 .and. kl1 == 3322)then
        if(isinel(141) == 0) goto 13
        the=pmas(pycomp(111),1)+pmas(pycomp(3334),1)
        tw(icp)=s1724(ss,ilo,0,the)/cskn
        if(ilo == 0) goto 13
        lc(icp,3)=111
        lc(icp,4)=3334
        lc(icp,5)=141
        goto 10
        endif

!       k0 + cascade-ba to pion0 + Omega-ba
        if(kl == 311 .and. kl1 == -3312)then
        if(isinel(144) == 0) goto 13
        the=pmas(pycomp(111),1)+pmas(pycomp(-3334),1)
        tw(icp)=s1724(ss,ilo,0,the)/cskn
        if(ilo == 0) goto 13
        lc(icp,3)=111
        lc(icp,4)=-3334
        lc(icp,5)=144
        goto 10
        endif

!       k0 + cascade0- to pion- + Omega-ba
        if(kl == 311 .and. kl1 == -3322)then
        if(isinel(146) == 0) goto 13
        the=pmas(pycomp(-211),1)+pmas(pycomp(-3334),1)
        tw(icp)=s1724(ss,ilo,0,the)/cskn
        if(ilo == 0) goto 13
        lc(icp,3)=-211
        lc(icp,4)=-3334
        lc(icp,5)=146
        goto 10
        endif

!       k0- + cascade- to pion0 + Omega-
        if(kl == -311 .and. kl1 == 3312)then
        if(isinel(140) == 0) goto 13
        the=pmas(pycomp(111),1)+pmas(pycomp(3334),1)
        tw(icp)=s1724(ss,ilo,0,the)/cskn
        if(ilo == 0) goto 13
        lc(icp,3)=111
        lc(icp,4)=3334
        lc(icp,5)=140
        goto 10
        endif

!       k0- + cascade0 to pion+ + Omega-
        if(kl == -311 .and. kl1 == 3322)then
        if(isinel(142) == 0) goto 13
        the=pmas(pycomp(211),1)+pmas(pycomp(3334),1)
        tw(icp)=s1724(ss,ilo,0,the)/cskn
        if(ilo == 0) goto 13
        lc(icp,3)=211
        lc(icp,4)=3334
        lc(icp,5)=142
        goto 10
        endif

!       follows are for reverse reactions
!       K+ + k-
        if((kl == 321 .and. kl1 == -321) &
         .or.(kl == -321 .and. kl1 == 321))then
        if(isinel(201) == 0)then
        fact1=0D0
        goto 601
        endif
        ik3=211
        ik4=-211
        ic3=201
        call srev(kl,kl1,ik3,ik4,ss,ilo1,fac, &
         0.5,0.5,0.,0.,1.,1.,0.,0.,1.)
!       K+ + k- to pi+ + pi-
        fact1=fac
        if(ilo1 == 0)fact1=0D0
601     if(isinel(202) == 0)then
        fact2=0D0
        goto 602
        endif
        ik5=111
        ik6=111
        ic5=202
        call srev(kl,kl1,ik5,ik6,ss,ilo2,fac, &
         0.5,0.5,0.,0.,1.,1.,0.,0.,0.5)
!       K+ + k- to pi0 + pi0
        fact2=fac
        if(ilo2 == 0)fact2=0D0
602     fact=fact1+fact2
        if( fact1 <= 1D-5 .AND. fact2 <= 1D-5 ) goto 13
!       if(ilo1 == 0 .and. ilo2 == 0) goto 13
        lc(icp,3)=ik3
        lc(icp,4)=ik4
        lc(icp,5)=ic3
        if(pyr(1) > fact1/fact)then
        lc(icp,3)=ik5
        lc(icp,4)=ik6
        lc(icp,5)=ic5
        endif
        tw(icp)=fact*cspipiKK/cspipi
        goto 10
        endif

!       K+ + k0-
        if((kl == 321 .and. kl1 == -311) &
         .or.(kl == -311 .and. kl1 == 321))then
        if(isinel(203) == 0) goto 13
        ik3=211
        ik4=111
        call srev(kl,kl1,ik3,ik4,ss,ilo1,fac, &
         0.5,0.5,0.,0.,1.,1.,0.,0.,1.)
!       K+ + k0- to pi+ + pi0
        if(ilo1 == 0) goto 13
        lc(icp,3)=ik3
        lc(icp,4)=ik4
        lc(icp,5)=203
        tw(icp)=fac*cspipiKK/cspipi
        goto 10
        endif

!       K- + k0
        if((kl == -321 .and. kl1 == 311) &
         .or.(kl == 311 .and. kl1 == -321))then
        if(isinel(204) == 0) goto 13
        ik3=-211
        ik4=111
        call srev(kl,kl1,ik3,ik4,ss,ilo1,fac, &
         0.5,0.5,0.,0.,1.,1.,0.,0.,1.)
!       K- + k0 to pi- + pi0
        if(ilo1 == 0) goto 13
        lc(icp,3)=ik3
        lc(icp,4)=ik4
        lc(icp,5)=204
        tw(icp)=fac*cspipiKK/cspipi
        goto 10
        endif

!       k0 + k0-
        if((kl == 311 .and. kl1 == -311) &
         .or.(kl == -311 .and. kl1 == 311))then
        if(isinel(205) == 0)then
        fact1=0D0
        goto 603
        endif
        ik3=211
        ik4=-211
        ic3=205
        call srev(kl,kl1,ik3,ik4,ss,ilo1,fac, &
         0.5,0.5,0.,0.,1.,1.,0.,0.,1.)
!       k0 + k0- to pi- + pi+
        fact1=fac
        if(ilo1 == 0)fact1=0D0
603     if(isinel(206) == 0)then
        fact2=0D0
        goto 604
        endif
        ik5=111
        ik6=111
        ic5=206
        call srev(kl,kl1,ik5,ik6,ss,ilo2,fac, &
         0.5,0.5,0.,0.,1.,1.,0.,0.,0.5)
!       k0 + k0- to pi0 + pi0
        fact2=fac
        if(ilo2 == 0)fact2=0D0
604     fact=fact1+fact2
        if( fact1 <= 1D-5 .AND. fact2 <= 1D-5 ) goto 13
!       if(ilo1 == 0 .and. ilo2 == 0) goto 13
        lc(icp,3)=ik3
        lc(icp,4)=ik4
        lc(icp,5)=ic3
        if(pyr(1) > fact1/fact)then
        lc(icp,3)=ik5
        lc(icp,4)=ik6
        lc(icp,5)=ic5
        endif
        tw(icp)=fact*cspipiKK/cspipi
        goto 10
        endif

!       K+ + sigma+ to pion+ + p
        if(kl == 321 .and. kl1 == 3222)then
!       3222 is the flavor code of sigma+
        if(isinel(207) == 0) goto 13
        the=pmas(pycomp(211),1)+pmas(pycomp(2212),1)
        ww=s1713(ss,ilo,0,the)/cskn/10D0
        if(ilo == 0) goto 13
        lc(icp,3)=211
        lc(icp,4)=2212
        lc(icp,5)=207
        ik3=211
        ik4=2212
        call srev(kl,kl1,ik3,ik4,ss,ilo2,fac, &
         0.5,1.,0.,0.5,1.,0.5,0.,0.5,1.)
        tw(icp)=ww*fac
        goto 10
        endif

!       K+ + sigma-
        if(kl == 321 .and. kl1 == 3112)then
        if(isinel(208) == 0)then
        si1=0D0
        goto 606
        endif
        ik1=-211
        ik2=2212
        call srev(kl,kl1,ik1,ik2,ss,ilo1,fac, &
         0.5,1.,0.,0.5,1.,0.5,0.,0.5,1.)
        the=pmas(pycomp(-211),1)+pmas(pycomp(2212),1)
        si1=s0715(ss,ilo1,0,the)*fac
606     if(isinel(209) == 0)then
        si2=0D0
        goto 607
        endif
        ik1=111
        ik2=2112
        call srev(kl,kl1,ik1,ik2,ss,ilo2,fac, &
         0.5,1.,0.,0.5,1.,0.5,0.,0.5,1.)
        the=pmas(pycomp(111),1)+pmas(pycomp(2112),1)
        si2=s2325(ss,ilo2,0,the)*fac
607     if(ilo1 == 0.and.ilo2 == 0) goto 13
        if(si1 < 1D-6.and.si2 < 1D-6) goto 13
        sit=si1+si2
        s1=si1/sit
        rlus=pyr(1)
        if(rlus <= s1)then
        ik1=-211
        ik2=2212
        ic=208
!       K+ + sigma- to pion- + p
        goto 609
        endif
        ik1=111
        ik2=2112
        ic=209
!       K+ + sigma- to pion0 + n
609     lc(icp,3)=ik1
        lc(icp,4)=ik2
        lc(icp,5)=ic
        tw(icp)=sit/cskn/10D0
        goto 10
        endif

!       K+ + sigma0
        if(kl == 321 .and. kl1 == 3212)then
        if(isinel(210) == 0)then
        si1=0D0
        goto 610
        endif
        ik1=211
        ik2=2112
        call srev(kl,kl1,ik1,ik2,ss,ilo1,fac, &
         0.5,1.,0.,0.5,1.,0.5,0.,0.5,1.)
        the=pmas(pycomp(211),1)+pmas(pycomp(2112),1)
        si1=s1724(ss,ilo1,0,the)*fac
610     if(isinel(211) == 0)then
        si2=0D0
        goto 611
        endif
        ik1=111
        ik2=2212
        call srev(kl,kl1,ik1,ik2,ss,ilo2,fac, &
         0.5,1.,0.,0.5,1.,0.5,0.,0.5,1.)
        the=pmas(pycomp(111),1)+pmas(pycomp(2212),1)
        si2=s2314(ss,ilo2,0,the)*fac
611     if(ilo1 == 0.and.ilo2 == 0) goto 13
        if(si1 < 1D-6.and.si2 < 1D-6) goto 13
        sit=si1+si2
        s1=si1/sit
        rlus=pyr(1)
        if(rlus <= s1)then
        ik1=211
        ik2=2112
        ic=210
!       K+ + sigma0 to pion+ + n
        goto 612
        endif
        ik1=111
        ik2=2212
        ic=211
!       K+ + sigma0 to pion0 + p
612     lc(icp,3)=ik1
        lc(icp,4)=ik2
        lc(icp,5)=ic
        tw(icp)=sit/cskn/10D0
        goto 10
        endif

!       K+ + lambda0
        if(kl == 321 .and. kl1 == 3122)then
        if(isinel(212) == 0)then
        si1=0D0
        goto 613
        endif
        ik1=211
        ik2=2112
        call srev(kl,kl1,ik1,ik2,ss,ilo1,fac, &
         0.5,0.,0.,0.5,1.,0.5,0.,0.5,1.)
        the=pmas(pycomp(211),1)+pmas(pycomp(2112),1)
        si1=s1727(ss,ilo1,0,the)*fac
613     if(isinel(213) == 0)then
        si2=0D0
        goto 614
        endif
        ik1=111
        ik2=2212
        call srev(kl,kl1,ik1,ik2,ss,ilo2,fac, &
         0.5,0.,0.,0.5,1.,0.5,0.,0.5,1.)
        the=pmas(pycomp(111),1)+pmas(pycomp(2212),1)
        si2=s2317(ss,ilo2,0,the)*fac
614     if(ilo1 == 0.and.ilo2 == 0) goto 13
        if(si1 < 1D-6.and.si2 < 1D-6) goto 13
        sit=si1+si2
        s1=si1/sit
        rlus=pyr(1)
        if(rlus <= s1)then
        ik1=211
        ik2=2112
        ic=212
!       K+ + lambda0 to pion+ + n
        goto 615
        endif
        ik1=111
        ik2=2212
        ic=213
!       K+ + lambda0 to pion0 + p
615     lc(icp,3)=ik1
        lc(icp,4)=ik2
        lc(icp,5)=ic
        tw(icp)=sit/cskn/10D0
        goto 10
        endif

!       k0 + sigma+
        if(kl == 311 .and. kl1 == 3222)then
        if(isinel(214) == 0)then
        si1=0D0
        goto 616
        endif
        ik1=211
        ik2=2112
        call srev(kl,kl1,ik1,ik2,ss,ilo1,fac, &
         0.5,1.,0.,0.5,1.,0.5,0.,0.5,1.)
        the=pmas(pycomp(211),1)+pmas(pycomp(2112),1)
        si1=s1724(ss,ilo1,0,the)*fac
616     if(isinel(215) == 0)then
        si2=0D0
        goto 617
        endif
        ik1=111
        ik2=2212
        call srev(kl,kl1,ik1,ik2,ss,ilo2,fac, &
         0.5,1.,0.,0.5,1.,0.5,0.,0.5,1.)
        the=pmas(pycomp(111),1)+pmas(pycomp(2212),1)
        si2=s1724(ss,ilo2,0,the)*fac
617     if(ilo1 == 0.and.ilo2 == 0) goto 13
        if(si1 < 1D-6.and.si2 < 1D-6) goto 13
        sit=si1+si2
        s1=si1/sit
        rlus=pyr(1)
        if(rlus <= s1)then
        ik1=211
        ik2=2112
        ic=214
!       K+ + sigma+ to pion+ + n
        goto 618
        endif
        ik1=111
        ik2=2212
        ic=215
!       K+ + sigma+ to pion0 + p
618     lc(icp,3)=ik1
        lc(icp,4)=ik2
        lc(icp,5)=ic
        tw(icp)=sit/cskn/10D0
        goto 10
        endif

!       k0 + sigma- to pion- +n
        if(kl == 311 .and. kl1 == 3112)then
!       3112 is the flavor code of sigma-
        if(isinel(216) == 0) goto 13
        the=pmas(pycomp(-211),1)+pmas(pycomp(2112),1)
        ww=s1724(ss,ilo,0,the)/cskn/10D0
        if(ilo == 0) goto 13
        lc(icp,3)=-211
        lc(icp,4)=2112
        lc(icp,5)=216
        ik3=-211
        ik4=2112
        call srev(kl,kl1,ik3,ik4,ss,ilo2,fac, &
         0.5,1.,0.,0.5,1.,0.5,0.,0.5,1.)
        tw(icp)=ww*fac
        goto 10
        endif

!       k0 + sigma0
        if(kl == 311 .and. kl1 == 3212)then
        if(isinel(217) == 0)then
        si1=0D0
        goto 619
        endif
        ik1=-211
        ik2=2212
        call srev(kl,kl1,ik1,ik2,ss,ilo1,fac, &
         0.5,1.,0.,0.5,1.,0.5,0.,0.5,1.)
        the=pmas(pycomp(-211),1)+pmas(pycomp(2212),1)
        si1=s07123(ss,ilo1,0,the)*fac
619     if(isinel(218) == 0)then
        si2=0D0
        goto 620
        endif
        ik1=111
        ik2=2112
        call srev(kl,kl1,ik1,ik2,ss,ilo2,fac, &
         0.5,1.,0.,0.5,1.,0.5,0.,0.5,1.)
        the=pmas(pycomp(111),1)+pmas(pycomp(2112),1)
        si2=s1724(ss,ilo2,0,the)*fac
620     if(ilo1 == 0.and.ilo2 == 0) goto 13
        if(si1 < 1D-6.and.si2 < 1D-6) goto 13
        sit=si1+si2
        s1=si1/sit
        rlus=pyr(1)
        if(rlus <= s1)then
        ik1=-211
        ik2=2212
        ic=217
!       k0 + sigma0 to pion- + p
        goto 621
        endif
        ik1=111
        ik2=2112
        ic=218
!       k0 + sigma0 to pion0 + n
621     lc(icp,3)=ik1
        lc(icp,4)=ik2
        lc(icp,5)=ic
        tw(icp)=sit/cskn/10D0
        goto 10
        endif

!       k0 + lambda0
        if(kl == 311 .and. kl1 == 3122)then
        if(isinel(219) == 0)then
        si1=0D0
        goto 622
        endif
        ik1=-211
        ik2=2212
        call srev(kl,kl1,ik1,ik2,ss,ilo1,fac, &
         0.5,0.,0.,0.5,1.,0.5,0.,0.5,1.)
        the=pmas(pycomp(-211),1)+pmas(pycomp(2212),1)
        si1=s07122(ss,ilo1,0,the)*fac
622     if(isinel(220) == 0)then
        si2=0D0
        goto 623
        endif
        ik1=111
        ik2=2112
        call srev(kl,kl1,ik1,ik2,ss,ilo2,fac, &
         0.5,0.,0.,0.5,1.,0.5,0.,0.5,1.)
        the=pmas(pycomp(111),1)+pmas(pycomp(2112),1)
        si2=s1724(ss,ilo2,0,the)*fac
623     if(ilo1 == 0.and.ilo2 == 0) goto 13
        if(si1 < 1D-6.and.si2 < 1D-6) goto 13
        sit=si1+si2
        s1=si1/sit
        rlus=pyr(1)
        if(rlus <= s1)then
        ik1=-211
        ik2=2212
        ic=219
!       k0 + lambda0 to pion- + p
        goto 624
        endif
        ik1=111
        ik2=2112
        ic=220
!       k0 + lambda0 to pion0 + n
624     lc(icp,3)=ik1
        lc(icp,4)=ik2
        lc(icp,5)=ic
        tw(icp)=sit/cskn/10D0
        goto 10
        endif

!       K- + sigma+ba to pion- + pba
        if(kl == -321 .and. kl1 == -3222)then
        if(isinel(221) == 0) goto 13
        the=pmas(pycomp(-211),1)+pmas(pycomp(-2212),1)
        ww=s1724(ss,ilo,0,the)/cskn/10D0
        if(ilo == 0) goto 13
        lc(icp,3)=-211
        lc(icp,4)=-2212
        lc(icp,5)=221
        ik3=-211
        ik4=-2212
        call srev(kl,kl1,ik3,ik4,ss,ilo1,fac, &
         0.5,1.,0.,0.5,1.,0.5,0.,0.5,1.)
        tw(icp)=ww*fac
        goto 10
        endif

!       K- + sigma-ba
        if(kl == -321 .and. kl1 == -3112)then
        if(isinel(222) == 0)then
        si1=0D0
        goto 625
        endif
        ik1=211
        ik2=-2212
        call srev(kl,kl1,ik1,ik2,ss,ilo1,fac, &
         0.5,1.,0.,0.5,1.,0.5,0.,0.5,1.)
        the=pmas(pycomp(211),1)+pmas(pycomp(-2212),1)
        si1=s1724(ss,ilo1,0,the)*fac
625     if(isinel(223) == 0)then
        si2=0D0
        goto 626
        endif
        ik1=111
        ik2=-2112
        call srev(kl,kl1,ik1,ik2,ss,ilo2,fac, &
         0.5,1.,0.,0.5,1.,0.5,0.,0.5,1.)
        the=pmas(pycomp(111),1)+pmas(pycomp(-2112),1)
        si2=s1724(ss,ilo2,0,the)*fac
626     if(ilo1 == 0.and.ilo2 == 0) goto 13
        if(si1 < 1D-6.and.si2 < 1D-6) goto 13
        sit=si1+si2
        s1=si1/sit
        rlus=pyr(1)
        if(rlus <= s1)then
        ik1=211
        ik2=-2212
        ic=222
!       K- + sigma-ba to pion+ + pba
        goto 627
        endif
        ik1=111
        ik2=-2112
        ic=223
!       K- + sigma-ba to pion0 + nba
627     lc(icp,3)=ik1
        lc(icp,4)=ik2
        lc(icp,5)=ic
        tw(icp)=sit/cskn/10D0
        goto 10
        endif

!       K- + sigma0ba
        if(kl == -321 .and. kl1 == -3212)then
        if(isinel(224) == 0)then
        si1=0D0
        goto 628
        endif
        ik1=-211
        ik2=-2112
        call srev(kl,kl1,ik1,ik2,ss,ilo1,fac, &
         0.5,1.,0.,0.5,1.,0.5,0.,0.5,1.)
        the=pmas(pycomp(-211),1)+pmas(pycomp(-2112),1)
        si1=s1724(ss,ilo1,0,the)*fac
628     if(isinel(225) == 0)then
        si2=0D0
        goto 629
        endif
        ik1=111
        ik2=-2212
        call srev(kl,kl1,ik1,ik2,ss,ilo2,fac, &
         0.5,1.,0.,0.5,1.,0.5,0.,0.5,1.)
        the=pmas(pycomp(111),1)+pmas(pycomp(-2212),1)
        si2=s1724(ss,ilo2,0,the)*fac
629     if(ilo1 == 0.and.ilo2 == 0) goto 13
        if(si1 < 1D-6.and.si2 < 1D-6) goto 13
        sit=si1+si2
        s1=si1/sit
        rlus=pyr(1)
        if(rlus <= s1)then
        ik1=-211
        ik2=-2112
        ic=224
!       K- + sigma0ba to pion- + nba
        goto 630
        endif
        ik1=111
        ik2=-2212
        ic=225
!       K- + sigma0ba to pion0 + pba
630     lc(icp,3)=ik1
        lc(icp,4)=ik2
        lc(icp,5)=ic
        tw(icp)=sit/cskn/10D0
        goto 10
        endif

!       K- + lambda0-
        if(kl == -321 .and. kl1 == -3122)then
        if(isinel(226) == 0)then
        si1=0D0
        goto 631
        endif
        ik1=-211
        ik2=-2112
        call srev(kl,kl1,ik1,ik2,ss,ilo1,fac, &
         0.5,0.,0.,0.5,1.,0.5,0.,0.5,1.)
        the=pmas(pycomp(-211),1)+pmas(pycomp(-2112),1)
        si1=s1724(ss,ilo1,0,the)*fac
631     if(isinel(227) == 0)then
        si2=0D0
        goto 632
        endif
        ik1=111
        ik2=-2212
        call srev(kl,kl1,ik1,ik2,ss,ilo2,fac, &
         0.5,0.,0.,0.5,1.,0.5,0.,0.5,1.)
        the=pmas(pycomp(111),1)+pmas(pycomp(-2212),1)
        si2=s1724(ss,ilo2,0,the)*fac
632     if(ilo1 == 0.and.ilo2 == 0) goto 13
        if(si1 < 1D-6.and.si2 < 1D-6) goto 13
        sit=si1+si2
        s1=si1/sit
        rlus=pyr(1)
        if(rlus <= s1)then
        ik1=-211
        ik2=-2112
        ic=226
!       K- + lambda0ba to pion- + nba
        goto 633
        endif
        ik1=111
        ik2=-2212
        ic=227
!       K- + lambda0ba to pion0 + pba
633     lc(icp,3)=ik1
        lc(icp,4)=ik2
        lc(icp,5)=ic
        tw(icp)=sit/cskn/10D0
        goto 10
        endif

!       k0- + sigma+ba
        if(kl == -311 .and. kl1 == -3222)then
        if(isinel(228) == 0)then
        si1=0D0
        goto 634
        endif
        ik1=-211
        ik2=-2112
        call srev(kl,kl1,ik1,ik2,ss,ilo1,fac, &
         0.5,1.,0.,0.5,1.,0.5,0.,0.5,1.)
        the=pmas(pycomp(-211),1)+pmas(pycomp(-2112),1)
        si1=s1724(ss,ilo1,0,the)*fac
634     if(isinel(229) == 0)then
        si2=0D0
        goto 635
        endif
        ik1=111
        ik2=-2212
        call srev(kl,kl1,ik1,ik2,ss,ilo2,fac, &
         0.5,1.,0.,0.5,1.,0.5,0.,0.5,1.)
        the=pmas(pycomp(111),1)+pmas(pycomp(-2212),1)
        si2=s1724(ss,ilo2,0,the)*fac
635     if(ilo1 == 0.and.ilo2 == 0) goto 13
        if(si1 < 1D-6.and.si2 < 1D-6) goto 13
        sit=si1+si2
        s1=si1/sit
        rlus=pyr(1)
        if(rlus <= s1)then
        ik1=-211
        ik2=-2112
        ic=228
!       k0ba + sigma+ba to pion- + nba
        goto 636
        endif
        ik1=111
        ik2=-2212
        ic=229
!       k0ba + sigma+ba to pion0 + pba
636     lc(icp,3)=ik1
        lc(icp,4)=ik2
        lc(icp,5)=ic
        tw(icp)=sit/cskn/10D0
        goto 10
        endif

!       k0- + sigma-ba to pion+ + nba
        if(kl == -311 .and. kl1 == -3112)then
        if(isinel(230) == 0) goto 13
        the=pmas(pycomp(211),1)+pmas(pycomp(-2112),1)
        ww=s1724(ss,ilo,0,the)/cskn/10D0
        if(ilo == 0) goto 13
        lc(icp,3)=211
        lc(icp,4)=-2112
        lc(icp,5)=230
        ik3=211
        ik4=-2112
        call srev(kl,kl1,ik3,ik4,ss,ilo1,fac, &
         0.5,1.,0.,0.5,1.,0.5,0.,0.5,1.)
        tw(icp)=ww*fac
        goto 10
        endif

!       k0- + sigma0-
        if(kl == -311 .and. kl1 == -3212)then
        if(isinel(231) == 0)then
        si1=0D0
        goto 637
        endif
        ik1=211
        ik2=-2212
        call srev(kl,kl1,ik1,ik2,ss,ilo1,fac, &
         0.5,1.,0.,0.5,1.,0.5,0.,0.5,1.)
        the=pmas(pycomp(211),1)+pmas(pycomp(-2212),1)
        si1=s1724(ss,ilo1,0,the)*fac
637     if(isinel(232) == 0)then
        si2=0D0
        goto 638
        endif
        ik1=111
        ik2=-2112
        call srev(kl,kl1,ik1,ik2,ss,ilo2,fac, &
         0.5,1.,0.,0.5,1.,0.5,0.,0.5,1.)
        the=pmas(pycomp(111),1)+pmas(pycomp(-2112),1)
        si2=s1724(ss,ilo2,0,the)*fac
638     if(ilo1 == 0.and.ilo2 == 0) goto 13
        if(si1 < 1D-6.and.si2 < 1D-6) goto 13
        sit=si1+si2
        s1=si1/sit
        rlus=pyr(1)
        if(rlus <= s1)then
        ik1=211
        ik2=-2212
        ic=231
!       k0- + sigma0- to pion+ + pba
        goto 639
        endif
        ik1=111
        ik2=-2112
        ic=232
!       k0- + sigma0- to pion0 + nba
639     lc(icp,3)=ik1
        lc(icp,4)=ik2
        lc(icp,5)=ic
        tw(icp)=sit/cskn/10D0
        goto 10
        endif

!       k0- + lambda0-
        if(kl == -311 .and. kl1 == -3122)then
        if(isinel(233) == 0)then
        si1=0D0
        goto 640
        endif
        ik1=211
        ik2=-2212
        call srev(kl,kl1,ik1,ik2,ss,ilo1,fac, &
         0.5,0.,0.,0.5,1.,0.5,0.,0.5,1.)
        the=pmas(pycomp(211),1)+pmas(pycomp(-2212),1)
        si1=s1724(ss,ilo1,0,the)*fac
640     if(isinel(234) == 0)then
        si2=0D0
        goto 641
        endif
        ik1=111
        ik2=-2112
        call srev(kl,kl1,ik1,ik2,ss,ilo2,fac, &
         0.5,0.,0.,0.5,1.,0.5,0.,0.5,1.)
        the=pmas(pycomp(111),1)+pmas(pycomp(-2112),1)
        si2=s1724(ss,ilo2,0,the)*fac
641     if(ilo1 == 0.and.ilo2 == 0) goto 13
        if(si1 < 1D-6.and.si2 < 1D-6) goto 13
        sit=si1+si2
        s1=si1/sit
        rlus=pyr(1)
        if(rlus <= s1)then
        ik1=211
        ik2=-2212
        ic=233
!       k0- + lambda0- to pion+ + pba
        goto 642
        endif
        ik1=111
        ik2=-2112
        ic=234
!       k0- + lambda0- to pion0 + nba
642     lc(icp,3)=ik1
        lc(icp,4)=ik2
        lc(icp,5)=ic
        tw(icp)=sit/cskn/10D0
        goto 10
        endif

!       K+ + cascade-
        if(kl == 321 .and. kl1 == 3312)then
        if(isinel(235) == 0)then
        si1=0D0
        goto 699
        endif
        ik1=211
        ik2=3112
        call srev(kl,kl1,ik1,ik2,ss,ilo1,fac, &
         0.5,0.5,0.,0.5,1.,1.,0.,0.5,1.)
        the=pmas(pycomp(211),1)+pmas(pycomp(3112),1)
        si1=s1724(ss,ilo1,0,the)*fac
699     if(isinel(236) == 0)then
        si2=0D0
        goto 643
        endif
        ik1=-211
        ik2=3222
        call srev(kl,kl1,ik1,ik2,ss,ilo2,fac, &
         0.5,0.5,0.,0.5,1.,1.,0.,0.5,1.)
        the=pmas(pycomp(-211),1)+pmas(pycomp(3222),1)
        si2=s1724(ss,ilo2,0,the)*fac
643     if(isinel(237) == 0)then
        si3=0D0
        goto 644
        endif
        ik1=111
        ik2=3122
        call srev(kl,kl1,ik1,ik2,ss,ilo3,fac, &
         0.5,0.5,0.,0.5,1.,0.,0.,0.5,1.)
        the=pmas(pycomp(111),1)+pmas(pycomp(3122),1)
        si3=s1724(ss,ilo3,0,the)*fac
644     if(isinel(238) == 0)then
        si4=0D0
        goto 6430
        endif
        ik1=111
        ik2=3212
        call srev(kl,kl1,ik1,ik2,ss,ilo4,fac, &
        0.5,0.5,0.,0.5,1.,1.,0.,0.5,1.)
        the=pmas(pycomp(111),1)+pmas(pycomp(3212),1)
        si4=s1724(ss,ilo4,0,the)*fac
6430    if(isinel(253) == 0)then
        si5=0D0
        goto 6440
        endif
        ik1=-321
        ik2=2212
        call srev(kl,kl1,ik1,ik2,ss,ilo5,fac, &
        0.5,0.5,0.,0.5,1.,0.,0.,0.5,1.)
        the=pmas(pycomp(321),1)+pmas(pycomp(2212),1)
        si5=s1724(ss,ilo5,0,the)*fac
6440    if(isinel(265) == 0)then
        si6=0D0
        goto 645
        endif
        ik1=-311
        ik2=2112
        call srev(kl,kl1,ik1,ik2,ss,ilo6,fac, &
        0.5,0.5,0.,0.5,1.,0.,0.,0.5,1.)
        the=pmas(pycomp(311),1)+pmas(pycomp(2112),1)
        si6=s1724(ss,ilo6,0,the)*fac
645     if(ilo1 == 0 .and. ilo2 == 0 .and. ilo3 == 0 .and. ilo4 == 0 &
         .and. ilo5 == 0 .and. ilo6 == 0) goto 13
        if(si1 < 1D-6.and.si2 < 1D-6.and.si3 < 1D-6 .and. &
         si4 < 1D-6.and.si5 < 1D-6.and.si6 < 1D-6) goto 13
        si12=si1+si2
        si13=si12+si3
        si14=si13+si4
        si15=si14+si5
        si16=si15+si6
        s1=si1/si16
        s2=si12/si16
        s3=si13/si16
        s4=si14/si16
        s5=si15/si16
        rlu1=pyr(1)
        if(rlu1 <= s1)then
        ik1=211
        ik2=3112
        ic=235
!       K+ + cascade- to pion+ + sigma-
        goto 646
        endif
        if(rlu1 > s1 .and. rlu1 <= s2)then
        ik1=-211
        ik2=3222
        ic=236
!       K+ + cascade- to pion- + sigma+
        goto 646
        endif
        if(rlu1 > s2 .and. rlu1 <= s3)then
        ik1=111
        ik2=3122
        ic=237
!       K+ + cascade- to pion0 + lambda0
        goto 646
        endif
        if(rlu1 > s3 .and. rlu1 <= s4)then
        ik1=111
        ik2=3212
        ic=238
!       K+ + cascade- to pion0 + sigma0
        goto 646
        endif
        if(rlu1 > s4 .and. rlu1 <= s5)then
        ik1=-321
        ik2=2212
        ic=253
!       K+ + cascade- to k- + p
        goto 646
        endif
        ik1=-311
        ik2=2112
        ic=265
!       K+ + cascade- to k0- + n
646     lc(icp,3)=ik1
        lc(icp,4)=ik2
        lc(icp,5)=ic
        tw(icp)=si16/cskn/10D0
        goto 10
        endif

!       K- + cascade-ba
        if(kl == -321 .and. kl1 == -3312)then
        if(isinel(239) == 0)then
        si1=0D0
        goto 647
        endif
        ik1=211
        ik2=-3222
        call srev(kl,kl1,ik1,ik2,ss,ilo1,fac, &
         0.5,0.5,0.,0.5,1.,1.,0.,0.5,1.)
        the=pmas(pycomp(211),1)+pmas(pycomp(-3222),1)
        si1=s1724(ss,ilo1,0,the)*fac
647     if(isinel(240) == 0)then
        si2=0D0
        goto 648
        endif
        ik1=-211
        ik2=-3112
        call srev(kl,kl1,ik1,ik2,ss,ilo2,fac, &
         0.5,0.5,0.,0.5,1.,1.,0.,0.5,1.)
        the=pmas(pycomp(-211),1)+pmas(pycomp(-3112),1)
        si2=s1724(ss,ilo2,0,the)*fac
648     if(isinel(241) == 0)then
        si3=0D0
        goto 649
        endif
        ik1=111
        ik2=-3122
        call srev(kl,kl1,ik1,ik2,ss,ilo3,fac, &
         0.5,0.5,0.,0.5,1.,0.,0.,0.5,1.)
        the=pmas(pycomp(111),1)+pmas(pycomp(-3122),1)
        si3=s1724(ss,ilo3,0,the)*fac
649     if(isinel(242) == 0)then
        si4=0D0
        goto 6500
        endif
        ik1=111
        ik2=-3212
        call srev(kl,kl1,ik1,ik2,ss,ilo4,fac, &
        0.5,0.5,0.,0.5,1.,1.,0.,0.5,1.)
        the=pmas(pycomp(111),1)+pmas(pycomp(-3212),1)
        si4=s1724(ss,ilo4,0,the)*fac
6500    if(isinel(270) == 0)then
        si5=0D0
        goto 6501
        endif
        ik1=321
        ik2=-2212
        call srev(kl,kl1,ik1,ik2,ss,ilo5,fac, &
        0.5,0.5,0.,0.5,1.,1.,0.,0.5,1.)
        the=pmas(pycomp(321),1)+pmas(pycomp(2212),1)
        si5=s1724(ss,ilo5,0,the)*fac
6501    if(isinel(282) == 0)then
        si6=0D0
        goto 650
        endif
        ik1=311
        ik2=-2112
        call srev(kl,kl1,ik1,ik2,ss,ilo6,fac, &
        0.5,0.5,0.,0.5,1.,1.,0.,0.5,1.)
        the=pmas(pycomp(311),1)+pmas(pycomp(2112),1)
        si6=s1724(ss,ilo6,0,the)*fac
650     if(ilo1 == 0 .and. ilo2 == 0 .and. ilo3 == 0 .and. ilo4 == 0 &
         .and. ilo5 == 0 .and. ilo6 == 0) goto 13
        if(si1 < 1D-6.and.si2 < 1D-6.and.si3 < 1D-6.and. &
         si4 < 1D-6.and.si5 < 1D-6.and.si6 < 1D-6) goto 13
        si12=si1+si2
        si13=si12+si3
        si14=si13+si4
        si15=si14+si5
        si16=si15+si6
        s1=si1/si16
        s2=si12/si16
        s3=si13/si16
        s4=si14/si16
        s5=si15/si16
        rlu1=pyr(1)
        if(rlu1 <= s1)then
        ik1=211
        ik2=-3222
        ic=239
!       K- + cascade-ba to pion+ + sigma+ba
        goto 651
        endif
        if(rlu1 > s1 .and. rlu1 <= s2)then
        ik1=-211
        ik2=-3112
        ic=240
!       K- + cascade-ba to pion- + sigma-ba
        goto 651
        endif
        if(rlu1 > s2 .and. rlu1 <= s3)then
        ik1=111
        ik2=-3122
        ic=241
!       K- + cascade-ba to pion0 + lambda0-
        goto 651
        endif
        if(rlu1 > s3 .and. rlu1 <= s4)then
        ik1=111
        ik2=-3212
        ic=242
!       K- + cascade-ba to pion0 + sigma0-
        goto 651
        endif
        if(rlu1 > s4 .and. rlu1 <= s5)then
        ik1=321
        ik2=-2212
        ic=270
!       K- + cascade-ba to k+ + pba
        goto 651
        endif
        ik1=311
        ik2=-2112
        ic=282
!       K- + cascade-ba to k0 +nba
651     lc(icp,3)=ik1
        lc(icp,4)=ik2
        lc(icp,5)=ic
        tw(icp)=si16/cskn/10D0
        goto 10
        endif

!       k0 + cascade-
        if(kl == 311 .and. kl1 == 3312)then
        if(isinel(243) == 0)then
        si1=0D0
        goto 652
        endif
        ik1=-211
        ik2=3122
        call srev(kl,kl1,ik1,ik2,ss,ilo1,fac, &
         0.5,0.5,0.,0.5,1.,0.,0.,0.5,1.)
        the=pmas(pycomp(-211),1)+pmas(pycomp(3122),1)
        si1=s1724(ss,ilo1,0,the)*fac
652     if(isinel(244) == 0)then
        si2=0D0
        goto 653
        endif
        ik1=-211
        ik2=3212
        call srev(kl,kl1,ik1,ik2,ss,ilo2,fac, &
         0.5,0.5,0.,0.5,1.,1.,0.,0.5,1.)
        the=pmas(pycomp(-211),1)+pmas(pycomp(3212),1)
        si2=s1724(ss,ilo2,0,the)*fac
653     if(isinel(245) == 0)then
        si3=0D0
        goto 6540
        endif
        ik1=111
        ik2=3112
        call srev(kl,kl1,ik1,ik2,ss,ilo3,fac, &
        0.5,0.5,0.,0.5,1.,1.,0.,0.5,1.)
        the=pmas(pycomp(111),1)+pmas(pycomp(3112),1)
        si3=s1724(ss,ilo3,0,the)*fac
6540    if(isinel(257) == 0)then
        si4=0D0
        goto 654
        endif
        ik1=-321
        ik2=2112
        call srev(kl,kl1,ik1,ik2,ss,ilo4,fac, &
        0.5,0.5,0.,0.5,1.,1.,0.,0.5,1.)
        the=pmas(pycomp(321),1)+pmas(pycomp(2112),1)
        si4=s1724(ss,ilo4,0,the)*fac
654     if(ilo1 == 0 .and. ilo2 == 0 .and. ilo3 == 0 &
         .and. ilo4 == 0) goto 13
        if(si1 < 1D-6.and.si2 < 1D-6.and.si3 < 1D-6 &
         .and.si4 < 1D-6) goto 13
        si12=si1+si2
        si13=si12+si3
        si14=si13+si4
        s1=si1/si14
        s2=si12/si14
        s3=si13/si14
        rlu1=pyr(1)
        if(rlu1 <= s1)then
        ik1=-211
        ik2=3122
        ic=243
!       k0 + cascade- to pion- + lambda0
        goto 655
        endif
        if(rlu1 > s1 .and. rlu1 <= s2)then
        ik1=-211
        ik2=3212
        ic=244
!       k0 + cascade- to pion- + sigma0
        goto 655
        endif
        if(rlu1 > s2 .and. rlu1 <= s3)then
        ik1=111
        ik2=3112
        ic=245
!       k0 + cascade- to pion0 + sigma-
        goto 655
        endif
        ik1=-321
        ik2=2112
        ic=257
!       k0 + cascade- to k- + n
655     lc(icp,3)=ik1
        lc(icp,4)=ik2
        lc(icp,5)=ic
        tw(icp)=si14/cskn/10D0
        goto 10
        endif

!       k0- + cascade-ba
        if(kl == -311 .and. kl1 == -3312)then
        if(isinel(246) == 0)then
        si1=0D0
        goto 656
        endif
        ik1=211
        ik2=-3122
        call srev(kl,kl1,ik1,ik2,ss,ilo1,fac, &
         0.5,0.5,0.,0.5,1.,0.,0.,0.5,1.)
        the=pmas(pycomp(211),1)+pmas(pycomp(-3122),1)
        si1=s1724(ss,ilo1,0,the)*fac
656     if(isinel(247) == 0)then
        si2=0D0
        goto 657
        endif
        ik1=211
        ik2=-3212
        call srev(kl,kl1,ik1,ik2,ss,ilo2,fac, &
         0.5,0.5,0.,0.5,1.,1.,0.,0.5,1.)
        the=pmas(pycomp(211),1)+pmas(pycomp(-3212),1)
        si2=s1724(ss,ilo2,0,the)*fac
657     if(isinel(248) == 0)then
        si3=0D0
        goto 6580
        endif
        ik1=111
        ik2=-3112
        call srev(kl,kl1,ik1,ik2,ss,ilo3,fac, &
        0.5,0.5,0.,0.5,1.,1.,0.,0.5,1.)
        the=pmas(pycomp(111),1)+pmas(pycomp(-3112),1)
        si3=s1724(ss,ilo3,0,the)*fac
6580    if(isinel(274) == 0)then
        si4=0D0
        goto 658
        endif
        ik1=321
        ik2=-2112
        call srev(kl,kl1,ik1,ik2,ss,ilo4,fac, &
        0.5,0.5,0.,0.5,1.,1.,0.,0.5,1.)
        the=pmas(pycomp(321),1)+pmas(pycomp(2112),1)
        si4=s1724(ss,ilo4,0,the)*fac
658     if(ilo1 == 0 .and. ilo2 == 0 .and. ilo3 == 0 &
         .and. ilo4 == 0) goto 13
        if(si1 < 1D-6.and.si2 < 1D-6.and.si3 < 1D-6 &
         .and.si4 < 1D-6) goto 13
        si12=si1+si2
        si13=si12+si3
        si14=si13+si4
        s1=si1/si14
        s2=si12/si14
        s3=si13/si14
        rlu1=pyr(1)
        if(rlu1 <= s1)then
        ik1=211
        ik2=-3122
        ic=246
!       k0- + cascade-ba to pion+ + lambda0-
        goto 659
        endif
        if(rlu1 > s1 .and. rlu1 <= s2)then
        ik1=211
        ik2=-3212
        ic=247
!       k0- + cascade-ba to pion+ + sigma0-
        goto 659
        endif
        if(rlu1 > s2 .and. rlu1 <= s3)then
        ik1=111
        ik2=-3112
        ic=248
!       k0- + cascade-ba to pion0 + sigma-ba
        goto 659
        endif
        ik1=321
        ik2=-2112
        ic=274
!       k0- + cascade-ba to k+ + nba
659     lc(icp,3)=ik1
        lc(icp,4)=ik2
        lc(icp,5)=ic
        tw(icp)=si14/cskn/10D0
        goto 10
        endif

!       pion+ + lambda0
        if(kl == 211 .and. kl1 == 3122)then
        if(isinel(99) == 0)then
        sigma1=0D0
        goto 2522
        endif
        the=pmas(pycomp(321),1)+pmas(pycomp(3322),1)
        sigma1=s1724(ss,ilo1,0,the)
!       pion+ + lambda0 to k+ + cascade
2522    if(isinel(252) == 0)then
        sigma2=0D0
        goto 2523
        endif
        ik3=-311
         ik4=2212
        call srev(kl,kl1,ik3,ik4,ss,ilo2,fac, &
         1.,0.,0.,0.5,0.5,0.5,0.,0.5,1.)
        the=pmas(pycomp(-311),1)+pmas(pycomp(2212),1)
        sigma2=fac*s1724(ss,ilo2,0,the)*10D0

!       cross section of     pion+ + lambda0 to k0- + p
2523    if(ilo1 == 0.and.ilo2 == 0) goto 13
        if(sigma1 < 1D-6 .and. sigma2 < 1D-6) goto 13
        ik1=321
        ik2=3322
        ic=99
!       pion+ + lambda0 to k+ + cascade0
        sigm12=sigma1+sigma2
        if(pyr(1) > sigma1/sigm12)then
        ik1=-311
        ik2=2212
        ic=252
!       cross section of     pion+ + lambda0 to k0- + p
        endif
        lc(icp,3)=ik1
        lc(icp,4)=ik2
        lc(icp,5)=ic
        tw(icp)=sigm12/cspin/10D0
        goto 10
        endif

!       pion- + lambda0-
        if(kl == -211 .and. kl1 == -3122)then
        if(isinel(107) == 0)then
        sigma1=0D0
        goto 3522
        endif
        the=pmas(pycomp(-321),1)+pmas(pycomp(-3322),1)
        sigma1=s1724(ss,ilo1,0,the)
!       pion- + lambda- to k- + cascade0-
3522    if(isinel(275) == 0)then
        sigma2=0D0
        goto 3523
        endif
        ik3=311
        ik4=-2212
        call srev(kl,kl1,ik3,ik4,ss,ilo2,fac, &
         1.,0.,0.,0.5,0.5,0.5,0.,0.5,1.)
        the=pmas(pycomp(311),1)+pmas(pycomp(-2212),1)
        sigma2=fac*s1724(ss,ilo2,0,the)*10D0
!       pion- + lambda0- to k0 + pba

3523    if(ilo1 == 0.and.ilo2 == 0) goto 13
        if(sigma1 < 1D-6 .and. sigma2 < 1D-6) goto 13
        ik1=-321
        ik2=-3322
        ic=107
!       pion- + lambda- to k- + cascade0-
        sigm12=sigma1+sigma2
        if(pyr(1) > sigma1/sigm12)then
        ik1=311
        ik2=-2212
        ic=275
!       pion- + lambda0- to k0 + pba
        endif
        lc(icp,3)=ik1
        lc(icp,4)=ik2
        lc(icp,5)=ic
        tw(icp)=sigm12/cspin/10D0
        goto 10
        endif

!       pion- + lambda0
        if(kl == -211 .and. kl1 == 3122)then
        if(isinel(36) == 0)then
        sigma1=0D0
        goto 1522
        endif
        the=pmas(pycomp(311),1)+pmas(pycomp(3312),1)
        sigma1=s1724(ss,ilo1,0,the)
!       pion- + lambda to k0 + cascade-
1522    if(isinel(258) == 0)then
        sigma2=0D0
        goto 1523
        endif
        ik3=-321
        ik4=2112
        call srev(kl,kl1,ik3,ik4,ss,ilo2,fac, &
         1.,0.,0.,0.5,0.5,0.5,0.,0.5,1.)
        the=pmas(pycomp(-321),1)+pmas(pycomp(2112),1)
        sigma2=fac*s1724(ss,ilo2,0,the)*10D0

!       cross section of   pion- + lambda0 to k- + n
1523    if(ilo1 == 0.and.ilo2 == 0) goto 13
        if(sigma1 < 1D-6 .and. sigma2 < 1D-6)  goto 13
        ik1=311
        ik2=3312
        ic=36
!       pion- + lambda to k0 + cascade-
        sigm12=sigma1+sigma2
        if(pyr(1) > sigma1/sigm12)then
        ik1=-321
        ik2=2112
        ic=258
!       cross section of pion- + lambda0 to k- + n
        endif
        lc(icp,3)=ik1
        lc(icp,4)=ik2
        lc(icp,5)=ic
        tw(icp)=sigm12/cspin/10D0
        goto 10
        endif
!       pion+ + lambda-
        if(kl == 211 .and. kl1 == -3122)then
        if(isinel(42) == 0)then
        sigma1=0D0
        goto 6522
        endif
        the=pmas(pycomp(-311),1)+pmas(pycomp(-3312),1)
        sigma1=s1724(ss,ilo1,0,the)
!       pion+ + lambdaba to k0- + cascade-ba
6522    if(isinel(269) == 0)then
        sigma2=0D0
        goto 6523
        endif
        ik3=321
        ik4=-2112
        call srev(kl,kl1,ik3,ik4,ss,ilo2,fac, &
         1.,0.,0.,0.5,0.5,0.5,0.,0.5,1.)
        the=pmas(pycomp(321),1)+pmas(pycomp(-2112),1)
        sigma2=fac*s1724(ss,ilo2,0,the)*10D0

!       pion+ + lambda0- to k+ + nba
6523    if(ilo1 == 0.and.ilo2 == 0) goto 13
        if(sigma1 < 1D-6 .and. sigma2 < 1D-6) goto 13
        ik1=-311
        ik2=-3312
        ic=42
!       pion+ + lambdaba to k0- + cascade-ba
        sigm12=sigma1+sigma2
        if(pyr(1) > sigma1/sigm12)then
        ik1=321
        ik2=-2112
        ic=269
!       pion+ + lambda0- to k+ + nba
        endif
        lc(icp,3)=ik1
        lc(icp,4)=ik2
        lc(icp,5)=ic
        tw(icp)=sigm12/cspin/10D0
        goto 10
        endif

!       pion0 + lambda0
        if(kl == 111 .and. kl1 == 3122)then
        if(isinel(263) == 0)then
        si1=0D0
        goto 669
        endif
        ik1=-321
        ik2=2212
        call srev(kl,kl1,ik1,ik2,ss,ilo1,fac, &
         1.,0.,0.,0.5,0.5,0.5,0.,0.5,1.)
        the=pmas(pycomp(-321),1)+pmas(pycomp(2212),1)
        si1=10D0*s1724(ss,ilo1,0,the)*fac
669     if(isinel(264) == 0)then
        si2=0D0
        goto 670
        endif
        ik1=-311
        ik2=2112
        call srev(kl,kl1,ik1,ik2,ss,ilo2,fac, &
         1.,0.,0.,0.5,0.5,0.5,0.,0.5,1.)
        the=pmas(pycomp(-311),1)+pmas(pycomp(2112),1)
        si2=10D0*s1724(ss,ilo2,0,the)*fac
670      if(isinel(39) == 0)then
        si3=0D0
        goto 1670
        endif
        the=pmas(pycomp(321),1)+pmas(pycomp(3312),1)
        si3=s1724(ss,ilo3,0,the)
1670      if(isinel(103) == 0)then
        si4=0D0
        goto 1671
        endif
        the=pmas(pycomp(311),1)+pmas(pycomp(3322),1)
        si4=s1724(ss,ilo4,0,the)
1671    if(ilo1 == 0.and.ilo2 == 0.and.ilo3 == 0.and.ilo4 == 0) goto 13
        if(si1 < 1D-6.and.si2 < 1D-6.and.si3 < 1D-6 &
         .and.si4 < 1D-6) goto 13
        si12=si1+si2
        si13=si12+si3
        si14=si13+si4
        s1=si1/si14
        s2=si12/si14
        s3=si13/si14
        rlus=pyr(1)
        if(rlus <= s1)then
        ik1=-321
        ik2=2212
        ic=263
!       pion0 + lambda0 to k- + p
        goto 671
        endif
        if(rlus > s1 .and. rlus <= s2)then
        ik1=-311
        ik2=2112
        ic=264
!       pion0 + lambda0 to k0- + n
        goto 671
        endif
        if(rlus > s2 .and. rlus <= s3)then
        ik1=321
        ik2=3312
        ic=39
!       pion0 + lambda0 to  k+ + cascade-
        goto 671
        endif
        ik1=311
        ik2=3322
        ic=103
!       pion0 + lambda to k0 + cascade0
671     lc(icp,3)=ik1
        lc(icp,4)=ik2
        lc(icp,5)=ic
        tw(icp)=si14/cspin/10D0
        goto 10
        endif
!       pion0 + lambda0-
        if(kl == 111 .and. kl1 == -3122)then
        if(isinel(280) == 0)then
        si1=0.
        goto 2669
        endif
        ik1=321
        ik2=-2212
        call srev(kl,kl1,ik1,ik2,ss,ilo1,fac, &
         1.,0.,0.,0.5,0.5,0.5,0.,0.5,1.)
        the=pmas(pycomp(321),1)+pmas(pycomp(-2212),1)
        si1=10D0*s1724(ss,ilo1,0,the)*fac
2669    if(isinel(281) == 0)then
        si2=0D0
        goto 2677
        endif
        ik1=311
        ik2=-2112
        call srev(kl,kl1,ik1,ik2,ss,ilo2,fac, &
         1.,0.,0.,0.5,0.5,0.5,0.,0.5,1.)
        the=pmas(pycomp(311),1)+pmas(pycomp(-2112),1)
        si2=10D0*s1724(ss,ilo2,0,the)*fac
2677      if(isinel(46) == 0)then
        si3=0D0
        goto 2670
        endif
        the=pmas(pycomp(-321),1)+pmas(pycomp(-3312),1)
        si3=s1724(ss,ilo3,0,the)
!       pion0 + lambda- to k- + cascade-ba
2670      if(isinel(110) == 0)then
        si4=0.
        goto 2671
        endif
        the=pmas(pycomp(-311),1)+pmas(pycomp(-3322),1)
        si4=s1724(ss,ilo4,0,the)
2671    if(ilo1 == 0.and.ilo2 == 0.and.ilo3 == 0.and.ilo4 == 0) goto 13
        if(si1 < 1D-6.and.si2 < 1D-6.and.si3 < 1D-6 &
         .and.si4 < 1D-6) goto 13
        si12=si1+si2
        si13=si12+si3
        si14=si13+si4
        s1=si1/si14
        s2=si12/si14
        s3=si13/si14
        rlus=pyr(1)
        if(rlus <= s1)then
        ik1=321
        ik2=-2212
        ic=280
!       pion0 + lambda0- to k+ + pbar
        goto 6711
        endif
        if(rlus > s1 .and. rlus <= s2)then
        ik1=311
        ik2=-2112
        ic=281
!       pion0 + lambda0- to k0 + nba
        goto 6711
        endif
        if(rlus > s2 .and. rlus <= s3)then
        ik1=-321
        ik2=-3312
        ic=46
!       pion0 + lambda- to k- + cascade-ba
        goto 6711
        endif
        ik1=-311
        ik2=-3322
        ic=110
!          pion0 + lambda- to k0- + cascade0-
6711    lc(icp,3)=ik1
        lc(icp,4)=ik2
        lc(icp,5)=ic
        tw(icp)=si14/cspin/10D0
        goto 10
        endif

!       pion+ + cascade-
        if(kl == 211 .and. kl1 == 3312)then
        if(isinel(283) == 0)then
        si1=0D0
        goto 684
        endif
        ik1=-321
        ik2=3222
        call srev(kl,kl1,ik1,ik2,ss,ilo1,fac, &
         1.,0.5,0.,0.5,0.5,1.,0.,0.5,1.)
        the=pmas(pycomp(-321),1)+pmas(pycomp(3222),1)
        si1=10D0*s1724(ss,ilo1,0,the)*fac
684     if(isinel(284) == 0)then
        si2=0D0
        goto 685
        endif
        ik1=-311
        ik2=3122
        call srev(kl,kl1,ik1,ik2,ss,ilo2,fac, &
         1.,0.5,0.,0.5,0.5,0.,0.,0.5,1.)
        the=pmas(pycomp(-311),1)+pmas(pycomp(3122),1)
        si2=10D0*s1724(ss,ilo2,0,the)*fac
685     if(isinel(285) == 0)then
        si3=0D0
        goto 686
        endif
        ik1=-311
        ik2=3212
        call srev(kl,kl1,ik1,ik2,ss,ilo3,fac, &
         1.,0.5,0.,0.5,0.5,1.,0.,0.5,1.)
        the=pmas(pycomp(-311),1)+pmas(pycomp(3212),1)
        si3=10D0*s1724(ss,ilo3,0,the)*fac
686     if(isinel(131) == 0)then
        si4=0D0
        goto 1686
        endif
        the=pmas(pycomp(321),1)+pmas(pycomp(3334),1)
        si4=s1724(ss,ilo4,0,the)

1686    if(ilo1 == 0 .and. ilo2 == 0 .and. ilo3 == 0 &
         .and.ilo4 == 0) goto 13
        if(si1 < 1D-6.and.si2 < 1D-6.and.si3 < 1D-6.and. &
         si4 < 1D-6) goto 13
        si12=si1+si2
        si13=si12+si3
        si14=si13+si4
        s1=si1/si14
        s2=si12/si14
        s3=si13/si14
        rlu1=pyr(1)
        if(rlu1 <= s1)then
        ik1=-321
        ik2=3222
        ic=283
!       pion+ + cascade- to k- + sigma+
        goto 687
        endif
        if(rlu1 > s1 .and. rlu1 <= s2)then
        ik1=-311
        ik2=3122
        ic=284
!       pion+ + cascade- to k0- + lambda0
        goto 687
        endif
        if(rlu1 > s2 .and. rlu1 <= s3)then
        ik1=-311
        ik2=3212
        ic=285
!       pion+ + cascade- to k0- + sigma0
        goto 687
        endif
        ik1=321
        ik2=3334
        ic=131
!       pion+ + cascade- to k+ + Omega-
687     lc(icp,3)=ik1
        lc(icp,4)=ik2
        lc(icp,5)=ic
        tw(icp)=si14/cspin/10D0
        goto 10
        endif

!       pion- + cascade- to k- + sigma-
        if(kl == -211 .and. kl1 == 3312)then
        if(isinel(286) == 0) goto 13
        the=pmas(pycomp(-321),1)+pmas(pycomp(3112),1)
        ww=s1724(ss,ilo,0,the)/cspin
        if(ilo == 0) goto 13
        lc(icp,3)=-321
        lc(icp,4)=3112
        lc(icp,5)=286
        ik3=-321
        ik4=3112
        call srev(kl,kl1,ik3,ik4,ss,ilo1,fac, &
         1.,0.5,0.,0.5,0.5,1.,0.,0.5,1.)
        tw(icp)=ww*fac
        goto 10
        endif

!       pion0 + cascade-
        if(kl == 111 .and. kl1 == 3312)then
        if(isinel(287) == 0)then
        si1=0D0
        goto 688
        endif
        ik1=-321
        ik2=3122
        call srev(kl,kl1,ik1,ik2,ss,ilo1,fac, &
         1.,0.5,0.,0.5,0.5,0.,0.,0.5,1.)
        the=pmas(pycomp(-321),1)+pmas(pycomp(3122),1)
        si1=10D0*s1724(ss,ilo1,0,the)*fac
688     if(isinel(288) == 0)then
        si2=0D0
        goto 689
        endif
        ik1=-321
        ik2=3212
        call srev(kl,kl1,ik1,ik2,ss,ilo2,fac, &
         1.,0.5,0.,0.5,0.5,1.,0.,0.5,1.)
        the=pmas(pycomp(-321),1)+pmas(pycomp(3212),1)
        si2=10D0*s1724(ss,ilo2,0,the)*fac
689     if(isinel(289) == 0)then
        si3=0D0
        goto 690
        endif
        ik1=-311
        ik2=3112
        call srev(kl,kl1,ik1,ik2,ss,ilo3,fac, &
         1.,0.5,0.,0.5,0.5,1.,0.,0.5,1.)
        the=pmas(pycomp(-311),1)+pmas(pycomp(3112),1)
        si3=10D0*s1724(ss,ilo3,0,the)*fac
690     if(isinel(132) == 0)then
        si4=0D0
        goto 1690
        endif
        the=pmas(pycomp(311),1)+pmas(pycomp(3334),1)
        si4=s1724(ss,ilo4,0,the)

!       pion0 + cascade- to k0 + Omega-
1690    if(ilo1 == 0 .and. ilo2 == 0 .and. ilo3 == 0.and. &
         ilo4 == 0) goto 13
        if(si1 < 1D-6.and.si2 < 1D-6.and.si3 < 1D-6.and. &
         si4 < 1D-6) goto 13
        si12=si1+si2
        si13=si12+si3
        si14=si13+si4
        s1=si1/si14
        s2=si12/si14
        s3=si13/si14
        rlu1=pyr(1)
        if(rlu1 <= s1)then
        ik1=-321
        ik2=3122
        ic=287
!       pion0 + cascade- to k- + lambda0
        goto 691
        endif
        if(rlu1 > s1 .and. rlu1 <= s2)then
        ik1=-321
        ik2=3212
        ic=288
!       pion0 + cascade- to k- + sigma0
        goto 691
        endif
        if(rlu1 > s2 .and. rlu1 <= s3)then
        ik1=-311
        ik2=3112
        ic=289
!       pion0 + cascade- to k0- + sigma-
        goto 691
        endif
        ik1=311
        ik2=3334
        ic=132
!       pion0 + cascade- to k0 + Omega-
691     lc(icp,3)=ik1
        lc(icp,4)=ik2
        lc(icp,5)=ic
        tw(icp)=si14/cspin/10D0
        goto 10
        endif

!       pion+ + cascade-ba to k+ + sigma-ba
        if(kl == 211 .and. kl1 == -3312)then
        if(isinel(290) == 0) goto 13
        the=pmas(pycomp(321),1)+pmas(pycomp(-3112),1)
        ww=s1724(ss,ilo,0,the)/cspin
        if(ilo == 0) goto 13
        lc(icp,3)=321
        lc(icp,4)=-3112
        lc(icp,5)=290
        ik3=321
        ik4=-3112
        call srev(kl,kl1,ik3,ik4,ss,ilo1,fac, &
         1.,0.5,0.,0.5,0.5,1.,0.,0.5,1.)

        tw(icp)=ww*fac
        goto 10
        endif

!       pion- + cascade-ba
        if(kl == -211 .and. kl1 == -3312)then
        if(isinel(291) == 0)then
        si1=0D0
        goto 692
        endif
        ik1=321
        ik2=-3222
        call srev(kl,kl1,ik1,ik2,ss,ilo1,fac, &
         1.,0.5,0.,0.5,0.5,1.,0.,0.5,1.)
        the=pmas(pycomp(321),1)+pmas(pycomp(-3222),1)
        si1=10D0*s1724(ss,ilo1,0,the)*fac
692     if(isinel(292) == 0)then
        si2=0D0
        goto 693
        endif
        ik1=311
        ik2=-3122
        call srev(kl,kl1,ik1,ik2,ss,ilo2,fac, &
         1.,0.5,0.,0.5,0.5,0.,0.,0.5,1.)
        the=pmas(pycomp(311),1)+pmas(pycomp(-3122),1)
        si2=10D0*s1724(ss,ilo2,0,the)*fac
693     if(isinel(293) == 0)then
        si3=0D0
        goto 694
        endif
        ik1=311
        ik2=-3212
        call srev(kl,kl1,ik1,ik2,ss,ilo3,fac, &
         1.,0.5,0.,0.5,0.5,1.,0.,0.5,1.)
        the=pmas(pycomp(311),1)+pmas(pycomp(-3212),1)
        si3=10D0*s1724(ss,ilo3,0,the)*fac
694     if(isinel(133) == 0)then
        si4=0D0
        goto 1694
        endif
        the=pmas(pycomp(-321),1)+pmas(pycomp(-3334),1)
        si4=s1724(ss,ilo4,0,the)
1694    if(ilo1 == 0 .and. ilo2 == 0 .and. ilo3 == 0.and. &
         ilo4 == 0) goto 13
        if(si1 < 1D-6.and.si2 < 1D-6.and.si3 < 1D-6.and. &
         si4 < 1D-6) goto 13
        si12=si1+si2
        si13=si12+si3
        si14=si13+si4
        s1=si1/si14
        s2=si12/si14
        s3=si13/si14
        rlu1=pyr(1)
        if(rlu1 <= s1)then
        ik1=321
        ik2=-3222
        ic=291
!       pion- + cascade-ba to k+ + sigma+-
        goto 695
        endif
        if(rlu1 > s1 .and. rlu1 <= s2)then
        ik1=311
        ik2=-3122
        ic=292
!       pion- + cascade-ba to k0 + lambda0-
        goto 695
        endif
        if(rlu1 > s2 .and. rlu1 <= s3)then
        ik1=311
        ik2=-3212
        ic=293
!       pion- + cascade-ba to k0 + sigma0-
        goto 695
        endif
        ik1=-321
        ik2=-3334
        ic=133
!       pion- + cascade-ba to k- + Omega-ba
695     lc(icp,3)=ik1
        lc(icp,4)=ik2
        lc(icp,5)=ic
        tw(icp)=si14/cspin/10D0
        goto 10
        endif

!       pion0 + cascade-ba
        if(kl == 111 .and. kl1 == -3312)then
        if(isinel(294) == 0)then
        si1=0.
        goto 696
        endif
        ik1=321
        ik2=-3122
        call srev(kl,kl1,ik1,ik2,ss,ilo1,fac, &
         1.,0.5,0.,0.5,0.5,0.,0.,0.5,1.)
        the=pmas(pycomp(321),1)+pmas(pycomp(-3122),1)
        si1=10D0*s1724(ss,ilo1,0,the)*fac
696     if(isinel(295) == 0)then
        si2=0D0
        goto 697
        endif
        ik1=321
        ik2=-3212
        call srev(kl,kl1,ik1,ik2,ss,ilo2,fac, &
         1.,0.5,0.,0.5,0.5,1.,0.,0.5,1.)
        the=pmas(pycomp(321),1)+pmas(pycomp(-3212),1)
        si2=10D0*s1724(ss,ilo2,0,the)*fac
697     if(isinel(296) == 0)then
        si3=0D0
        goto 698
        endif
        ik1=311
        ik2=-3112
        call srev(kl,kl1,ik1,ik2,ss,ilo3,fac, &
         1.,0.5,0.,0.5,0.5,1.,0.,0.5,1.)
        the=pmas(pycomp(311),1)+pmas(pycomp(-3112),1)
        si3=10D0*s1724(ss,ilo3,0,the)*fac
698     if(isinel(134) == 0)then
        si4=0D0
        goto 1698
        endif
         the=pmas(pycomp(-311),1)+pmas(pycomp(-3334),1)
        si4=s1724(ss,ilo4,0,the)
1698    if(ilo1 == 0 .and. ilo2 == 0 .and. ilo3 == 0.and. &
         ilo4 == 0) goto 13
        if(si1 < 1D-6.and.si2 < 1D-6.and.si3 < 1D-6.and. &
         si4 < 1D-6) goto 13
        si12=si1+si2
        si13=si12+si3
        si14=si13+si4
        s1=si1/si14
        s2=si12/si14
        s3=si13/si14
        rlu1=pyr(1)
        if(rlu1 <= s1)then
        ik1=321
        ik2=-3122
        ic=294
!       pion0 + cascade-ba to k+ + lambda0-
        goto 700
        endif
        if(rlu1 > s1 .and. rlu1 <= s2)then
        ik1=321
        ik2=-3212
        ic=295
!       pion0 + cascade-ba to k+ + sigma0-
        goto 700
        endif
        if(rlu1 > s2 .and. rlu1 <= s3)then
        ik1=311
        ik2=-3112
        ic=296
!       pion0 + cascade-ba to k0 + sigma-ba
        goto 700
        endif
        ik1=-311
        ik2=-3334
        ic=134
!       pion0 + cascade-ba to k0- + Omega-ba
700     lc(icp,3)=ik1
        lc(icp,4)=ik2
        lc(icp,5)=ic
        tw(icp)=si14/cspin/10D0
        goto 10
        endif

!       K+ + cascade0
        if(kl == 321 .and. kl1 == 3322)then
        if(isinel(298) == 0)then
        si1=0D0
        goto 768
        endif
        ik1=211
        ik2=3212
        call srev(kl,kl1,ik1,ik2,ss,ilo1,fac, &
         0.5,0.5,0.,0.5,1.,1.,0.,0.5,1.)
        the=pmas(pycomp(211),1)+pmas(pycomp(3212),1)
        si1=s1724(ss,ilo1,0,the)*fac
768     if(isinel(299) == 0)then
        si2=0D0
        goto 769
        endif
        ik1=211
        ik2=3122
        call srev(kl,kl1,ik1,ik2,ss,ilo2,fac, &
         0.5,0.5,0.,0.5,1.,0.,0.,0.5,1.)
        the=pmas(pycomp(211),1)+pmas(pycomp(3122),1)
        si2=s1724(ss,ilo2,0,the)*fac
769     if(isinel(301) == 0)then
        si3=0D0
        goto 770
        endif
        ik1=111
        ik2=3222
        call srev(kl,kl1,ik1,ik2,ss,ilo3,fac, &
         0.5,0.5,0.,0.5,1.,1.,0.,0.5,1.)
        the=pmas(pycomp(111),1)+pmas(pycomp(3222),1)
        si3=s1724(ss,ilo3,0,the)*fac
770     if(isinel(326) == 0)then
        si4=0D0
        goto 771
        endif
        ik1=-311
        ik2=2212
        call srev(kl,kl1,ik1,ik2,ss,ilo4,fac, &
         0.5,0.5,0.,0.5,0.5,0.5,0.,0.5,1.)
        the=pmas(pycomp(-311),1)+pmas(pycomp(2212),1)
        si4=s1724(ss,ilo4,0,the)*fac
771     if(ilo1 == 0 .and. ilo2 == 0 .and. ilo3 == 0 .and. ilo4 == 0 &
         ) goto 13
        if(si1 < 1D-6.and.si2 < 1D-6.and.si3 < 1D-6 &
         .and.si4 < 1D-6) goto 13
        si12=si1+si2
        si13=si12+si3
        si14=si13+si4
        s1=si1/si14
        s2=si12/si14
        s3=si13/si14
        rlu1=pyr(1)
        if(rlu1 <= s1)then
        ik1=211
        ik2=3212
        ic=298
!       K+ + cascade0 to pion+ + sigma0
        goto 772
        endif
        if(rlu1 > s1 .and. rlu1 <= s2)then
        ik1=211
        ik2=3122
        ic=299
!       K+ + cascade0 to pion+ + lambda
        goto 772
        endif
        if(rlu1 > s2 .and. rlu1 <= s3)then
        ik1=111
        ik2=3222
        ic=301
!       K+ + cascade0 to pion0 + sigma+
        goto 772
        endif
        ik1=-311
        ik2=2212
        ic=326
!       K+ + cascade0 to k0- + p
772     lc(icp,3)=ik1
        lc(icp,4)=ik2
        lc(icp,5)=ic
        tw(icp)=si14/cskn/10D0
        goto 10
        endif

!       K- + cascad0-
        if(kl == -321 .and. kl1 == -3322)then
        if(isinel(306) == 0)then
        si1=0D0
        goto 773
        endif
        ik1=-211
        ik2=-3212
        call srev(kl,kl1,ik1,ik2,ss,ilo1,fac, &
         0.5,0.5,0.,0.5,1.,1.,0.,0.5,1.)
        the=pmas(pycomp(-211),1)+pmas(pycomp(-3212),1)
        si1=s1724(ss,ilo1,0,the)*fac
773     if(isinel(307) == 0)then
        si2=0D0
        goto 774
        endif
        ik1=-211
        ik2=-3122
        call srev(kl,kl1,ik1,ik2,ss,ilo2,fac, &
         0.5,0.5,0.,0.5,1.,0.,0.,0.5,1.)
        the=pmas(pycomp(-211),1)+pmas(pycomp(3122),1)
        si2=s1724(ss,ilo2,0,the)*fac
774     if(isinel(308) == 0)then
        si3=0D0
        goto 775
        endif
        ik1=111
        ik2=-3222
        call srev(kl,kl1,ik1,ik2,ss,ilo3,fac, &
         0.5,0.5,0.,0.5,1.,1.,0.,0.5,1.)
        the=pmas(pycomp(111),1)+pmas(pycomp(-3222),1)
        si3=s1724(ss,ilo3,0,the)*fac
775     if(isinel(329) == 0)then
        si4=0D0
        goto 776
        endif
        ik1=311
        ik2=-2212
        call srev(kl,kl1,ik1,ik2,ss,ilo4,fac, &
         0.5,0.5,0.,0.5,0.5,0.5,0.,0.5,1.)
        the=pmas(pycomp(311),1)+pmas(pycomp(-2212),1)
        si4=s1724(ss,ilo4,0,the)*fac
776     if(ilo1 == 0 .and. ilo2 == 0 .and. ilo3 == 0 .and. ilo4 == 0 &
         ) goto 13
        if(si1 < 1D-6.and.si2 < 1D-6.and.si3 < 1D-6 &
         .and.si4 < 1D-6) goto 13
        si12=si1+si2
        si13=si12+si3
        si14=si13+si4
        s1=si1/si14
        s2=si12/si14
        s3=si13/si14
        rlu1=pyr(1)
        if(rlu1 <= s1)then
        ik1=-211
        ik2=-3212
        ic=306
!       K- + cascade0- to pion- + sigma0-
        goto 777
        endif
        if(rlu1 > s1 .and. rlu1 <= s2)then
        ik1=-211
        ik2=-3122
        ic=307
!       K- + cascade0- to pion- + lambda-
        goto 777
        endif
        if(rlu1 > s2 .and. rlu1 <= s3)then
        ik1=111
        ik2=-3222
        ic=308
!       K- + cascade0- to pion0 + sigma+ba
        goto 777
        endif
        ik1=311
        ik2=-2212
        ic=329
!       K- + cascade0- to k0 + p-
777     lc(icp,3)=ik1
        lc(icp,4)=ik2
        lc(icp,5)=ic
        tw(icp)=si14/cskn/10D0
        goto 10
        endif

!       k0 + cascad0
        if(kl == 311 .and. kl1 == 3322)then
        if(isinel(297) == 0)then
        si1=0D0
        goto 778
        endif
        ik1=211
        ik2=3112
        call srev(kl,kl1,ik1,ik2,ss,ilo1,fac, &
         0.5,0.5,0.,0.5,1.,1.,0.,0.5,1.)
        the=pmas(pycomp(211),1)+pmas(pycomp(3112),1)
        si1=s1724(ss,ilo1,0,the)*fac
778     if(isinel(300) == 0)then
        si2=0D0
        goto 779
        endif
        ik1=-211
        ik2=3222
        call srev(kl,kl1,ik1,ik2,ss,ilo2,fac, &
         0.5,0.5,0.,0.5,1.,1.,0.,0.5,1.)
        the=pmas(pycomp(-211),1)+pmas(pycomp(3222),1)
        si2=s1724(ss,ilo2,0,the)*fac
779     if(isinel(302) == 0)then
        si3=0D0
        goto 780
        endif
        ik1=111
        ik2=3212
        call srev(kl,kl1,ik1,ik2,ss,ilo3,fac, &
         0.5,0.5,0.,0.5,1.,1.,0.,0.5,1.)
        the=pmas(pycomp(111),1)+pmas(pycomp(3212),1)
        si3=s1724(ss,ilo3,0,the)*fac
780     if(isinel(303) == 0)then
        si4=0D0
        goto 781
        endif
        ik1=111
        ik2=3122
        call srev(kl,kl1,ik1,ik2,ss,ilo4,fac, &
         0.5,0.5,0.,0.5,1.,0.,0.,0.5,1.)
        the=pmas(pycomp(111),1)+pmas(pycomp(3122),1)
        si4=s1724(ss,ilo4,0,the)*fac
781     if(isinel(325) == 0)then
        si5=0D0
        goto 782
        endif
        ik1=-321
        ik2=2212
        call srev(kl,kl1,ik1,ik2,ss,ilo5,fac, &
         0.5,0.5,0.,0.5,0.5,0.5,0.,0.5,1.)
        the=pmas(pycomp(-321),1)+pmas(pycomp(2212),1)
        si5=s1724(ss,ilo5,0,the)*fac
782     if(isinel(327) == 0)then
        si6=0D0
        goto 783
        endif
        ik1=-311
        ik2=2112
        call srev(kl,kl1,ik1,ik2,ss,ilo6,fac, &
         0.5,0.5,0.,0.5,0.5,0.5,0.,0.5,1.)
        the=pmas(pycomp(-311),1)+pmas(pycomp(2112),1)
        si6=s1724(ss,ilo6,0,the)*fac
783     if(ilo1 == 0 .and. ilo2 == 0 .and. ilo3 == 0 .and. ilo4 == 0 &
         .and. ilo5 == 0.and. ilo6 == 0) goto 13
        if(si1 < 1D-6.and.si2 < 1D-6.and.si3 < 1D-6 &
         .and.si4 < 1D-6.and.si5 < 1D-6.and.si6 < 1D-6) goto 13
        si12=si1+si2
        si13=si12+si3
        si14=si13+si4
        si15=si14+si5
        si16=si15+si6
        s1=si1/si16
        s2=si12/si16
        s3=si13/si16
        s4=si14/si16
        s5=si15/si16
        rlu1=pyr(1)
        if(rlu1 <= s1)then
        ik1=211
        ik2=3112
        ic=297
!       k0 + cascade0 to pion+ + sigma-
        goto 784
        endif
        if(rlu1 > s1 .and. rlu1 <= s2)then
        ik1=-211
        ik2=3222
        ic=300
!       k0 + cascade0 to pion- + sigma+
        goto 784
        endif
        if(rlu1 > s2 .and. rlu1 <= s3)then
        ik1=111
        ik2=3212
        ic=302
!       k0 + cascade0 to pion0 + sigma0
        goto 784
        endif
        if(rlu1 > s3 .and. rlu1 <= s4)then
        ik1=111
        ik2=3122
        ic=303
!       k0 + cascade0 to pion0 + lambda
        goto 784
        endif
        if(rlu1 > s4 .and. rlu1 <= s5)then
        ik1=-321
        ik2=2212
        ic=325
!       k0 + cascade0 to k- + p
        goto 784
        endif
        ik1=-311
        ik2=2112
        ic=327
!       k0 + cascade0 to k0- + n
784     lc(icp,3)=ik1
        lc(icp,4)=ik2
        lc(icp,5)=ic
        tw(icp)=si16/cskn/10D0
        goto 10
        endif

!       k0- + cascad0-
        if(kl == -311 .and. kl1 == -3322)then
        if(isinel(304) == 0)then
        si1=0D0
        goto 785
        endif
        ik1=211
        ik2=-3222
        call srev(kl,kl1,ik1,ik2,ss,ilo1,fac, &
         0.5,0.5,0.,0.5,1.0,1.0,0.,0.5,1.)
        the=pmas(pycomp(211),1)+pmas(pycomp(-3222),1)
        si1=s1724(ss,ilo1,0,the)*fac
785     if(isinel(305) == 0)then
        si2=0D0
        goto 786
        endif
        ik1=-211
        ik2=-3112
        call srev(kl,kl1,ik1,ik2,ss,ilo2,fac, &
         0.5,0.5,0.,0.5,1.0,1.0,0.,0.5,1.)
        the=pmas(pycomp(-211),1)+pmas(pycomp(-3112),1)
        si2=s1724(ss,ilo2,0,the)*fac
786     if(isinel(309) == 0)then
        si3=0D0
        goto 787
        endif
        ik1=111
        ik2=-3212
        call srev(kl,kl1,ik1,ik2,ss,ilo3,fac, &
         0.5,0.5,0.,0.5,1.0,1.0,0.,0.5,1.)
        the=pmas(pycomp(111),1)+pmas(pycomp(-3212),1)
        si3=s1724(ss,ilo3,0,the)*fac
787     if(isinel(310) == 0)then
        si4=0D0
        goto 788
        endif
        ik1=111
        ik2=-3122
        call srev(kl,kl1,ik1,ik2,ss,ilo4,fac, &
         0.5,0.5,0.,0.5,1.0,0.,0.,0.5,1.)
        the=pmas(pycomp(111),1)+pmas(pycomp(-3122),1)
        si4=s1724(ss,ilo4,0,the)*fac
788     if(isinel(328) == 0)then
        si5=0D0
        goto 789
        endif
        ik1=321
        ik2=-2212
        call srev(kl,kl1,ik1,ik2,ss,ilo5,fac, &
         0.5,0.5,0.,0.5,0.5,0.5,0.,0.5,1.)
        the=pmas(pycomp(321),1)+pmas(pycomp(-2212),1)
        si5=s1724(ss,ilo5,0,the)*fac
789     if(isinel(330) == 0)then
        si6=0D0
        goto 790
        endif
        ik1=311
        ik2=-2112
        call srev(kl,kl1,ik1,ik2,ss,ilo6,fac, &
         0.5,0.5,0.,0.5,0.5,0.5,0.,0.5,1.)
        the=pmas(pycomp(311),1)+pmas(pycomp(-2112),1)
        si6=s1724(ss,ilo6,0,the)*fac
790     if(ilo1 == 0 .and. ilo2 == 0 .and. ilo3 == 0 .and. ilo4 == 0 &
         .and. ilo5 == 0.and. ilo6 == 0) goto 13
        if(si1 < 1D-6.and.si2 < 1D-6.and.si3 < 1D-6 &
         .and.si4 < 1D-6.and.si5 < 1D-6.and.si6 < 1D-6) goto 13
        si12=si1+si2
        si13=si12+si3
        si14=si13+si4
        si15=si14+si5
        si16=si15+si6
        s1=si1/si16
        s2=si12/si16
        s3=si13/si16
        s4=si14/si16
        s5=si15/si16
        rlu1=pyr(1)
        if(rlu1 <= s1)then
        ik1=211
        ik2=-3222
        ic=304
!       K-0 + cascade0- to pion+ + sigma-ba
        goto 791
        endif
        if(rlu1 > s1 .and. rlu1 <= s2)then
        ik1=-211
        ik2=-3112
        ic=305
!       k0- + cascade0- to pion- + sigma-ba
        goto 791
        endif
        if(rlu1 > s2 .and. rlu1 <= s3)then
        ik1=111
        ik2=-3212
        ic=309
!       k0- + cascade0- to pion0 + sigma0-
        goto 791
        endif
        if(rlu1 > s3 .and. rlu1 <= s4)then
        ik1=111
        ik2=-3122
        ic=310
!       k0- + cascade0- to pion0 + lambda-
        goto 791
        endif
        if(rlu1 > s4 .and. rlu1 <= s5)then
        ik1=321
        ik2=-2212
        ic=328
!       k0- + cascade0- to k+ + p-
        goto 791
        endif
        ik1=311
        ik2=-2112
        ic=330
!       k0- + cascade0- to k0 + n-
791     lc(icp,3)=ik1
        lc(icp,4)=ik2
        lc(icp,5)=ic
        tw(icp)=si16/cskn/10D0
        goto 10
        endif

!       K+ + Omega-
        if(kl == 321 .and. kl1 == 3334)then
        if(isinel(331) == 0)then
        si1=0D0
        goto 792
        endif
        ik1=211
        ik2=3312
        call srev(kl,kl1,ik1,ik2,ss,ilo1,fac, &
         0.5,0.,0.,1.5,1.0,0.5,0.,0.5,1.)
        the=pmas(pycomp(211),1)+pmas(pycomp(3312),1)
        si1=s1724(ss,ilo1,0,the)*fac
792     if(isinel(336) == 0)then
        si2=0D0
        goto 793
        endif
        ik1=111
        ik2=3322
        call srev(kl,kl1,ik1,ik2,ss,ilo2,fac, &
         0.5,0.,0.,1.5,1.0,0.5,0.,0.5,1.)
        the=pmas(pycomp(111),1)+pmas(pycomp(3322),1)
        si2=s1724(ss,ilo2,0,the)*fac
793     if(ilo1 == 0.and.ilo2 == 0) goto 13
        if(si1 < 1D-6.and.si2 < 1D-6) goto 13
        sit=si1+si2
        s1=si1/sit
        rlus=pyr(1)
        if(rlus <= s1)then
        ik1=211
        ik2=3312
        ic=331
!       K+ + Omega- to pion+ + cascade-
        goto 794
        endif
        ik1=111
        ik2=3322
        ic=336
!       K+ + Omega- to pion0 + cascade0
794     lc(icp,3)=ik1
        lc(icp,4)=ik2
        lc(icp,5)=ic
        tw(icp)=sit/cskn/10D0
        goto 10
        endif

!       K- + Omega-ba
        if(kl == -321 .and. kl1 == -3334)then
        if(isinel(333) == 0)then
        si1=0D0
        goto 795
        endif
        ik1=-211
        ik2=-3312
        call srev(kl,kl1,ik1,ik2,ss,ilo1,fac, &
         0.5,0.,0.,1.5,1.0,0.5,0.,0.5,1.)
        the=pmas(pycomp(-211),1)+pmas(pycomp(-3312),1)
        si1=s1724(ss,ilo1,0,the)*fac
795     if(isinel(338) == 0)then
        si2=0D0
        goto 796
        endif
        ik1=111
        ik2=-3322
        call srev(kl,kl1,ik1,ik2,ss,ilo2,fac, &
         0.5,0.,0.,1.5,1.0,0.5,0.,0.5,1.)
        the=pmas(pycomp(111),1)+pmas(pycomp(-3322),1)
        si2=s1724(ss,ilo2,0,the)*fac
796     if(ilo1 == 0.and.ilo2 == 0) goto 13
        if(si1 < 1D-6.and.si2 < 1D-6) goto 13
        sit=si1+si2
        s1=si1/sit
        rlus=pyr(1)
        if(rlus <= s1)then
        ik1=-211
        ik2=-3312
        ic=333
!       K- + Omega-ba to pion- + cascade-ba
        goto 797
        endif
        ik1=111
        ik2=-3322
        ic=338
!       K- + Omega-ba to pion0 + cascade0-
797     lc(icp,3)=ik1
        lc(icp,4)=ik2
        lc(icp,5)=ic
        tw(icp)=sit/cskn/10D0
        goto 10
        endif

!       k0 + Omega-
        if(kl == 311 .and. kl1 == 3334)then
        if(isinel(332) == 0)then
        si1=0D0
        goto 798
        endif
        ik1=111
        ik2=3312
        call srev(kl,kl1,ik1,ik2,ss,ilo1,fac, &
         0.5,0.,0.,1.5,1.0,0.5,0.,0.5,1.)
        the=pmas(pycomp(111),1)+pmas(pycomp(3312),1)
        si1=s1724(ss,ilo1,0,the)*fac
798     if(isinel(335) == 0)then
        si2=0D0
        goto 799
        endif
        ik1=-211
        ik2=3322
        call srev(kl,kl1,ik1,ik2,ss,ilo2,fac, &
         0.5,0.,0.,1.5,1.0,0.5,0.,0.5,1.)
        the=pmas(pycomp(-211),1)+pmas(pycomp(3322),1)
        si2=s1724(ss,ilo2,0,the)*fac
799     if(ilo1 == 0.and.ilo2 == 0) goto 13
        if(si1 < 1D-6.and.si2 < 1D-6) goto 13
        sit=si1+si2
        s1=si1/sit
        rlus=pyr(1)
        if(rlus <= s1)then
        ik1=111
        ik2=3312
        ic=332
!       k0 + Omega- to pion0 + cascade-
        goto 800
        endif
        ik1=-211
        ik2=3322
        ic=335
!       k0 + Omega- to pion- + cascade0
800     lc(icp,3)=ik1
        lc(icp,4)=ik2
        lc(icp,5)=ic
        tw(icp)=sit/cskn/10D0
        goto 10
        endif

!       k0- + Omega-ba
        if(kl == -311 .and. kl1 == -3334)then
        if(isinel(334) == 0)then
        si1=0D0
        goto 801
        endif
        ik1=111
        ik2=-3312
        call srev(kl,kl1,ik1,ik2,ss,ilo1,fac, &
         0.5,0.,0.,1.5,1.0,0.5,0.,0.5,1.)
        the=pmas(pycomp(111),1)+pmas(pycomp(-3312),1)
        si1=s1724(ss,ilo1,0,the)*fac
801     if(isinel(337) == 0)then
        si2=0D0
        goto 802
        endif
        ik1=211
        ik2=-3322
        call srev(kl,kl1,ik1,ik2,ss,ilo2,fac, &
         0.5,0.,0.,1.5,1.0,0.5,0.,0.5,1.)
        the=pmas(pycomp(211),1)+pmas(pycomp(-3322),1)
        si2=s1724(ss,ilo2,0,the)*fac
802     if(ilo1 == 0.and.ilo2 == 0) goto 13
        if(si1 < 1D-6.and.si2 < 1D-6) goto 13
        sit=si1+si2
        s1=si1/sit
        rlus=pyr(1)
        if(rlus <= s1)then
        ik1=111
        ik2=-3312
        ic=334
!       k0- + Omega-ba to pion0 + cascade-ba
        goto 803
        endif
        ik1=211
        ik2=-3322
        ic=337
!       k0- + Omega-ba to pion+ + cascade0-
803     lc(icp,3)=ik1
        lc(icp,4)=ik2
        lc(icp,5)=ic
        tw(icp)=sit/cskn/10D0
        goto 10
        endif

!       pion+ + cascade0 to k0- + sigma+
        if(kl == 211 .and. kl1 == 3322)then
        if(isinel(314) == 0) goto 13
        the=pmas(pycomp(-311),1)+pmas(pycomp(3222),1)
        ww=s1724(ss,ilo,0,the)/cspin
        if(ilo == 0) goto 13
        lc(icp,3)=-311
        lc(icp,4)=3222
        lc(icp,5)=314
        ik3=-311
        ik4=3222
        call srev(kl,kl1,ik3,ik4,ss,ilo1,fac, &
         1.0,0.5,0.,0.5,0.5,1.0,0.,0.5,1.)
        tw(icp)=ww*fac
        goto 10
        endif

!       pion+ + cascade0-
        if(kl == 211 .and. kl1 == -3322)then
        if(isinel(319) == 0)then
        si1=0D0
        goto 804
        endif
        ik1=321
        ik2=-3212
        call srev(kl,kl1,ik1,ik2,ss,ilo1,fac, &
         1.0,0.5,0.,0.5,0.5,1.0,0.,0.5,1.)
        the=pmas(pycomp(321),1)+pmas(pycomp(-3212),1)
        si1=10D0*s1724(ss,ilo1,0,the)*fac
804     if(isinel(320) == 0)then
        si2=0D0
        goto 805
        endif
        ik1=321
        ik2=-3122
        call srev(kl,kl1,ik1,ik2,ss,ilo2,fac, &
         1.0,0.5,0.,0.5,0.5,0.0,0.,0.5,1.)
        the=pmas(pycomp(321),1)+pmas(pycomp(-3122),1)
        si2=10D0*s1724(ss,ilo2,0,the)*fac
805     if(isinel(322) == 0)then
        si3=0D0
        goto 806
        endif
        ik1=311
        ik2=-3112
        call srev(kl,kl1,ik1,ik2,ss,ilo3,fac, &
         1.0,0.5,0.,0.5,0.5,1.0,0.,0.5,1.)
        the=pmas(pycomp(311),1)+pmas(pycomp(-3112),1)
        si3=10D0*s1724(ss,ilo3,0,the)*fac
806    if(isinel(137) == 0)then
        si4=0D0
        goto 1806
        endif
        the=pmas(pycomp(-311),1)+pmas(pycomp(-3334),1)
        si4=s1724(ss,ilo4,0,the)
1806    if(ilo1 == 0 .and. ilo2 == 0 .and. ilo3 == 0.and. &
         ilo4 == 0) goto 13
        if(si1 < 1D-6.and.si2 < 1D-6.and.si3 < 1D-6.and. &
         si4 < 1D-6) goto 13
        si12=si1+si2
        si13=si12+si3
        si14=si13+si4
        s1=si1/si14
        s2=si12/si14
        s3=si13/si14
        rlu1=pyr(1)
        if(rlu1 <= s1)then
        ik1=321
        ik2=-3212
        ic=319
!       pion+ + cascade0- to k+ + sigma0-
        goto 807
        endif
        if(rlu1 > s1 .and. rlu1 <= s2)then
        ik1=321
        ik2=-3122
        ic=320
!       pion+ + cascade0- to k+ + lambda-
        goto 807
        endif
       if(rlu1 > s2 .and. rlu1 <= s3)then
        ik1=311
        ik2=-3112
        ic=322
!       pion+ + cascade0- to k0 + sigma-ba
        goto 807
        endif
        ik1=-311
        ik2=-3334
        ic=137
!       pion+ + cascade0- to k0- + Omega-ba
807     lc(icp,3)=ik1
        lc(icp,4)=ik2
        lc(icp,5)=ic
        tw(icp)=si14/cspin/10D0
        goto 10
        endif

!       pion- + cascade0
        if(kl == -211 .and. kl1 == 3322)then
        if(isinel(312) == 0)then
        si1=0D0
        goto 808
        endif
        ik1=-321
        ik2=3212
        call srev(kl,kl1,ik1,ik2,ss,ilo1,fac, &
         1.0,0.5,0.,0.5,0.5,1.0,0.,0.5,1.)
        the=pmas(pycomp(-321),1)+pmas(pycomp(3212),1)
        si1=10D0*s1724(ss,ilo1,0,the)*fac
808     if(isinel(313) == 0)then
        si2=0D0
        goto 809
        endif
        ik1=-321
        ik2=3122
        call srev(kl,kl1,ik1,ik2,ss,ilo2,fac, &
         1.0,0.5,0.,0.5,0.5,0.0,0.,0.5,1.)
        the=pmas(pycomp(-321),1)+pmas(pycomp(3122),1)
        si2=10D0*s1724(ss,ilo2,0,the)*fac
809     if(isinel(315) == 0)then
        si3=0D0
        goto 810
        endif
        ik1=-311
        ik2=3112
        call srev(kl,kl1,ik1,ik2,ss,ilo3,fac, &
         1.0,0.5,0.,0.5,0.5,1.0,0.,0.5,1.)
        the=pmas(pycomp(-311),1)+pmas(pycomp(3112),1)
        si3=10D0*s1724(ss,ilo3,0,the)*fac
810     if(isinel(135) == 0)then
        si4=0D0
        goto 1810
        endif
        the=pmas(pycomp(311),1)+pmas(pycomp(3334),1)
        si4=s1724(ss,ilo4,0,the)

1810    if(ilo1 == 0 .and. ilo2 == 0 .and. ilo3 == 0.and. &
         ilo4 == 0) goto 13
        if(si1 < 1D-6.and.si2 < 1D-6.and.si3 < 1D-6.and. &
         si4 < 1D-6) goto 13
        si12=si1+si2
        si13=si12+si3
        si14=si13+si4
        s1=si1/si14
        s2=si12/si14
        s3=si13/si14
        rlu1=pyr(1)
        if(rlu1 <= s1)then
        ik1=-321
        ik2=3212
        ic=312
!       pion- + cascade0 to k- + sigma0
        goto 811
        endif
        if(rlu1 > s1 .and. rlu1 <= s2)then
        ik1=-321
        ik2=3122
        ic=313
!       pion- + cascade0 to k- + lambda
        goto 811
        endif
        if(rlu1 > s2 .and. rlu1 <= s3)then
        ik1=-311
        ik2=3112
        ic=315
!       pion- + cascade0 to k0- + sigma-
        goto 811
        endif
        ik1=311
        ik2=3334
        ic=135
!       pion- + cascade0 to k0 + Omega-
811     lc(icp,3)=ik1
        lc(icp,4)=ik2
        lc(icp,5)=ic
        tw(icp)=si14/cspin/10D0
        goto 10
        endif

!       pion- + cascade0- to k0 + sigma+ba
        if(kl == -211 .and. kl1 == -3322)then
        if(isinel(321) == 0) goto 13
        the=pmas(pycomp(311),1)+pmas(pycomp(-3222),1)
        ww=s1724(ss,ilo,0,the)/cspin
        if(ilo == 0) goto 13
        lc(icp,3)=311
        lc(icp,4)=-3222
        lc(icp,5)=321
        ik3=311
        ik4=-3222
        call srev(kl,kl1,ik3,ik4,ss,ilo1,fac, &
         1.0,0.5,0.,0.5,0.5,1.0,0.,0.5,1.)
        tw(icp)=ww*fac
        goto 10
        endif

!       pion0 + cascade0
        if(kl == 111 .and. kl1 == 3322)then
        if(isinel(311) == 0)then
        si1=0D0
        goto 812
        endif
        ik1=-321
        ik2=3222
        call srev(kl,kl1,ik1,ik2,ss,ilo1,fac, &
         1.0,0.5,0.,0.5,0.5,1.0,0.,0.5,1.)
        the=pmas(pycomp(-321),1)+pmas(pycomp(3222),1)
        si1=10D0*s1724(ss,ilo1,0,the)*fac
812     if(isinel(316) == 0)then
        si2=0D0
        goto 813
        endif
        ik1=-311
        ik2=3212
        call srev(kl,kl1,ik1,ik2,ss,ilo2,fac, &
         1.0,0.5,0.,0.5,0.5,1.0,0.,0.5,1.)
        the=pmas(pycomp(-311),1)+pmas(pycomp(3212),1)
        si2=10D0*s1724(ss,ilo2,0,the)*fac
813     if(isinel(317) == 0)then
        si3=0D0
        goto 814
        endif
        ik1=-311
        ik2=3122
        call srev(kl,kl1,ik1,ik2,ss,ilo3,fac, &
         1.0,0.5,0.,0.5,0.5,0.0,0.,0.5,1.)
        the=pmas(pycomp(-311),1)+pmas(pycomp(3122),1)
        si3=10D0*s1724(ss,ilo3,0,the)*fac
814     if(isinel(136) == 0)then
        si4=0D0
        goto 1814
        endif
        the=pmas(pycomp(111),1)+pmas(pycomp(2214),1)
        si4=s1724(ss,ilo4,0,the)

1814    if(ilo1 == 0 .and. ilo2 == 0 .and. ilo3 == 0 &
         .and.ilo4 == 0) goto 13
        if(si1 < 1D-6.and.si2 < 1D-6.and.si3 < 1D-6.and. &
         si4 < 1D-6) goto 13
        si12=si1+si2
        si13=si12+si3
        si14=si13+si4
        s1=si1/si14
        s2=si12/si14
        s3=si13/si14
        rlu1=pyr(1)
        if(rlu1 <= s1)then
        ik1=-321
        ik2=3222
        ic=311
!       pion0 + cascade0 to k- + sigma+
        goto 815
        endif
        if(rlu1 > s1 .and. rlu1 <= s2)then
        ik1=-311
        ik2=3212
        ic=316
!       pion0 + cascade0 to k0- + sigma0
        goto 815
        endif
        if(rlu1 > s2 .and. rlu1 <= s3)then
        ik1=-311
        ik2=3122
        ic=317
!       pion0 + cascade0 to k0- + lambda
        goto 815
        endif
        ik1=321
        ik2=3334
        ic=136
!       pion0 + cascade0 to k+ + Omega-
815     lc(icp,3)=ik1
        lc(icp,4)=ik2
        lc(icp,5)=ic
        tw(icp)=si14/cspin/10D0
        goto 10
        endif

!       pion0 + cascade0-
        if(kl == 111 .and. kl1 == -3322)then
        if(isinel(318) == 0)then
        si1=0D0
        goto 816
        endif
        ik1=321
        ik2=-3222
        call srev(kl,kl1,ik1,ik2,ss,ilo1,fac, &
         1.0,0.5,0.,0.5,0.5,1.0,0.,0.5,1.)
        the=pmas(pycomp(321),1)+pmas(pycomp(-3222),1)
        si1=10D0*s1724(ss,ilo1,0,the)*fac
816     if(isinel(323) == 0)then
        si2=0D0
        goto 817
        endif
        ik1=311
        ik2=-3212
        call srev(kl,kl1,ik1,ik2,ss,ilo2,fac, &
         1.0,0.5,0.,0.5,0.5,1.0,0.,0.5,1.)
        the=pmas(pycomp(311),1)+pmas(pycomp(-3212),1)
        si2=10D0*s1724(ss,ilo2,0,the)*fac
817     if(isinel(324) == 0)then
        si3=0D0
        goto 818
        endif
        ik1=311
        ik2=-3122
        call srev(kl,kl1,ik1,ik2,ss,ilo3,fac, &
         1.0,0.5,0.,0.5,0.5,0.0,0.,0.5,1.)
        the=pmas(pycomp(311),1)+pmas(pycomp(-3122),1)
        si3=10D0*s1724(ss,ilo3,0,the)*fac
818     if(isinel(138) == 0)then
        si4=0D0
        goto 1818
        endif
        the=pmas(pycomp(-321),1)+pmas(pycomp(-3334),1)
        si4=s1724(ss,ilo4,0,the)

1818    if(ilo1 == 0 .and. ilo2 == 0 .and. ilo3 == 0.and. &
         ilo4 == 0) goto 13
        if(si1 < 1D-6.and.si2 < 1D-6.and.si3 < 1D-6.and. &
         si4 < 1D-6) goto 13
        si12=si1+si2
        si13=si12+si3
        si14=si13+si4
        s1=si1/si14
        s2=si12/si14
        s3=si13/si14
        rlu1=pyr(1)
        if(rlu1 <= s1)then
        ik1=321
        ik2=-3222
        ic=318
!       pion0 + cascade0- to k+ + sigma+ba
        goto 819
        endif
        if(rlu1 > s1 .and. rlu1 <= s2)then
        ik1=311
        ik2=-3212
        ic=323
!       pion0 + cascade0- to k0 + sigma0-
        goto 819
        endif
        if(rlu1 > s2 .and. rlu1 <= s3)then
        ik1=311
        ik2=-3122
        ic=324
!       pion0 + cascade0- to k0 + lambda-
        goto 819
        endif
        ik1=-321
        ik2=-3334
        ic=138
!       pion0 + cascade0- to k- + Omega-ba
819     lc(icp,3)=ik1
        lc(icp,4)=ik2
        lc(icp,5)=ic
        tw(icp)=si14/cspin/10D0
        goto 10
        endif

!       pion+ + Omega- to k0- + cascade0
        if(kl == 211 .and. kl1 == 3334)then
        if(isinel(342) == 0) goto 13
        the=pmas(pycomp(-311),1)+pmas(pycomp(3322),1)
        ww=s1724(ss,ilo,0,the)/cspin
        if(ilo == 0) goto 13
        lc(icp,3)=-311
        lc(icp,4)=3322
        lc(icp,5)=342
        ik3=-311
        ik4=3322
        call srev(kl,kl1,ik3,ik4,ss,ilo1,fac, &
         1.0,0.0,0.,1.5,0.5,0.5,0.,0.5,1.)
        tw(icp)=ww*fac
        goto 10
        endif

!       pion+ + Omega-ba to k+ + cascade-ba
        if(kl == 211 .and. kl1 == -3334)then
        if(isinel(343) == 0) goto 13
        the=pmas(pycomp(321),1)+pmas(pycomp(-3312),1)
        ww=s1724(ss,ilo,0,the)/cspin
        if(ilo == 0) goto 13
        lc(icp,3)=321
        lc(icp,4)=-3312
        lc(icp,5)=343
        ik3=321
        ik4=-3312
        call srev(kl,kl1,ik3,ik4,ss,ilo1,fac, &
         1.0,0.0,0.,1.5,0.5,0.5,0.,0.5,1.)
        tw(icp)=ww*fac
        goto 10
        endif

!       pion- + Omega- to k- + cascade-
        if(kl == -211 .and. kl1 == 3334)then
        if(isinel(339) == 0) goto 13
        the=pmas(pycomp(-321),1)+pmas(pycomp(3312),1)
        ww=s1724(ss,ilo,0,the)/cspin
        if(ilo == 0) goto 13
        lc(icp,3)=-321
        lc(icp,4)=3312
        lc(icp,5)=339
        ik3=-321
        ik4=3312
        call srev(kl,kl1,ik3,ik4,ss,ilo1,fac, &
         1.0,0.0,0.,1.5,0.5,0.5,0.,0.5,1.)
        tw(icp)=ww*fac
        goto 10
        endif

!       pion- + Omega-ba to k0 + cascade0-
        if(kl == -211 .and. kl1 == -3334)then
        if(isinel(346) == 0) goto 13
        the=pmas(pycomp(311),1)+pmas(pycomp(-3322),1)
        ww=s1724(ss,ilo,0,the)/cspin
        if(ilo == 0) goto 13
        lc(icp,3)=311
        lc(icp,4)=-3322
        lc(icp,5)=346
        ik3=311
        ik4=-3322
        call srev(kl,kl1,ik3,ik4,ss,ilo1,fac, &
         1.0,0.0,0.,1.5,0.5,0.5,0.,0.5,1.)
        tw(icp)=ww*fac
        goto 10
        endif

!       pion0 + Omega-
        if(kl == 111 .and. kl1 == 3334)then
        if(isinel(340) == 0)then
        si1=0D0
        goto 820
        endif
        ik1=-311
        ik2=3312
        call srev(kl,kl1,ik1,ik2,ss,ilo1,fac, &
         1.0,0.0,0.,1.5,0.5,0.5,0.,0.5,1.)
        the=pmas(pycomp(-311),1)+pmas(pycomp(3312),1)
        si1=s1724(ss,ilo1,0,the)*fac
820     if(isinel(341) == 0)then
        si2=0D0
        goto 821
        endif
        ik1=-321
        ik2=3322
        call srev(kl,kl1,ik1,ik2,ss,ilo2,fac, &
         1.0,0.0,0.,1.5,0.5,0.5,0.,0.5,1.)
        the=pmas(pycomp(-321),1)+pmas(pycomp(3322),1)
        si2=s1724(ss,ilo2,0,the)*fac
821     if(ilo1 == 0 .and. ilo2 == 0) goto 13
        if(si1 < 1D-6.and.si2 < 1D-6) goto 13
        si12=si1+si2
        s1=si1/si12
        rlu1=pyr(1)
        if(rlu1 <= s1)then
        ik1=-311
        ik2=3312
        ic=340
!       pion0 + Omega- to k0- + cascade-
        goto 822
        endif
        ik1=-321
        ik2=3322
        ic=341
!       pion0 + Omega- to k- + cascade0
822     lc(icp,3)=ik1
        lc(icp,4)=ik2
        lc(icp,5)=ic
        tw(icp)=si12/cspin   ! 280224 Lei si13 -> si12
        goto 10
        endif

!       pion0 + Omega-ba
        if(kl == 111 .and. kl1 == -3334)then
        if(isinel(344) == 0)then
        si1=0D0
        goto 823
        endif
        ik1=311
        ik2=-3312
        call srev(kl,kl1,ik1,ik2,ss,ilo1,fac, &
         1.0,0.0,0.,1.5,0.5,0.5,0.,0.5,1.)
        the=pmas(pycomp(311),1)+pmas(pycomp(-3312),1)
        si1=s1724(ss,ilo1,0,the)*fac
823     if(isinel(345) == 0)then
        si2=0D0
        goto 824
        endif
        ik1=321
        ik2=-3322
        call srev(kl,kl1,ik1,ik2,ss,ilo2,fac, &
         1.0,0.0,0.,1.5,0.5,0.5,0.,0.5,1.)
        the=pmas(pycomp(321),1)+pmas(pycomp(-3322),1)
        si2=s1724(ss,ilo2,0,the)*fac
824     if(ilo1 == 0 .and. ilo2 == 0) goto 13
        if(si1 < 1D-6.and.si2 < 1D-6) goto 13
        si12=si1+si2
        s1=si1/si12
        rlu1=pyr(1)
        if(rlu1 <= s1)then
        ik1=311
        ik2=-3312
        ic=344
!       pion0 + Omega-ba to k0 + cascade-ba
        goto 825
        endif
        ik1=321
        ik2=-3322
        ic=345
!       pion0 + Omega-ba to k+ + cascade0-
825     lc(icp,3)=ik1
        lc(icp,4)=ik2
        lc(icp,5)=ic
        tw(icp)=si12/cspin
        goto 10
        endif
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!       pion+ + delta+ to pi+ + p
        if((kl == 211 .and. kl1 == 2214).or. &
         (kl1 == 211 .and. kl == 2214))then
        if(isinel(351) == 0) goto 13
        ik3=211
        ik4=2212
        call srev(kl,kl1,ik3,ik4,ss,ilo,fac, &
         1.0,1.5,0.,1.5,1.0,0.5,0.,0.5,1.)
        if(ilo == 0) goto 13
        lc(icp,3)=211
        lc(icp,4)=2212
        lc(icp,5)=351
        ww=sdelta(ss,ilo1,1,0.d0)/cspin/10D0
        tw(icp)=ww*fac
        goto 10
        endif
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!       pion0 + delta++ to pi+ + p
        if((kl == 111 .and. kl1 == 2224).or. &
         (kl1 == 111 .and. kl == 2224))then
        if(isinel(352) == 0) goto 13
        ik3=211
        ik4=2212
        call srev(kl,kl1,ik3,ik4,ss,ilo,fac, &
         1.0,1.5,0.,1.5,1.0,0.5,0.,0.5,1.)
        if(ilo == 0) goto 13
        lc(icp,3)=211
        lc(icp,4)=2212
        lc(icp,5)=352
        ww=sdelta(ss,ilo1,1,0.d0)/cspin/10D0
        tw(icp)=ww*fac
        goto 10
        endif
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!       pion+ + delta-
        if((kl == 211 .and. kl1 == 1114).or. &
         (kl1 == 211 .and. kl == 1114))then
        if(isinel(347) == 0)then
        si1=0D0
        goto 1000
        endif
        ik1=-211
        ik2=2212
        call srev(kl,kl1,ik1,ik2,ss,ilo1,fac, &
         1.0,1.5,0.,1.5,1.0,0.5,0.,0.5,1.)

        if(ilo1 == 0)then
        si1=0D0
        else
        si1=sdelta(ss,ilo1,1,0.d0)*WEIGH(19)*fac
        endif
!       cross section of pion+ + delta- to pi- + p

1000    if(isinel(348) == 0)then
        si2=0D0
        goto 1002
        endif
        ik1=111
        ik2=2112
        call srev(kl,kl1,ik1,ik2,ss,ilo2,fac, &
         1.0,1.5,0.,1.5,1.0,0.5,0.,0.5,1.)
        if(ilo2 == 0)then
        si2=0D0
        else
! since threshold energy is check in srev() we set the=0D0 here in
! sdelta()
        si2=sdelta(ss,ilo2,1,0.d0)*WEIGH(19)*fac
        endif
!       cross section of pion+ + delta- to pi0 + n

1002    if(ilo1 == 0 .and. ilo2 == 0) goto 13
        if(si1 < 1D-6.and.si2 < 1D-6) goto 13
        si12=si1+si2
        s1=si1/si12
        rlu1=pyr(1)
        if(rlu1 <= s1)then
        ik1=-211
        ik2=2212
        ic=347
!       cross section of pion+ + delta- to pi- + p
        goto 1004
        endif
        ik1=111
        ik2=2112
        ic=348
!       cross section of pion+ + delta- to pi0 + n
1004    lc(icp,3)=ik1
        lc(icp,4)=ik2
        lc(icp,5)=ic
        tw(icp)=si12/cspin/10D0
        goto 10
        endif
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!       pion+ + delta0
        if((kl == 211 .and. kl1 == 2114).or. &
         (kl1 == 211 .and. kl == 2114))then
        if(isinel(349) == 0)then
        si1=0D0
        goto 1006
        endif
        ik1=211
        ik2=2112
        call srev(kl,kl1,ik1,ik2,ss,ilo1,fac, &
         1.0,1.5,0.,1.5,1.0,0.5,0.,0.5,1.)

        if(ilo1 == 0)then
        si1=0D0
        else
        si1=sdelta(ss,ilo1,1,0.d0)*WEIGH(19)*fac
        endif
!       cross section of pion+ + delta0 to pi+ + n

1006    if(isinel(350) == 0)then
        si2=0D0
        goto 1008
        endif
        ik1=111
        ik2=2212
        call srev(kl,kl1,ik1,ik2,ss,ilo2,fac, &
         1.0,1.5,0.,1.5,1.0,0.5,0.,0.5,1.)
        if(ilo2 == 0)then
        si2=0D0
        else
! since threshold energy is check in srev() we set the=0D0 here in
! sdelta()
        si2=sdelta(ss,ilo2,1,0.d0)*WEIGH(19)*fac
        endif
!       cross section of pion+ + delta0 to pi0 + p

1008    if(ilo1 == 0 .and. ilo2 == 0) goto 13
        if(si1 < 1D-6.and.si2 < 1D-6) goto 13
        si12=si1+si2
        s1=si1/si12
        rlu1=pyr(1)
        if(rlu1 <= s1)then
        ik1=211
        ik2=2112
        ic=349
!       cross section of pion+ + delta0 to pi+ + n
        goto 1010
        endif
        ik1=111
        ik2=2212
        ic=350
!       cross section of pion+ + delta0 to pi0 + p
1010    lc(icp,3)=ik1
        lc(icp,4)=ik2
        lc(icp,5)=ic
        tw(icp)=si12/cspin/10D0
        goto 10
        endif

!       pion0 + delta+
        if((kl == 111 .and. kl1 == 2214).or. &
         (kl1 == 111 .and. kl == 2214))then
        if(isinel(353) == 0)then
        si1=0D0
        goto 1012
        endif
        ik1=111
        ik2=2212
        call srev(kl,kl1,ik1,ik2,ss,ilo1,fac, &
         1.0,1.5,0.,1.5,1.0,0.5,0.,0.5,1.)

        if(ilo1 == 0)then
        si1=0D0
        else
        si1=sdelta(ss,ilo1,1,0.d0)*WEIGH(19)*fac
        endif
!       cross section of pion0 + delta+ to pi0 + p

1012    if(isinel(354) == 0)then
        si2=0D0
        goto 1014
        endif
        ik1=211
        ik2=2112
        call srev(kl,kl1,ik1,ik2,ss,ilo2,fac, &
         1.0,1.5,0.,1.5,1.0,0.5,0.,0.5,1.)
        if(ilo2 == 0)then
        si2=0D0
        else
! since threshold energy is check in srev() we set the=0D0 here in
! sdelta()
        si2=sdelta(ss,ilo2,1,0.d0)*WEIGH(19)*fac
        endif
!       cross section of pion0 + delta+ to pi+ + n

1014    if(ilo1 == 0 .and. ilo2 == 0) goto 13
        if(si1 < 1D-6.and.si2 < 1D-6) goto 13
        si12=si1+si2
        s1=si1/si12
        rlu1=pyr(1)
        if(rlu1 <= s1)then
        ik1=111
        ik2=2212
        ic=353
!       cross section of pion0 + delta+ to pi0 + p
        goto 1016
        endif
        ik1=211
        ik2=2112
        ic=354
!       cross section of pion0 + delta+ to pi+ + n
1016    lc(icp,3)=ik1
        lc(icp,4)=ik2
        lc(icp,5)=ic
        tw(icp)=si12/cspin/10D0
        goto 10
        endif
!       pion0 + delta0
        if((kl == 111 .and. kl1 == 2114).or. &
         (kl1 == 111 .and. kl == 2114))then
        if(isinel(355) == 0)then
        si1=0D0
        goto 1018
        endif
        ik1=111
        ik2=2112
        call srev(kl,kl1,ik1,ik2,ss,ilo1,fac, &
         1.0,1.5,0.,1.5,1.0,0.5,0.,0.5,1.)

        if(ilo1 == 0)then
        si1=0D0
        else
        si1=sdelta(ss,ilo1,1,0.d0)*WEIGH(19)*fac
        endif
!       cross section of pion0 + delta0 to pi0 + n

1018    if(isinel(356) == 0)then
        si2=0D0
        goto 1020
        endif
        ik1=-211
        ik2=2212
        call srev(kl,kl1,ik1,ik2,ss,ilo2,fac, &
         1.0,1.5,0.,1.5,1.0,0.5,0.,0.5,1.)
        if(ilo2 == 0)then
        si2=0D0
        else
! since threshold energy is check in srev() we set the=0D0 here in
! sdelta()
        si2=sdelta(ss,ilo2,1,0.d0)*WEIGH(19)*fac
        endif
!       cross section of pion0 + delta0 to pi- + p

1020    if(ilo1 == 0 .and. ilo2 == 0) goto 13
        if(si1 < 1D-6.and.si2 < 1D-6) goto 13
        si12=si1+si2
        s1=si1/si12
        rlu1=pyr(1)
        if(rlu1 <= s1)then
        ik1=111
        ik2=2112
        ic=355
!       cross section of pion0 + delta0 to pi0 + n
        goto 1022
        endif
        ik1=-211
        ik2=2212
        ic=356
!       cross section of pion0 + delta0 to pi- + p
1022    lc(icp,3)=ik1
        lc(icp,4)=ik2
        lc(icp,5)=ic
        tw(icp)=si12/cspin/10D0
        goto 10
        endif
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!       pion0 + delta- to pi- + n
        if((kl == 111 .and. kl1 == 1114).or. &
         (kl1 == 111 .and. kl == 1114))then
        if(isinel(357) == 0) goto 13
        ik3=-211
        ik4=2112
        call srev(kl,kl1,ik3,ik4,ss,ilo,fac, &
         1.0,1.5,0.,1.5,1.0,0.5,0.,0.5,1.)
        if(ilo == 0) goto 13
        lc(icp,3)=-211
        lc(icp,4)=2112
        lc(icp,5)=357
        ww=sdelta(ss,ilo1,1,0.d0)/cspin/10D0
        tw(icp)=ww*fac
        goto 10
        endif
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!       pion- + delta++
        if((kl == -211 .and. kl1 == 2224).or. &
         (kl1 == -211 .and. kl == 2224))then
        if(isinel(358) == 0)then
        si1=0D0
        goto 1024
        endif
        ik1=111
        ik2=2212
        call srev(kl,kl1,ik1,ik2,ss,ilo1,fac, &
         1.0,1.5,0.,1.5,1.0,0.5,0.,0.5,1.)

        if(ilo1 == 0)then
        si1=0D0
        else
        si1=sdelta(ss,ilo1,1,0.d0)*WEIGH(19)*fac
        endif
!       cross section of pion- + delta++ to pi0 + p

1024    if(isinel(359) == 0)then
        si2=0D0
        goto 1026
        endif
        ik1=211
        ik2=2112
        call srev(kl,kl1,ik1,ik2,ss,ilo2,fac, &
         1.0,1.5,0.,1.5,1.0,0.5,0.,0.5,1.)
        if(ilo2 == 0)then
        si2=0D0
        else
! since threshold energy is check in srev() we set the=0D0 here in
! sdelta()
        si2=sdelta(ss,ilo2,1,0.d0)*WEIGH(19)*fac
        endif
!       cross section of pion- + delta++ to pi+ + n

1026    if(ilo1 == 0 .and. ilo2 == 0) goto 13
        if(si1 < 1D-6.and.si2 < 1D-6) goto 13
        si12=si1+si2
        s1=si1/si12
        rlu1=pyr(1)
        if(rlu1 <= s1)then
        ik1=111
        ik2=2212
        ic=358
!       cross section of pion- + delta++ to pi0 + p
        goto 1028
        endif
        ik1=211
        ik2=2112
        ic=359
!       cross section of pion- + delta++ to pi+ + n
1028    lc(icp,3)=ik1
        lc(icp,4)=ik2
        lc(icp,5)=ic
        tw(icp)=si12/cspin/10D0
        goto 10
        endif
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!       pion- + delta+
        if((kl == -211 .and. kl1 == 2214).or. &
         (kl1 == -211 .and. kl == 2214))then
        if(isinel(360) == 0)then
        si1=0D0
        goto 1030
        endif
        ik1=-211
        ik2=2212
        call srev(kl,kl1,ik1,ik2,ss,ilo1,fac, &
         1.0,1.5,0.,1.5,1.0,0.5,0.,0.5,1.)

        if(ilo1 == 0)then
        si1=0D0
        else
        si1=sdelta(ss,ilo1,1,0.d0)*WEIGH(19)*fac
        endif
!       cross section of pion- + delta+ to pi- + p

1030    if(isinel(361) == 0)then
        si2=0D0
        goto 1032
        endif
        ik1=111
        ik2=2112
        call srev(kl,kl1,ik1,ik2,ss,ilo2,fac, &
         1.0,1.5,0.,1.5,1.0,0.5,0.,0.5,1.)
        if(ilo2 == 0)then
        si2=0D0
        else
! since threshold energy is check in srev() we set the=0D0 here in
! sdelta()
        si2=sdelta(ss,ilo2,1,0.d0)*WEIGH(19)*fac
        endif
!       cross section of pion- + delta+ to pi0 + n

1032    if(ilo1 == 0 .and. ilo2 == 0) goto 13
        if(si1 < 1D-6.and.si2 < 1D-6) goto 13
        si12=si1+si2
        s1=si1/si12
        rlu1=pyr(1)
        if(rlu1 <= s1)then
        ik1=-211
        ik2=2212
        ic=360
!       cross section of pion- + delta+ to pi- + p
        goto 1034
        endif
        ik1=111
        ik2=2112
        ic=361
!       cross section of pion- + delta+ to pi0 + n
1034    lc(icp,3)=ik1
        lc(icp,4)=ik2
        lc(icp,5)=ic
        tw(icp)=si12/cspin/10D0
        goto 10
        endif
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!       pion- + delta0 to pi- + n
        if((kl == -211 .and. kl1 == 2114).or. &
         (kl1 == -211 .and. kl == 2114))then
        if(isinel(362) == 0) goto 13
        ik3=-211
        ik4=2112
        call srev(kl,kl1,ik3,ik4,ss,ilo,fac, &
         1.0,1.5,0.,1.5,1.0,0.5,0.,0.5,1.)
        if(ilo == 0) goto 13
        lc(icp,3)=-211
        lc(icp,4)=2112
        lc(icp,5)=362
        ww=sdelta(ss,ilo1,1,0.d0)/cspin/10D0
        tw(icp)=ww*fac
        goto 10
        endif
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!       prho0 + n
        if((kl == 113 .and. kl1 == 2112).or. &
         (kl1 == 113 .and. kl == 2112))then
        if(isinel(363) == 0)then
        si1=0D0
        goto 1036
        endif
        ik1=-211
        ik2=2212
        call srev(kl,kl1,ik1,ik2,ss,ilo1,fac, &
         1.0,0.5,1.0,0.5,1.0,0.5,0.,0.5,1.)

        if(ilo1 == 0)then
        si1=0D0
        else
        si1=srho(ss,ilo1,1,0.d0)*WEIGH(19)*fac
        endif
!       cross section of rho0 + n to pi- + p

1036    if(isinel(364) == 0)then
        si2=0D0
        goto 1038
        endif
        ik1=111
        ik2=2112
        call srev(kl,kl1,ik1,ik2,ss,ilo2,fac, &
         1.0,0.5,1.0,0.5,1.0,0.5,0.,0.5,1.)
        if(ilo2 == 0)then
        si2=0D0
        else
! since threshold energy is check in srev() we set the=0D0 here in
! sdelta()
        si2=srho(ss,ilo2,1,0.d0)*WEIGH(19)*fac
        endif
!       cross section of rho0 + n to pi0 + n

1038    if(ilo1 == 0 .and. ilo2 == 0) goto 13
        if(si1 < 1D-6.and.si2 < 1D-6) goto 13
        si12=si1+si2
        s1=si1/si12
        rlu1=pyr(1)
        if(rlu1 <= s1)then
        ik1=-211
        ik2=2212
        ic=363
!       cross section of rho0 + n to pi- + p
        goto 1040
        endif
        ik1=111
        ik2=2112
        ic=364
!       cross section of rho0 + n to pi0 + n
1040    lc(icp,3)=ik1
        lc(icp,4)=ik2
        lc(icp,5)=ic
        tw(icp)=si12/cspin/10D0
        goto 10
        endif
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!       prho+ + n
        if((kl == 213 .and. kl1 == 2112).or. &
         (kl1 == 213 .and. kl == 2112))then
        if(isinel(366) == 0)then
        si1=0D0
        goto 1042
        endif
        ik1=111
        ik2=2212
        call srev(kl,kl1,ik1,ik2,ss,ilo1,fac, &
         1.0,0.5,1.0,0.5,1.0,0.5,0.,0.5,1.)

        if(ilo1 == 0)then
        si1=0D0
        else
        si1=srho(ss,ilo1,1,0.d0)*WEIGH(19)*fac
        endif
!       cross section of rho+ + n to pi0 + p

1042    if(isinel(367) == 0)then
        si2=0D0
        goto 1044
        endif
        ik1=211
        ik2=2112
        call srev(kl,kl1,ik1,ik2,ss,ilo2,fac, &
         1.0,0.5,1.0,0.5,1.0,0.5,0.,0.5,1.)
        if(ilo2 == 0)then
        si2=0D0
        else
! since threshold energy is check in srev() we set the=0D0 here in
! sdelta()
        si2=srho(ss,ilo2,1,0.d0)*WEIGH(19)*fac
        endif
!       cross section of rho+ + n to pi+ + n

1044    if(ilo1 == 0 .and. ilo2 == 0) goto 13
        if(si1 < 1D-6.and.si2 < 1D-6) goto 13
        si12=si1+si2
        s1=si1/si12
        rlu1=pyr(1)
        if(rlu1 <= s1)then
        ik1=111
        ik2=2212
        ic=366
!       cross section of rho+ + n to pi0 + p
        goto 1046
        endif
        ik1=211
        ik2=2112
        ic=367
!       cross section of rho+ + n to pi+ + n
1046    lc(icp,3)=ik1
        lc(icp,4)=ik2
        lc(icp,5)=ic
        tw(icp)=si12/cspin/10D0
        goto 10
        endif
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!       rho0 + p
        if((kl == 113 .and. kl1 == 2212).or. &
         (kl1 == 113 .and. kl == 2212))then
        if(isinel(368) == 0)then
        si1=0D0
        goto 1048
        endif
        ik1=111
        ik2=2212
        call srev(kl,kl1,ik1,ik2,ss,ilo1,fac, &
         1.0,0.5,1.0,0.5,1.0,0.5,0.,0.5,1.)

        if(ilo1 == 0)then
        si1=0D0
        else
        si1=srho(ss,ilo1,1,0.d0)*WEIGH(19)*fac
        endif
!       cross section of rho0 + p to pi0 + p

1048    if(isinel(369) == 0)then
        si2=0D0
        goto 1050
        endif
        ik1=211
        ik2=2112
        call srev(kl,kl1,ik1,ik2,ss,ilo2,fac, &
         1.0,0.5,1.0,0.5,1.0,0.5,0.,0.5,1.)
        if(ilo2 == 0)then
        si2=0D0
        else
! since threshold energy is check in srev() we set the=0D0 here in
! sdelta()
        si2=srho(ss,ilo2,1,0.d0)*WEIGH(19)*fac
        endif
!       cross section of rho0 + p to pi+ + n

1050    if(ilo1 == 0 .and. ilo2 == 0) goto 13
        if(si1 < 1D-6.and.si2 < 1D-6) goto 13
        si12=si1+si2
        s1=si1/si12
        rlu1=pyr(1)
        if(rlu1 <= s1)then
        ik1=111
        ik2=2212
        ic=368
!       cross section of rho0 + p to pi0 + p
        goto 1052
        endif
        ik1=211
        ik2=2112
        ic=369
!       cross section of rho0 + p to pi+ + n
1052    lc(icp,3)=ik1
        lc(icp,4)=ik2
        lc(icp,5)=ic
        tw(icp)=si12/cspin/10D0
        goto 10
        endif
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!       rho- + p
        if((kl == -213 .and. kl1 == 2212).or. &
         (kl1 == -213 .and. kl == 2212))then
        if(isinel(370) == 0)then
        si1=0D0
        goto 1054
        endif
        ik1=111
        ik2=2112
        call srev(kl,kl1,ik1,ik2,ss,ilo1,fac, &
         1.0,0.5,1.0,0.5,1.0,0.5,0.,0.5,1.)

        if(ilo1 == 0)then
        si1=0D0
        else
        si1=srho(ss,ilo1,1,0.d0)*WEIGH(19)*fac
        endif
!       cross section of rho- + p to pi0 + n

1054    if(isinel(371) == 0)then
        si2=0D0
        goto 1056
        endif
        ik1=-211
        ik2=2212
        call srev(kl,kl1,ik1,ik2,ss,ilo2,fac, &
         1.0,0.5,1.0,0.5,1.0,0.5,0.,0.5,1.)
        if(ilo2 == 0)then
        si2=0D0
        else
! since threshold energy is check in srev() we set the=0D0 here in
! sdelta()
        si2=srho(ss,ilo2,1,0.d0)*WEIGH(19)*fac
        endif
!       cross section of rho- + p to pi- + p

1056    if(ilo1 == 0 .and. ilo2 == 0) goto 13
        if(si1 < 1D-6.and.si2 < 1D-6) goto 13
        si12=si1+si2
        s1=si1/si12
        rlu1=pyr(1)
        if(rlu1 <= s1)then
        ik1=111
        ik2=2112
        ic=370
!       cross section of rho- + p to pi0 + n
        goto 1058
        endif
        ik1=-211
        ik2=2212
        ic=371
!       cross section of rho- + p to pi- + p
1058    lc(icp,3)=ik1
        lc(icp,4)=ik2
        lc(icp,5)=ic
        tw(icp)=si12/cspin/10D0
        goto 10
        endif
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!       rho+ + p to pi+ + p
        if((kl == 213 .and. kl1 == 2212).or. &
         (kl1 == 213 .and. kl == 2212))then
        if(isinel(372) == 0) goto 13
        ik3=211
        ik4=2212
        call srev(kl,kl1,ik3,ik4,ss,ilo,fac, &
         1.0,0.5,1.0,0.5,1.0,0.5,0.,0.5,1.)
        if(ilo == 0) goto 13
        lc(icp,3)=211
        lc(icp,4)=2212
        lc(icp,5)=372
        ww=srho(ss,ilo1,1,0.d0)/cspin/10D0
        tw(icp)=ww*fac
        goto 10
        endif
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!       rho- + n to pi- + n
        if((kl == -213 .and. kl1 == 2112).or. &
         (kl1 == -213 .and. kl == 2112))then
        if(isinel(365) == 0) goto 13
        ik3=-211
        ik4=2112
        call srev(kl,kl1,ik3,ik4,ss,ilo,fac, &
         1.0,0.5,1.0,0.5,1.0,0.5,0.,0.5,1.)
        if(ilo == 0) goto 13
        lc(icp,3)=-211
        lc(icp,4)=2112
        lc(icp,5)=365
        ww=srho(ss,ilo1,1,0.d0)/cspin/10D0
        tw(icp)=ww*fac
        goto 10
        endif
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!       delta++ + n to p + p
        if((kl == 2224 .and. kl1 == 2112).or. &
         (kl1 == 2224 .and. kl == 2112))then
        if(isinel(373) == 0) goto 13
        ik3=2212
        ik4=2212
        call srev(kl,kl1,ik3,ik4,ss,ilo,fac, &
         1.5,0.5,1.5,0.5,0.5,0.5,0.5,0.5,0.5)
        if(ilo == 0) goto 13
        lc(icp,3)=2212
        lc(icp,4)=2212
        lc(icp,5)=373
        ww=snn(ss,ilo1,1,0.d0)/csnn/10D0
        tw(icp)=ww*fac
        goto 10
        endif
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!       delta+ + n to p + n
        if((kl == 2214 .and. kl1 == 2112).or. &
         (kl1 == 2214 .and. kl == 2112))then
        if(isinel(374) == 0) goto 13
        ik3=2212
        ik4=2112
        call srev(kl,kl1,ik3,ik4,ss,ilo,fac, &
         1.5,0.5,1.5,0.5,0.5,0.5,0.5,0.5,1.0)
        if(ilo == 0) goto 13
        lc(icp,3)=2212
        lc(icp,4)=2112
        lc(icp,5)=374
        ww=snn(ss,ilo1,1,0.d0)/csnn/10D0
        tw(icp)=ww*fac
        goto 10
        endif
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!       delta+ + p to p + p
        if((kl == 2214 .and. kl1 == 2212).or. &
         (kl1 == 2214 .and. kl == 2212))then
        if(isinel(375) == 0) goto 13
        ik3=2212
        ik4=2212
        call srev(kl,kl1,ik3,ik4,ss,ilo,fac, &
         1.5,0.5,1.5,0.5,0.5,0.5,0.5,0.5,0.5)
        if(ilo == 0) goto 13
        lc(icp,3)=2212
        lc(icp,4)=2212
        lc(icp,5)=375
        ww=snn(ss,ilo1,1,0.d0)/csnn/10D0
        tw(icp)=ww*fac
        goto 10
        endif
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!       delta0 + p to p + n
        if((kl == 2114 .and. kl1 == 2212).or. &
         (kl1 == 2114 .and. kl == 2212))then
        if(isinel(376) == 0) goto 13
        ik3=2212
        ik4=2112
        call srev(kl,kl1,ik3,ik4,ss,ilo,fac, &
         1.5,0.5,1.5,0.5,0.5,0.5,0.5,0.5,1.0)
        if(ilo == 0) goto 13
        lc(icp,3)=2212
        lc(icp,4)=2112
        lc(icp,5)=376
        ww=snn(ss,ilo1,1,0.d0)/csnn/10D0
        tw(icp)=ww*fac
        goto 10
        endif
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!       delta0 + n to n + n
        if((kl == 2114 .and. kl1 == 2112).or. &
         (kl1 == 2114 .and. kl == 2112))then
        if(isinel(377) == 0) goto 13
        ik3=2112
        ik4=2112
        call srev(kl,kl1,ik3,ik4,ss,ilo,fac, &
         1.5,0.5,1.5,0.5,0.5,0.5,0.5,0.5,0.5)
        if(ilo == 0) goto 13
        lc(icp,3)=2112
        lc(icp,4)=2112
        lc(icp,5)=377
        ww=snn(ss,ilo1,1,0.d0)/csnn/10D0
        tw(icp)=ww*fac
        goto 10
        endif
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!       delta- + p to n + n
        if((kl == 1114 .and. kl1 == 2212).or. &
         (kl1 == 1114 .and. kl == 2212))then
        if(isinel(378) == 0) goto 13
        ik3=2112
        ik4=2112
        call srev(kl,kl1,ik3,ik4,ss,ilo,fac, &
         1.5,0.5,1.5,0.5,0.5,0.5,0.5,0.5,0.5)
        if(ilo == 0) goto 13
        lc(icp,3)=2112
        lc(icp,4)=2112
        lc(icp,5)=378
        ww=snn(ss,ilo1,1,0.d0)/csnn/10D0
        tw(icp)=ww*fac
        goto 10
        endif

!       lambda- p annihilation
        if((kl == -3122.and.kl1 == 2212) &
         .or.(kl1 == -3122.and.kl == 2212))then
        if(isinel(593) == 0) goto 13
        lc(icp,3)=323
        lc(icp,4)=223
        lc(icp,5)=593
          tw(icp)=para13/5./csnn
        goto 10
        endif

!       lambda- n annihilation
        if((kl == -3122.and.kl1 == 2112) &
         .or.(kl1 == -3122.and.kl == 2112))then
        if(isinel(594) == 0) goto 13
        tw(icp)=para13/5./csnn
!       if(ilo == 0) goto 13
        lc(icp,3)=313
        lc(icp,4)=223
        lc(icp,5)=594
        goto 10
        endif

!       sigma0- p annihilation
        if((kl == -3212.and.kl1 == 2212) &
         .or.(kl1 == -3212.and.kl == 2212))then
        if(isinel(595) == 0) goto 13
        tw(icp)=para13/5./csnn
!       if(ilo == 0) goto 13
        lc(icp,3)=323
        lc(icp,4)=223
        lc(icp,5)=595
        goto 10
        endif

!       sigma0- n annihilation
        if((kl == -3212.and.kl1 == 2112) &
         .or.(kl1 == -3212.and.kl == 2112))then
        if(isinel(596) == 0) goto 13
        tw(icp)=para13/5./csnn
!       if(ilo == 0) goto 13
        lc(icp,3)=313
        lc(icp,4)=223
        lc(icp,5)=596
        goto 10
        endif

!       p- p annihilation
        if((kl == -2212.and.kl1 == 2212) &
         .or.(kl1 == -2212.and.kl == 2212))then
        if(isinel(597) == 0) goto 13
        tw(icp)=para13/csnn
!       if(ilo == 0) goto 13
        lc(icp,3)=113
        lc(icp,4)=223
        lc(icp,5)=597
        goto 10
        endif

!       p- n annihilation
        if((kl == -2212.and.kl1 == 2112) &
         .or.(kl1 == -2212.and.kl == 2112))then
        if(isinel(598) == 0) goto 13
        tw(icp)=para13/csnn
!       if(ilo == 0) goto 13
        lc(icp,3)=-213
        lc(icp,4)=223
        lc(icp,5)=598
        goto 10
        endif

!       n- p annihilation
        if((kl == -2112.and.kl1 == 2212) &
         .or.(kl1 == -2112.and.kl == 2212))then
        if(isinel(599) == 0) goto 13
        tw(icp)=para13/csnn
!       if(ilo == 0) goto 13
        lc(icp,3)=213
        lc(icp,4)=223
        lc(icp,5)=599
        goto 10
        endif

!       n- n annihilation
        if((kl == -2112.and.kl1 == 2112) &
         .or.(kl1 == -2112.and.kl == 2112))then
        if(isinel(600) == 0) goto 13
        tw(icp)=para13/csnn
!       if(ilo == 0) goto 13
        lc(icp,3)=113
        lc(icp,4)=223
        lc(icp,5)=600
        goto 10
        endif

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!       J/psi interacts with baryon or meson
        if(kl == 443 .and. kl1 == 2112)then
        if(iabsb == 0) goto 13
        call jpsin(ss,icp,ioo,1)
        if(ioo == 0) goto 13
        goto 10
        endif

        if(kl == 443 .and. kl1 == 2212)then
        if(iabsb == 0) goto 13
        call jpsip(ss,icp,ioo,1)
        if(ioo == 0) goto 13
        goto 10
        endif

        if((kl == 443 .and. kl1 == 211) &
         .or. (kl == 211 .and. kl1 == 443))then
        if(iabsm == 0) goto 13
        call jpsip1(ss,icp,ioo,1)
        if(ioo == 0) goto 13
        goto 10
        endif

        if((kl == 443 .and. kl1 == 111) &
         .or. (kl == 111 .and. kl1 == 443))then
        if(iabsm == 0) goto 13
        call jpsip0(ss,icp,ioo,1)
        if(ioo == 0) goto 13
        goto 10
        endif

        if((kl == 443 .and. kl1 == -211) &
         .or. (kl == -211 .and. kl1 == 443))then
        if(iabsm == 0) goto 13
        call jpsip2(ss,icp,ioo,1)
        if(ioo == 0) goto 13
        goto 10
        endif

        if((kl == 443 .and. kl1 == 213) &
         .or. (kl == 213 .and. kl1 == 443))then
        if(iabsm == 0) goto 13
        call jpsir1(ss,icp,ioo,1)
        if(ioo == 0) goto 13
        goto 10
        endif

        if((kl == 443 .and. kl1 == 113) &
         .or. (kl == 113 .and. kl1 == 443))then
        if(iabsm == 0) goto 13
        call jpsir0(ss,icp,ioo,1)
        if(ioo == 0) goto 13
        goto 10
        endif

        if((kl == 443 .and. kl1 == -213) &
         .or. (kl == -213 .and. kl1 == 443))then
        if(iabsm == 0) goto 13
        call jpsir2(ss,icp,ioo,1)
        if(ioo == 0) goto 13
        goto 10
        endif

!       psi' interacts with baryon or meson
        if(kl == 100443 .and. kl1 == 2112)then
        if(iabsb == 0) goto 13
        call jpsin(ss,icp,ioo,2)
!       argument '2' here refers to the psi' induced reaction
!       corresponding onse
        if(ioo == 0) goto 13
        icp5=lc(icp,5)
        if(icp5 <= 186)then
        lc(icp,5)=icp5+14
        else
        lc(icp,5)=icp5+192
        endif
!       14 and 192 are the distance (in order number) of inelastic reaction
!        between J/psi induced reaction and psi, induced corresponding
!        reaction
        goto 10
        endif

        if(kl == 100443 .and. kl1 == 2212)then
        if(iabsb == 0) goto 13
        call jpsip(ss,icp,ioo,2)
        if(ioo == 0) goto 13
        icp5=lc(icp,5)
        if(icp5 <= 186)then
        lc(icp,5)=icp5+14
        else
        lc(icp,5)=icp5+192
        endif
        goto 10
        endif

        if((kl == 100443 .and. kl1 == 211) &
         .or. (kl == 211 .and. kl1 == 100443))then
        if(iabsm == 0) goto 13
        call jpsip1(ss,icp,ioo,2)
        if(ioo == 0) goto 13
        icp5=lc(icp,5)
        if(icp5 <= 186)then
        lc(icp,5)=icp5+14
        else
        lc(icp,5)=icp5+192
        endif
        goto 10
        endif

        if((kl == 100443 .and. kl1 == 111) &
         .or. (kl == 111 .and. kl1 == 100443))then
        if(iabsm == 0) goto 13
        call jpsip0(ss,icp,ioo,2)
        if(ioo == 0) goto 13
        icp5=lc(icp,5)
        if(icp5 <= 186)then
        lc(icp,5)=icp5+14
        else
        lc(icp,5)=icp5+192
        endif
        goto 10
        endif

        if((kl == 100443 .and. kl1 == -211) &
         .or. (kl == -211 .and. kl1 == 100443))then
        if(iabsm == 0) goto 13
        call jpsip2(ss,icp,ioo,2)
        if(ioo == 0) goto 13
        icp5=lc(icp,5)
        if(icp5 <= 186)then
        lc(icp,5)=icp5+14
        else
        lc(icp,5)=icp5+192
        endif
        goto 10
        endif

        if((kl == 100443 .and. kl1 == 213) &
         .or. (kl == 213 .and. kl1 == 100443))then
        if(iabsm == 0) goto 13
        call jpsir1(ss,icp,ioo,2)
        if(ioo == 0) goto 13
        icp5=lc(icp,5)
        if(icp5 <= 186)then
        lc(icp,5)=icp5+14
        else
        lc(icp,5)=icp5+192
        endif
        goto 10
        endif

        if((kl == 100443 .and. kl1 == 113) &
         .or. (kl == 113 .and. kl1 == 100443))then
        if(iabsm == 0) goto 13
        call jpsir0(ss,icp,ioo,2)
        if(ioo == 0) goto 13
        icp5=lc(icp,5)
        if(icp5 <= 186)then
        lc(icp,5)=icp5+14
        else
        lc(icp,5)=icp5+192
        endif
        goto 10
        endif

        if((kl == 100443 .and. kl1 == -213) &
         .or. (kl == -213 .and. kl1 == 100443))then
        if(iabsm == 0) goto 13
        call jpsir2(ss,icp,ioo,2)
        if(ioo == 0) goto 13
        icp5=lc(icp,5)
        if(icp5 <= 186)then
        lc(icp,5)=icp5+14
        else
        lc(icp,5)=icp5+192
        endif
        goto 10
        endif

!-------------------------------------------------------------------------------
!       D meson + pion/rho ( pion/rho + D meson )
        if( ( (ABS(kl)  == 411 .OR. ABS(kl)  == 421  .OR.  &
               ABS(kl)  == 413 .OR. ABS(kl)  == 423) .AND. &
              (ABS(kl1) == 211 .OR. ABS(kl1) == 111  .OR.  &
               ABS(kl1) == 213 .OR. ABS(kl1) == 113) )     &
        .OR.( (ABS(kl1) == 411 .OR. ABS(kl1) == 421  .OR.  &
               ABS(kl1) == 413 .OR. ABS(kl1) == 423) .AND. &
              (ABS(kl)  == 211 .OR. ABS(kl)  == 111  .OR.  &
               ABS(kl)  == 213 .OR. ABS(kl)  == 113) ) )then
            call DMeson_coll(kl,kl1,ss,icp,ioo)
            ! Treateda as elastic collision.
            if( ioo == 0 ) goto 13
            ! Inelastic collision.
            goto 10
!       D meson + N (the order has been specified when calling "prod")
        else if( (ABS(kl)  == 411  .OR. ABS(kl)  == 421  .OR.  &
                  ABS(kl)  == 413  .OR. ABS(kl)  == 423) .AND. &
                 (ABS(kl1) == 2212 .OR. ABS(kl1) == 2112) )then
            call Dpn_coll(kl,kl1,ss,icp,ioo)
            ! Treateda as elastic collision.
            if( ioo == 0 ) goto 13
            ! Inelastic collision.
            goto 10
        end if

        !debug
        i = l
        j = l1
        ! write(9,*) "No inel. channel is found in hadcas prod(), " &
        !         // "treated as elastic collision, " &
        !         // "l, l1, kl, kl1=", l, l1, kl, kl1
        !debug
13      do m=3,5,1
            lc(icp,m) = 0
        end do
        tw(icp) = 0D0
        return

10      if( lc(icp,5) < 1593 )then
            nchargef = PYCHGE( lc(icp,3) ) + PYCHGE( lc(icp,4) )
            if( nchargei/3 /= nchargef/3 )then
              write(9,*) "Warning, Charges were not conserved, " &
                      // "initial, final, lc(icp,5)=", &
                          nchargei/3D0, nchargef/3D0, lc(icp,5)
            end if
        end if
!-------------------------------------------------------------------------------


        return
        end



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine spipi(lc3,lc4,ss,ilo)
!       A part of prod()
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        am1 = PYMASS(lc3)
        am2 = PYMASS(lc4)
        ilo=1
        if( ss < (am1 + am2) ) ilo = 0
        return
        end



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine srev(kl,kl1,lc3,lc4,ss,ilo,fac,xii1,xii2,xsi1,xsi2, &
         xif1,xif2,xsf1,xsf2,pauli)
!!      Calculates the reverse reaction factor by the detailed balance.
!!      kl,  kl1: KF codes of two particles in initial state
!!      lc3, lc4: KF codes of two particles in final state
!       ss: CM energy
!       ilo: = 0, the CM energy is lower than threshold energy
!            = 1, succeeded.
!       xii1, xii2: the isospin of of two particles in initial state
!       xsi1, xsi2: the spin of two particles in initial state
!       xif1, xif2: the isospin of of two particles in final state
!       xsf1, xsf2: the spin of two particles in final state
!       paul1: =1 if two particles are different in final state
!              =0.5 if identical two particles in final state
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        real*4 xii1,xii2,xsi1,xsi2,xif1,xif2,xsf1,xsf2,pauli


        x_isospin_initial_1 = xii1
        x_isospin_initial_2 = xii2
        x_isospin_final_1 = xif1
        x_isospin_final_2 = xif2
        am1  = PYMASS(kl)
        am2  = PYMASS(kl1)
        am3  = PYMASS(lc3)
        am4  = PYMASS(lc4)
        am12 = am1*am1
        am22 = am2*am2
        am32 = am3*am3
        am42 = am4*am4
        ss2  = ss*ss
        ss4  = ss2*ss2
        ilo  = 1

        if( ss < (am3 + am4) )then
            ilo = 0
            fac = 0D0
        else
            pfp = ss4 - 2D0 * ss2 * (am32 + am42) + (am32 - am42)*(am32 - am42)
            pfp = pfp / 4D0 / ss2
            if( pfp < 0D0 ) pfp = 1D-10
!           pfp = SQRT(pfp) / 2D0 / ss
            pip = ss4 - 2D0 * ss2 * (am12 + am22) + (am12 - am22)*(am12 - am22)
            pip = pip / 4D0 / ss2
            if( pip < 0D0 ) pip = 1D-10
!           pip = SQRT(pip) / 2D0 / ss
!           phase = pauli * (2D0*xif1 + 1D0) * (2D0*xif2 + 1D0) &
!                         * (2D0*xsf1 + 1D0) * (2D0*xsf2 + 1D0) &
!                         / (2D0*xii1 + 1D0) * (2D0*xii2 + 1D0) &
!                         / (2D0*xsi1 + 1D0) * (2D0*xsi2 + 1D0)
            phase = pauli * (2D0*xsf1 + 1D0) * (2D0*xsf2 + 1D0) &
                          / (2D0*xsi1 + 1D0) / (2D0*xsi2 + 1D0)
            fac = phase * pfp / pip
        end if


        return
        end



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine updpli( l, l1, icp, pii, pjj, time, iia )
!!      Updates particle list after inelastic collision and
!!       truncate collision list correspondingly.
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        PARAMETER (KSZJ=80000,NSIZE=750000)
        common/sa12/ppsa(5),nchan,nsjp,sjp,taup,taujp
        common/sa1_h/nsa,non1,ksa(kszj,5),psa(kszj,5),vsa(kszj,5)
!       Note the name of the arrays in 'sa1_h' in this subroutine.
        common/sa8_h/tau(kszj),ishp(kszj)
        common/sa9_h/kfmax,kfaco(100),numb(100),non9,disbe(100,100)
        common/sa19_h/coor(3)
        common/sa20_h/t0,cspipiKK,dep,ddt,edipi,epin,ecsnn,ekn,ecspsn,ecspsm, &
         rnt,rnp,ecsspn,ecsspm
        common/ctllist_h/nctl,noinel(600),nctl0,noel
        common/ctllist_t/ lc(nsize,5), tc(nsize), tw(nsize)
        dimension pii(4),pjj(4),pp(4)


        kf1 = ksa(l,2)
        kf2 = ksa(l1,2)
!       do i=1,3
!           vsa(l,i)  = ( vsa(l,i) + vsa(l1,i) )*0.5D0
!           vsa(l1,i) = vsa(l,i)
!       end do
        if(iia < 1593) goto 1400

!       Updates the collision list for annihilation case.
        ic = l
        jc = l1
        j  = 0
        do i=1,nctl,1
            i1 = lc(i,1)
            j1 = lc(i,2)
            if(i1 == ic .or. i1 == jc) goto 4000
            if(j1 == ic .or. j1 == jc) goto 4000
!       Avoid the particle collideing with itself.
            if(i1 == j1) goto 4000
            if( (tc(i)-time) <= ddt ) goto 4000
!       Throws away the pairs with tc <= time.
            j = j + 1
            tc(j) = tc(i)
            tw(j) = tw(i)
            do m=1,5,1
                lc(j,m) = lc(i,m)
            end do
4000        continue
        end do
        do i = j+1, nctl, 1
            tc(i) = 0D0
            tw(i) = 0D0
            do m=1,5,1
                lc(i,m) = 0
            end do
        end do
        nctl = j
        goto 1200

1400    ik1 = lc(icp,3)
        ik2 = lc(icp,4)
!       Puts the scattered (produced) particles into particle list
!        & turncateS collision list correspondingly.
        ll  = l
        ll1 = l1
        kf  = ik1
        do i=1,4,1
            pp(i) = pii(i)
        end do

        do 500 i=1,2,1
!           nnn = 0
            do 600 j=1,kfmax,1
!               if( kf /= kfaco(j) ) goto 900
                if( kf /= kfaco(j) ) goto 600
                jj = numb(j) + 1
!           Updates the particle list.
!!1000          do m = nsa+n1, jj, -1
! 1000            do m = nsa, jj, -1
                do m = nsa, jj, -1
                    mm = m + 1
                    ksa(mm,2) = ksa(m,2)
                    ksa(mm,1) = 1
                    ksa(mm,3) = ksa(m,3)
                    do m1=1,5,1
                        psa(mm,m1) = psa(m,m1)
                        vsa(mm,m1) = vsa(m,m1)
                    end do
                    ishp(mm) = ishp(m)
                    tau(mm)  = tau(m)
                end do
                if( ll  >= jj ) ll  = ll  + 1
                if( ll1 >= jj ) ll1 = ll1 + 1
!           Updates the values of lc(m,1-2).
                do m=1,nctl,1
                    lc1 = lc(m,1)
                    if(lc1 >= jj) lc(m,1) = lc1 + 1
                    lc2 = lc(m,2)
                    if(lc2 >= jj) lc(m,2) = lc2 + 1
                end do
!           ll instead of jj all over the collision list.
                do m=1,nctl,1
                    lc1 = lc(m,1)
                    lc2 = lc(m,2)
                    if(lc1 == ll) lc(m,1) = jj
                    if(lc2 == ll) lc(m,2) = jj
                end do
!           Gives proper values to particle jj.
                ksa(jj,2) = kf
                ksa(jj,1) = 1
                ksa(jj,3) = 0
                do m=1,4,1
                    psa(jj,m) = pp(m)
                    vsa(jj,m) = vsa(ll,m)
                end do
                psa(jj,4) = pp(4)
                psa(jj,5) = PYMASS(kf)
                ishp(jj)  = ishp(ll)
!               Adds photon.
                if( kf == 22 .or. kf == 44 .or. kf == 55 .or. kf == 66 )then
                    tau(jj) = time
                    goto 301
                end if
!               Give zero formation time to particles produced from
!                rescatttering
!               tau(jj) = t0 * psa(jj,4) / psa(jj,5)
                tau(jj) = time + t0*taup * psa(jj,4) / psa(jj,5)
301             lc(icp,i) = jj
!               if(nnn == 1) goto 300
                do m = j, kfmax, 1
                    numb(m) = numb(m) + 1
                end do
                goto 2000
!!300           goto 200
!!900           if(j < kfmax) goto 600
!               Scattered particle with new flavor.
!               nnn = 1
!               jj = numb(kfmax) + 1
!               kfmax = kfmax + 1
!               numb(kfmax) = jj
!               kfaco(kfmax) = kf
!               goto 1000
600         continue
!           Gives proper values to particle not considered in 'kfaco'.
            jj = nsa + 1
            ksa(jj,2) = kf
            ksa(jj,1) = 1
            ksa(jj,3) = 0
            do m=1,4,1
                psa(jj,m) = pp(m)
                vsa(jj,m) = vsa(ll,m)
            end do
            psa(jj,5) = PYMASS(kf)
            ishp(jj)  = ishp(ll)
!           Gives zero formation time to particles produced from rescatttering.
!           tau(jj) = t0 * psa(jj,4) / psa(jj,5)
            tau(jj) = time + t0*taup * psa(jj,4) / psa(jj,5)
            lc(icp,i) = jj
!!200       continue
2000        nsa = nsa + 1
            if(i == 2) goto 500
            ll2 = ll
            ll  = ll1
            ll1 = ll2
            kf  = ik2
            do j=1,4,1
                pp(j) = pjj(j)
            end do
500     continue

        l  = ll1
        l1 = ll
!       Removes the scattering particles from particle list &
!        truncates the collision list correspondingly.
1200    kf = kf1
        ll = l
        do 700 i=1,2,1
            if(ll == nsa)then   !
                do i1=1,kfmax,1
                    if( kf /= kfaco(i1) ) goto 400
                    do m = i1, kfmax, 1
                        numb(m) = numb(m) - 1
                    end do
                    if(i1 > 1)then
                        numba = numb(i1)
                        do m = 1, i1-1, 1
                            if( numb(m) == numba ) numb(m) = numb(m) - 1
                        end do
                    end if
                    goto 100
400             end do
            end if   !
!           do j = ll+1, nsa+n1, 1
            do j = ll+1, nsa, 1
                jj = j - 1
                ksa(jj,2) = ksa(j,2)
                ksa(jj,1) = 1
                ksa(jj,3) = ksa(j,3)
                do m=1,5,1
                    psa(jj,m) = psa(j,m)
                    vsa(jj,m) = vsa(j,m)
                end do
                ishp(jj) = ishp(j)
                tau(jj)  = tau(j)
            end do
            do m=1,nctl,1
                lc1 = lc(m,1)
                lc2 = lc(m,2)
                if(lc1 > ll) lc(m,1) = lc1 - 1
                if(lc2 > ll) lc(m,2) = lc2 - 1
            end do
            do 800 j=1,kfmax,1
                if( kf /= kfaco(j) ) goto 800
                do m=j,kfmax,1
                    numb(m) = numb(m) - 1
                end do
                goto 100
800         continue

100         continue
            nsa = nsa - 1
            if(l1 > ll) l1 = l1 - 1
            if(i == 2) goto 700
            ll = l1
            kf = kf2
700     continue


        return
        end



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine coelas_h(ic,jc,eij,pi,pj)
!!      Performs elastic scattering.
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        PARAMETER (KSZJ=80000)
        common/sa1_h/nsa,non1,ksa(kszj,5),psa(kszj,5),vsa(kszj,5)
        dimension pi(4),pj(4)


        iic = ksa(ic,2)
        jjc = ksa(jc,2)
        d   = 3.65D0 * ( eij - PYMASS(iic) - PYMASS(jjc) )
        if( d < 1D-10 ) return
        ! changed from 0.2 to 0.5 on 100111
        ! pt = 0.2D0
        pt = 0.5D0
        A  = MIN( 10.3D0, 1.D0 / ( 1.12D0 * pt ) / ( 1.12D0 * pt ) )
        d6 = d**6
        B  = d6 * A / ( 1D0 + d6 )
        if( B < 1D-20 ) B = 1D-20

        ! Approximation
        ! pt = 0.5D0
        ! A  = MIN( 10.3D0, 1.D0 / ( 1.12D0 * pt ) / ( 1.12D0 * pt ) )
        ! B = A

        pm2 = pi(1)**2 + pi(2)**2 + pi(3)**2
        pm  = SQRT(pm2)
        t0  = - 4D0 * pm2
        if( ABS(t0) < 1D-20 )then
            cctas = 1D0
!       maximum of cos(theta_s) = 1
            goto 100
        end if
        if( ABS( B*t0 ) < 0.0001D0 )then
            ! abt = EXP( B*t0 ), B*t0 < 0
            abt = 1D0
!       else if( B*t0 < -50D0 )then
!           abt = 0D0
        else
            abt = EXP( MAX( -7.0D2, B*t0 ) )
        end if
        cc = PYR(1)
        tt1 = LOG( cc + ( 1D0 - cc ) * abt )
        if( ABS(tt1) < 1D-30 .AND. B <= 1D-20 )then
            cctas = 1D0
            goto 100
        end if
        tt = tt1 / B
        if( ABS(tt) < 1D-20 )then
            cctas = 1D0
            goto 100
        end if
        cctas = 1D0 - tt * 2D0 / t0
        if( ABS(cctas) > 1D0 )then
            cctas = SIGN( 1D0, cctas )
        end if
100     continue
        sctas = SQRT( 1D0 - cctas**2 )
        fis   = 2D0 * 3.141592653589793D0 * PYR(1)
        cfis  = COS(fis)
        sfis  = SIN(fis)
        call rotate_h( cctas, sctas, cfis, sfis, pm, pi, pj )


        return
        end



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine rotate_h(cctas,sctas,cfis,sfis,pp3,pi,pj)
!!      Performs the rotation after elastic scattering.
!       pi,pj: input, four momentum of colliding pair before scattering
!              output,four momentum of scattered particles after rotation
!       pp3: momentum modulus of pi or pj, both are equal in their cms,
!        after scattering
!       cctas,sctas,cfis,sfis: direction cosines of momentum of one of
!        scattered particle relative to the momentum
!        of corresponding particle before scattering.
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        dimension pi(4),pj(4)


        fi1   = PYANGL( pi(1), pi(2) )
        cta1  = PYANGL( pi(3), SQRT( pi(1)**2 + pi(2)**2 ) )
        cfi1  = COS(fi1)
        sfi1  = SIN(fi1)
        ccta1 = COS(cta1)
        scta1 = SIN(cta1)
        pi(1) = cfi1 * ( ccta1*sctas*cfis + scta1*cctas ) - sfi1*sctas*sfis
        pi(2) = sfi1 * ( ccta1*sctas*cfis + scta1*cctas ) + cfi1*sctas*sfis
        pi(3) = ccta1*cctas - scta1*sctas*cfis
        pi(1) = pp3 * pi(1)
        pi(2) = pp3 * pi(2)
        pi(3) = pp3 * pi(3)
        do i=1,3,1
            pj(i) = 0D0 - pi(i)
        end do


        return
        end



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine updple_h(ic,jc,b,pi,pj)
!!      Updates the particle list after elastic scattering.
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        PARAMETER (KSZJ=80000)
        common/sa1_h/nsa,non1,ksa(kszj,5),psa(kszj,5),vsa(kszj,5)
!       note the name of the arrays in 'sa2' in this subroutine
        common/sa8_h/tau(kszj),ishp(kszj)
        common/sa19_h/coor(3)
        common/sa20_h/t0,cspipiKK,dep,ddt,edipi,epin,ecsnn,ekn,ecspsn,ecspsm, &
         rnt,rnp,ecsspn,ecsspm
        dimension pi(4),pj(4),b(3)


        ilo = 1
!       ilo=1 for inverse Lorentz transformation
        call lorntz(ilo,b,pi,pj)
        do i=1,4,1
            psa(ic,i) = pi(i)
            psa(jc,i) = pj(i)
        end do


        return
        end



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine updatl_h(ic,jc,time)
!!      Updates the collision time list.
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        PARAMETER (NSIZE=750000)
        common/sa9_h/kfmax,kfaco(100),numb(100),non9,disbe(100,100)
        common/sa20_h/t0,cspipiKK,dep,ddt,edipi,epin,ecsnn,ekn,ecspsn,ecspsm, &
         rnt,rnp,ecsspn,ecsspm
        common/ctllist_h/nctl,noinel(600),nctl0,noel
        common/ctllist_t/ lc(nsize,5), tc(nsize), tw(nsize)


        time_resolution = ddt
        nctl_new = 0

!       Loops over the old colliding pairs.
        do i = 1, nctl, 1
            i1 = lc(i,1)
            j1 = lc(i,2)
            if( i1 == ic .OR. i1 == jc .OR. j1 == ic .OR. j1 == jc &
                .OR. i1 == j1 ) cycle
            ! Throws away the pairs which have tc <= time.
            if( ABS( tc(i) - time ) <= time_resolution ) cycle
            nctl_new = nctl_new + 1
            tc(nctl_new) = tc(i)
            tw(nctl_new) = tw(i)
            do m=1,5,1
                lc(nctl_new,m) = lc(i,m)
            end do
        end do

!       Loops over the particle list.
        nctl = nctl_new + 1
        n_max_considered = numb(kfmax)

        j = ic
        do i_scattered_particle = 1,2,1
            if( j <= 0 .OR. j > n_max_considered )then
                j = jc
                cycle
            end if
            do i = 1, n_max_considered, 1
            ! Forbids scattered particles colliding with each other immediately.
                if( i == ic .AND. j == jc ) cycle
                if( i == jc .AND. j == ic ) cycle
                ! Avoid the particle collideing with itself.
                if( i == j ) cycle
                call rsfilt_h( i, j, iflag )
                if( iflag == 0 ) cycle
                tc(nctl) = 0D0
                tw(nctl) = 0D0
                call tcolij_h( i, j, time, nctl )
                tci = tc(nctl)
                if( tci <= 1D-7 ) cycle
                do jj = 1, nctl-1, 1
                    tcj = tc(jj)
                    if( ABS(tcj - tci) < time_resolution )then
                        nctl = nctl - 1
                        exit
                    end if
                end do
                nctl = nctl + 1
            end do
            j = jc
        end do
        nctl = nctl - 1


        return
        end



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine rsfilt_h( l, l1, iflag )
!!      Plays the role of first range filter and guarantee the collision list
!!       is composed according to the entrance channels of considered
!!       inelastic reactions.
!!      Subroutine intdis plays the role of second range filter
!!      Collision pairs not interested can not filter through both of rsfilt
!!       and intdis.
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        PARAMETER (KSZJ=80000)
        common/sa1_h/nsa,non1,ksa(kszj,5),psa(kszj,5),vsa(kszj,5)
        common/sa8_h/tau(kszj),ishp(kszj)
        common/sa10_h/csnn,cspin,cskn,cspipi,cspsn,cspsm,rcsit,ifram, &
         iabsb,iabsm,i_sigma_AQM,csspn,csspm


        iflag = 0
        kl    = ksa(l,2)
        kl1   = ksa(l1,2)
        if( l == l1 ) return
        if( ishp(l) == 0 .or. ishp(l1) == 0 ) return
        ! Hadrons only.
        if( ABS(kl) < 100 .OR. ABS(kl1) < 100 ) return

!       Considers channels that are not well defined and treats them as elastic.
        if( i_sigma_AQM == 2 .OR. i_sigma_AQM == 4 )then
            iflag = 1
            return
        end if

!       Constraints on the direct reactions.
!       pi+ + pi-/pi0/p/n/Sigma/Lambda/Xi
        if( kl == 211 .and. (kl1 == -211 .or. kl1 == 111 .or.           &
            abs(kl1) == 2212 .or. abs(kl1) == 2112 .or. kl1 ==          &
            3112 .or. kl1 == -3122 .or. kl1 == -3222 .or. kl1           &
            == -3212 .or. kl1 == 3212 .or. kl1 == 3122 .or. kl1 == 3312 &
            .or. kl1 == -3322) ) iflag = 1
!       pi-/pi0/p/n/Sigma/Lambda/Xi + pi+
        if( kl1 == 211 .and. (kl == -211 .or. kl == 111 .or.         &
            abs(kl) == 2212 .or. abs(kl) == 2112 .or. kl ==          &
            3112 .or. kl == -3122 .or. kl == -3222 .or. kl           &
            == -3212 .or. kl == 3212 .or. kl == 3122 .or. kl == 3312 &
            .or. kl == -3322) ) iflag = 1
!       pi- + pi0/p/n/Sigma/Lambda/Xi
        if( kl == -211 .and. (kl1 == 111 .or.                       &
            abs(kl1) == 2212 .or. abs(kl1) == 2112 .or. kl1 ==      &
            -3112 .or. kl1 == 3122 .or. kl1 == 3222 .or. kl1        &
            == 3212 .or. kl1 == -3212 .or. kl1 == -3122 .or. kl1 == &
            -3312 .or. kl1 == 3322) ) iflag = 1
!       pi0/p/n/Sigma/Lambda/Xi + pi-
        if( kl1 == -211 .and. (kl == 111 .or.                    &
            abs(kl) == 2212 .or. abs(kl) == 2112 .or. kl ==      &
            -3112 .or. kl == 3122 .or. kl == 3222 .or. kl        &
            == 3212 .or. kl == -3212 .or. kl == -3122 .or. kl == &
            -3312 .or. kl == 3322) ) iflag = 1
!       pi0 + pi0/p/n/Sigma/Lambda/Xi
        if( kl == 111 .and. (kl1 == 111 .or. abs(kl1) == 2212            &
            .or. abs(kl1) == 2112 .or. abs(kl1) == 3122 .or.             &
            abs(kl1) == 3112 .or. abs(kl1) == 3212 .or. abs(kl1) == 3222 &
            .or. abs(kl1) == 3312 .or. abs(kl1) == 3322) ) iflag = 1
!       pi0/p/n/Sigma/Lambda/Xi + pi0
        if( kl1 == 111 .and. (kl == 111 .or. abs(kl) == 2212          &
            .or. abs(kl) == 2112 .or. abs(kl) == 3122 .or.            &
            abs(kl) == 3112 .or. abs(kl) == 3212 .or. abs(kl) == 3222 &
            .or. abs(kl) == 3312 .or. abs(kl) == 3322) ) iflag = 1
!       K+ + p/n/Sigma/Lambda/Xi
        if( kl == 321 .and. (kl1 == -2212 .or. kl1 == -2112 .or. &
            kl1 == -3122 .or. kl1 == -3222 .or. kl1 == -3112     &
            .or. kl1 == -3212 .or.kl1 == -3312 .or.kl1 == -3322) ) iflag = 1
!       p/n/Sigma/Lambda/Xi + K+
        if( kl1 == 321 .and. (kl == -2212 .or. kl == -2112 .or. &
            kl == -3122 .or. kl == -3222 .or. kl == -3112       &
            .or. kl == -3212 .or.kl == -3312 .or.kl == -3322) ) iflag = 1
!       K- + p/n/Sigma/Lambda/Xi
        if( kl == -321 .and. (kl1 == 2212 .or. kl1 == 2112 .or. &
            kl1 == 3122 .or. kl1 == 3222 .or. kl1 == 3112       &
            .or. kl1 == 3212 .or.kl1 == 3312 .or. kl1 == 3322) ) iflag = 1
!       p/n/Sigma/Lambda/Xi + K-
        if( kl1 == -321 .and. (kl == 2212 .or. kl == 2112 .or. &
            kl == 3122 .or. kl == 3222 .or. kl == 3112         &
            .or. kl == 3212 .or.kl == 3312 .or. kl == 3322) ) iflag = 1
!       K0 + p/n/Sigma/Lambda/Xi
        if( kl == 311 .and. (kl1 == -2212 .or. kl1 == -2112 .or. kl1 == &
            -3122 .or. kl1 == -3212 .or. kl1 == -3112 .or.kl1 == -3222  &
            .or. kl1 == -3312 .or. kl1 == -3322) ) iflag = 1
!       p/n/Sigma/Lambda/Xi + K0
        if( kl1 == 311 .and. (kl == -2212 .or. kl == -2112 .or. kl == &
            -3122 .or. kl == -3212 .or. kl == -3112 .or.kl == -3222   &
            .or. kl == -3312 .or. kl == -3322) ) iflag = 1
!       Kbar0 + p/n/Sigma/Lambda/Xi
        if( kl == -311 .and. (kl1 == 2212 .or. kl1 == 2112 .or. kl1 == &
            3122 .or. kl1 == 3212 .or. kl1 == 3112 .or.kl1 == 3222     &
            .or.kl1 == 3312 .or.kl1 == 3322) ) iflag = 1
!       p/n/Sigma/Lambda/Xi + Kbar0
        if( kl1 == -311 .and. (kl == 2212 .or. kl == 2112 .or. kl == &
            3122 .or. kl == 3212 .or. kl == 3112 .or.kl == 3222      &
            .or.kl == 3312 .or.kl == 3322) ) iflag = 1

!       Constraints on the NN scattering.
!       n0 + n0/p
        if( kl == 2112 .and. (kl1 == 2112 .or. kl1 == 2212) ) iflag = 1
!       p + n0/p
        if( kl == 2212 .and. (kl1 == 2112 .or. kl1 == 2212) ) iflag = 1

!       Constraints on the annihilation reactions.
!       Sigmabar0 + p/n0
        if( kl  == -3212 .and. (kl1 == 2212 .or. kl1 == 2112) ) iflag = 1
!       Lambdabar0 + p/n0
        if( kl  == -3122 .and. (kl1 == 2212 .or. kl1 == 2112) ) iflag = 1
!       p/n0 + Sigmabar0
        if( kl1 == -3212 .and. (kl  == 2212 .or. kl  == 2112) ) iflag = 1
!       p/n0 + Lambdabar0
        if( kl1 == -3122 .and. (kl  == 2212 .or. kl  == 2112) ) iflag = 1
!       pbar + p/n0
        if( kl  == -2212 .and. (kl1 == 2212 .or. kl1 == 2112) ) iflag = 1
!       p/n0 + pbar
        if( kl1 == -2212 .and. (kl  == 2212 .or. kl  == 2112) ) iflag = 1
!       n0bar + p/n0
        if( kl  == -2112 .and. (kl1 == 2212 .or. kl1 == 2112) ) iflag = 1
!       p/n0 + n0bar
        if( kl1 == -2112 .and. (kl  == 2212 .or. kl  == 2112) ) iflag = 1

!       Constraints on the reverse reactions.
!       pi+ + Sigma/Lambda/Xi/Omega/Delta
        if( kl == 211 .and. (kl1 == 3112 .or. abs(kl1) == 3122 .or. kl1 &
            == -3222 .or. abs(kl1) == 3212 .or. abs(kl1) == 3312 .or.   &
            abs(kl1) == 3322 .or. abs(kl1) == 3334.or.kl1 == 1114       &
            .or.kl1 == 2114.or.kl1 == 2214) ) iflag = 1
!       Sigma/Lambda/Xi/Omega/Delta + pi+
        if( kl1 == 211 .and. (kl == 3112 .or. abs(kl) == 3122 .or. kl &
            == -3222 .or. abs(kl) == 3212 .or. abs(kl) == 3312 .or.   &
            abs(kl) == 3322 .or. abs(kl) == 3334.or.kl == 1114        &
            .or.kl == 2114.or.kl == 2214) ) iflag = 1
!       pi- + Sigma/Lambda/Xi/Omega/Delta
        if( kl == -211 .and. (kl1 == 3222 .or. abs(kl1) == 3122 .or. kl1 &
            == -3112 .or. abs(kl1) == 3212 .or. abs(kl1) == 3312 .or.    &
            abs(kl1) == 3322 .or. abs(kl1) == 3334.or.                   &
            kl1 == 2224.or.kl1 == 2114.or.kl1 == 2214) ) iflag = 1
!       Sigma/Lambda/Xi/Omega/Delta + pi-
        if( kl1 == -211 .and. (kl == 3222 .or. abs(kl) == 3122 .or. kl &
            == -3112 .or. abs(kl) == 3212 .or. abs(kl) == 3312 .or.    &
            abs(kl) == 3322 .or. abs(kl) == 3334.or.                   &
            kl == 2224.or.kl == 2114.or.kl == 2214) ) iflag = 1
!       pi0 + Sigma/Lambda/Xi/Omega/Delta
        if( kl == 111 .and. (abs(kl1) == 3112 .or. abs(kl1) == 3122 .or. &
            abs(kl1) == 3222 .or. abs(kl1) == 3212 .or. abs(kl1) == 3312 &
            .or. abs(kl1) == 3322 .or. abs(kl1) == 3334.or.kl1 == 2224   &
            .or.kl1 == 2114.or.kl1 == 2214.or.kl1 == 1114) ) iflag = 1
!       Sigma/Lambda/Xi/Omega/Delta + pi0
        if( kl1 == 111 .and. (abs(kl) == 3112 .or. abs(kl) == 3122 .or. &
            abs(kl) == 3222 .or. abs(kl) == 3212 .or. abs(kl) == 3312   &
            .or. abs(kl) == 3322 .or. abs(kl) == 3334.or.kl == 2224     &
            .or.kl == 2114.or.kl == 2214.or.kl == 1114) ) iflag = 1
!       K+ + K-/Kbar0/Sigma/Lambda/Xi/Omega
        if( kl == 321 .and. (kl1 == -321 .or. kl1 == -311 .or. kl1 ==    &
            3222 .or. kl1 == 3212 .or. kl1 == 3112 .or. kl1 == 3122 .or. &
            kl1 == 3312 .or. kl1 == 3322 .or. kl1 == 3334) ) iflag = 1
!       K-/Kbar0/Sigma/Lambda/Xi/Omega + K+
        if( kl1 == 321 .and. (kl == -321 .or. kl == -311 .or. kl ==   &
            3222 .or. kl == 3212 .or. kl == 3112 .or. kl == 3122 .or. &
            kl == 3312 .or. kl == 3322 .or. kl == 3334) ) iflag = 1
!       K- + K+/K0/Sigma/Lambda/Xi/Omega
        if( kl == -321 .and. (kl1 == 311 .or. kl1 == -3222 .or. kl1 ==  &
            -3212 .or. kl1 == -3112 .or. kl1 == -3122 .or. kl1 == -3312 &
            .or. kl1 == -3322 .or. kl1 == -3334) ) iflag = 1
!       K+/K0/Sigma/Lambda/Xi/Omega + K-
        if( kl1 == -321 .and. (kl == 311 .or. kl == -3222 .or. kl == &
            -3212 .or. kl == -3112 .or. kl == -3122 .or. kl == -3312 &
            .or. kl == -3322 .or. kl == -3334) ) iflag = 1
!       K0 + Kbar0/Sigma/Lambda/Xi/Omega
        if( kl == 311 .and. (kl1 == -311 .or. kl1 == 3222 .or. kl1 == &
            3212 .or. kl1 == 3112 .or. kl1 == 3122 .or.               &
            kl1 == 3312 .or. kl1 == 3322 .or. kl1 == 3334) ) iflag = 1
!       Kbar0/Sigma/Lambda/Xi/Omega + K0
        if( kl1 == 311 .and. (kl == -311 .or. kl == 3222 .or. kl == &
            3212 .or. kl == 3112 .or. kl == 3122 .or.               &
            kl == 3312 .or. kl == 3322 .or. kl == 3334) ) iflag = 1
!       Kbar0 + Sigma/Lambda/Xi/Omega
        if( kl == -311 .and. (kl1 == -3222 .or. kl1 == -3212 .or. kl1 &
            == -3112 .or. kl1 == -3122 .or. kl1 == -3312 .or.         &
            kl1 == -3322 .or. kl1 == -3334) ) iflag = 1
!       Sigma/Lambda/Xi/Omega + Kbar0
        if( kl1 == -311 .and. (kl == -3222 .or. kl == -3212 .or. kl &
            == -3112 .or. kl == -3122 .or. kl == -3312 .or.         &
            kl == -3322 .or. kl == -3334) ) iflag = 1
!       Delta/rho + p
        if( kl1 == 2212.and.(kl == 1114.or.kl == 2114.or. &
            kl == 2214.or.abs(kl) == 213.or.kl == 113) ) iflag = 1
!       Delta/rho + n
        if( kl1 == 2112.and.(kl == 2224.or.kl == 2114.or. &
            kl == 2214.or.abs(kl) == 213.or.kl == 113) ) iflag = 1
!       p + Delta/rho
        if( kl == 2212.and.(kl1 == 1114.or.kl1 == 2114.or. &
            kl1 == 2214.or.abs(kl1) == 213.or.kl1 == 113) ) iflag = 1
!       n + Delta/rho
        if( kl == 2112.and.(kl1 == 2224.or.kl1 == 2114.or. &
            kl1 == 2214.or.abs(kl1) == 213.or.kl1 == 113 ) ) iflag = 1

!       Constraints on the J/psi (psi') induced reactions.
!       J/psi or psi' + p/n/pi/rho
        if( (      kl  == 443  .or. kl  == 100443)               &
            .and. (kl1 == 2212 .or. kl1 == 2112                  &
            .or.   kl1 == 211  .or. kl1 == 111 .or. kl1 == -211  &
            .or.   kl1 == 213  .or. kl1 == 113 .or. kl1 == -213) ) iflag = 1
!       p/n/pi/rho + J/psi or psi'
        if( (      kl1 == 443  .or. kl1 == 100443)              &
            .and. (kl  == 2212 .or. kl  == 2112                 &
            .or.   kl  == 211  .or. kl  == 111 .or. kl == -211  &
            .or.   kl  == 213  .or. kl  == 113 .or. kl == -213) ) iflag = 1

!       Constraints on the D meson induced reactions.
!       D + pion/rho/p/n
        if( (ABS(kl)  == 411  .OR. ABS(kl)  == 421  .OR.  &
             ABS(kl)  == 413  .OR. ABS(kl)  == 423) .AND. &
            (ABS(kl1) == 211  .OR. ABS(kl1) == 111  .OR.  &
             ABS(kl1) == 213  .OR. ABS(kl1) == 113  .OR.  &
             ABS(kl1) == 2212 .OR. ABS(kl1) == 2112) ) iflag = 1
!       pion/rho/p/n + D
        if( (ABS(kl1) == 411  .OR. ABS(kl1) == 421  .OR.  &
             ABS(kl1) == 413  .OR. ABS(kl1) == 423) .AND. &
            (ABS(kl)  == 211  .OR. ABS(kl)  == 111  .OR.  &
             ABS(kl)  == 213  .OR. ABS(kl)  == 113  .OR.  &
             ABS(kl)  == 2212 .OR. ABS(kl)  == 2112) ) iflag = 1


        return
        end



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine tcolij_h( l, l1, time, icp )
!!      Calculates the collision time & fill up lc(i,1-2), tc(i).
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        PARAMETER (KSZJ=80000,NSIZE=750000)
        common/sa1_h/nsa,non1,ksa(kszj,5),psa(kszj,5),vsa(kszj,5)
        common/sa8_h/tau(kszj),ishp(kszj)
        common/sa10_h/csnn,cspin,cskn,cspipi,cspsn,cspsm,rcsit,ifram, &
         iabsb,iabsm,i_sigma_AQM,csspn,csspm
        common/sa19_h/coor(3)
        common/sa20_h/t0,cspipiKK,dep,ddt,edipi,epin,ecsnn,ekn,ecspsn,ecspsm, &
         rnt,rnp,ecsspn,ecsspm
        common/sa24/adj1(40),nnstop,non24,zstop
        common/ctllist_t/ lc(nsize,5), tc(nsize), tw(nsize)
        dimension dr(3),db(3),pi(4),pj(4),vi(3),vj(3)
        dimension ri(4),rj(4),rfi(4),rfj(4),b(3)


        drmax = adj1(28)
        KF1   = ksa(l,2)
        KF2   = ksa(l1,2)
        pel   = psa(l,4)
        pel1  = psa(l1,4)
        if( pel  < 1D-10 ) pel  = 1D-10
        if( pel1 < 1D-10 ) pel1 = 1D-10
        pi(4) = pel
        pj(4) = pel1
        do i=1,3,1
            pi(i) = psa(l,i)
            pj(i) = psa(l1,i)
            b(i)  = ( pi(i) + pj(i) ) / ( pi(4) + pj(4) )
        end do
        ilo = 0
!       Performs Lorentz transf. to CMS frame.
        call lorntz(ilo,b,pi,pj)
        bta = SQRT( b(1)**2 + b(2)**2 + b(3)**2 )
!       If boost is too violent,put particles on mass shell by hand.
        bmi = PYMASS( KF1 )
        bmj = PYMASS( KF2 )
        if( bta > 0.99999D0 )then
            pi(4) = SQRT( bmi**2 + pi(1)**2 + pi(2)**2 + pi(3)**2 )
            pj(4) = SQRT( bmj**2 + pj(1)**2 + pj(2)**2 + pj(3)**2 )
        end if
        ss = pi(4) + pj(4)
        ! P. Yang
        ecut0 = 0.02D0
        ecut  = ss - bmi - bmj
        if( ecut <= ecut0 ) return
!       if ss < threshold collision may not happen

!       For isothermal direct reactions
!       pipi->KK
        if( (     ( abs(KF1) == 211 .or. KF1 == 111 )   &
            .and. ( abs(KF2) == 211 .or. KF2 == 111 ) ) &
            .and. ss <= 2.*0.5 ) return
!       piN->KY
        if( ( (     abs(KF1) == 211  .or.     KF1  == 111 )    &
            .and. ( abs(KF2) == 2112 .or. abs(KF2) == 2212 ) ) &
            .and. ss <= ( 0.5 + 1.2 ) ) return
        if( ( (     abs(KF2) == 211  .or.     KF2 == 111 )     &
            .and. ( abs(KF1) == 2112 .or. abs(KF1) == 2212 ) ) &
            .and. ss <= ( 0.5 + 1.2 ) ) return
!       piY->K(Cascade)
        if( ( (     abs(KF1) == 211  .or.     KF1  == 111 )    &
            .and. ( abs(KF2) == 3122 .or. abs(KF2) == 3212     &
            .or.        KF2  == 3112 .or.     KF2  == 3222 ) ) &
            .and. ss <= (0.5 + 1.322) ) return
        if( ( (     abs(KF2) == 211  .or.     KF2  == 111 )    &
            .and. ( abs(KF1) == 3122 .or. abs(KF1) == 3212     &
            .or.        KF1  == 3112 .or.     KF1  == 3222 ) ) &
            .and. ss <= (0.5 + 1.322) ) return
!       pi(Cascade)->K(Omega)
        if( ( (     abs(KF1) == 211  .or.     KF1  == 111 )    &
            .and. ( abs(KF2) == 3312 .or. abs(KF2) == 3322 ) ) &
            .and. ss <= (0.5 + 1.7) ) return
        if( ( (     abs(KF2) == 211  .or.     KF2  == 111 )    &
            .and. ( abs(KF1) == 3312 .or. abs(KF1) == 3322 ) ) &
            .and. ss <= (0.5 + 1.7) ) return
!       KN->K(Cascade)
        if( ( (     abs(KF1) == 321  .or. abs(KF1) == 311 )    &
            .and. ( abs(KF2) == 2112 .or. abs(KF2) == 2212 ) ) &
            .and. ss <= (0.5 + 1.322) ) return
        if( ( (     abs(KF2) == 321  .or. abs(KF2) == 311)     &
            .and. ( abs(KF1) == 2112 .or. abs(KF1) == 2212 ) ) &
            .and. ss <= (0.5 + 1.322) ) return
!       piN->pi(Delta)
        if( ( (     abs(KF1) == 211  .or. KF1 == 111 )         &
            .and. ( abs(KF2) == 2112 .or. abs(KF2) == 2212 ) ) &
            .and. ss <= (0.14 + 1.234) ) return
        if( ( (     abs(KF2) == 211  .or. KF2 == 111 )         &
            .and. ( abs(KF1) == 2112 .or. abs(KF1) == 2212 ) ) &
            .and. ss <= (0.14 + 1.234) ) return
!       piN->(rho)N
        if( ( (     abs(KF1) == 211  .or. KF1 == 111 )         &
            .and. ( abs(KF2) == 2112 .or. abs(KF2) == 2212 ) ) &
            .and. ss <= 1.71 ) return
        if( ( (     abs(KF2) == 211  .or. KF2 == 111 )         &
            .and. ( abs(KF1) == 2112 .or. abs(KF1) == 2212 ) ) &
            .and. ss <= 1.71 ) return
!       NN->N(Delta)
        if( ( (     abs(KF1) == 2112  .or. abs(KF1) == 2212 )  &
            .and. ( abs(KF2) == 2112 .or. abs(KF2) == 2212 ) ) &
            .and. ss <= 2.172 ) return
!       K(Cascade)->pi(Omega)
        if( ( (     abs(KF1) == 321  .or. abs(KF1) == 311 )    &
            .and. ( abs(KF2) == 3312 .or. abs(KF2) == 3322 ) ) &
            .and. ss <= 1.84 ) return
        if( ( (     abs(KF2) == 321  .or. abs(KF2) == 311 )    &
            .and. ( abs(KF1) == 3312 .or. abs(KF1) == 3322 ) ) &
            .and. ss <= 1.84 ) return

!       for isothermal reverse reactions
!       piY->KN
        if( ( (     abs(KF1) == 211  .or.     KF1  == 111 )    &
            .and. ( abs(KF2) == 3122 .or. abs(KF2) == 3212     &
            .or.        KF2  == 3112 .or.     KF2  == 3222 ) ) &
            .and. ss <= 1.44 ) return
        if( ( (     abs(KF2) == 211  .or.     KF2  == 111 )    &
            .and. ( abs(KF1) == 3122 .or. abs(KF1) == 3212     &
            .or.        KF1  == 3112 .or.      KF1 == 3222 ) ) &
            .and. ss <= 1.44 ) return
!       pi(Cascade)->KY
        if( ( (     abs(KF1) == 211  .or.     KF1  == 111 )    &
            .and. ( abs(KF2) == 3312 .or. abs(KF2) == 3322 ) ) &
            .and. ss <= 1.7 ) return
        if( ( (     abs(KF2) == 211  .or.     KF2  == 111 )    &
            .and. ( abs(KF1) == 3312 .or. abs(KF1) == 3322 ) ) &
            .and. ss <= 1.7 ) return

        do i=1,4,1
            ri(i) = vsa(l,i)
            rj(i) = vsa(l1,i)
        end do
!       ri(4) = time
!       rj(4) = time
!       Lorentz transf. to CMS frame
        call lorntz(ilo,b,ri,rj)
        rb = 0D0
        bb = 0D0
        rr = 0D0
        rtai = 0D0
        kflag = 0
        do ik=1,3,1
            vi(ik) = pi(ik) / pi(4)
            vj(ik) = pj(ik) / pj(4)
        enddo

        do i=1,3
            rfi(i) = vsa(l,i)  + ( tau(l)  - time ) * ( psa(l,i)  / psa(l,4)  )
            rfj(i) = vsa(l1,i) + ( tau(l1) - time ) * ( psa(l1,i) / psa(l1,4) )
        end do
        rfi(4) = tau(l)
        rfj(4) = tau(l1)
        call lorntz(ilo,b,rfi,rfj)
!       gamli = psa(l,4)  / psa(l,5)
!       gamlj = psa(l1,4) / psa(l1,5)
        ctaui = rfi(4)
        ctauj = rfj(4)
        tcol  = ctaui
        if( ctaui < ctauj ) tcol = ctauj
        do ik=1,3,1
            db(ik) = ( vi(ik) - vj(ik) ) * tcol
            dr(ik) = ri(ik) - rj(ik) - ( vi(ik)*ri(4) - vj(ik)*rj(4) ) + db(ik)
            rtai = rtai + dr(ik)*dr(ik)
        end do
        dot = 0D0
        do ik=1,3,1
            dot = dr(ik) * pi(ik) + dot
        end do
!       dot = -1
        if( dot >= 0D0 )then
            kflag = 1
            if( tcol <= ri(4) ) return
            if( tcol <= rj(4) ) return
        else
            rtai = 0D0
            do ik=1,3,1
                dr(ik) = ri(ik) - rj(ik) - ( vi(ik)*ri(4) - vj(ik)*rj(4) )
                db(ik) = vi(ik) - vj(ik)
                rb = rb + dr(ik)*db(ik)
                bb = bb + db(ik)*db(ik)
                rr = rr + dr(ik)*dr(ik)
            end do
            if( bb  <=  1D-10 ) return
            tcol = 0D0 - rb/bb
            if( tcol <= ri(4) ) return
            if( tcol <= rj(4) ) return
            if( tcol-ctaui  <=  0D0 ) return
            if( tcol-ctauj  <=  0D0 ) return
!       For collision occurs, time must one step ahead.
!Tai
            do ik=1,3,1
                dr(ik) = ri(ik)-rj(ik)-(vi(ik)*ri(4)-vj(ik)*rj(4)) + tcol*db(ik)
                rtai = rtai + dr(ik)*dr(ik)
            end do
!           gamai = pi(4) / PYMASS( KF1 )
!           gamaj = pj(4) / PYMASS( KF2 )
!TAIAN
!       When collision happens, particles should already be produced.
!       We give a zero formation time for particles produced from rescatttering.
!Tai
        end if
        sg = rtai

        dmin = sqrt(sg)
        call intdis_h( l, l1, rsig )
!       'intdis': calculate the interaction distance between particles l & l1.
        if( dmin > rsig ) return
!       Distance between the two particles should be smaller than rsig.
        do ik=1,3,1
            ri(ik) = ri(ik) + vi(ik) * ( tcol - ri(4) )
            rj(ik) = rj(ik) + vj(ik) * ( tcol - rj(4) )
        end do
!       Moves along Newton trajectory in CMS.
        ri(4) = tcol
        rj(4) = tcol
        ilo = 1
!       Transforms back to Lab.
        call lorntz(ilo,b,ri,rj)
        tcol1 = ri(4)
        tcol2 = rj(4)
        if( kflag == 0 )then
            if( tcol1 - tau(l)  < 0D0 ) return
            if( tcol2 - tau(l1) < 0D0 ) return
        else
            if( tcol1 - tau(l)  < -1D-4 ) return
            if( tcol2 - tau(l1) < -1D-4 ) return
        end if
        if( ri(4) > rj(4) ) ri(4) = rj(4)
        tcol = ri(4)
        if( tcol <= time ) return
!       Collision happens in the future.
        do i=1,3,1
            ri(i) = vsa(l,i)  + psa(l,i)  * ( tcol - time ) / pel  - coor(i)
            rj(i) = vsa(l1,i) + psa(l1,i) * ( tcol - time ) / pel1 - coor(i)
        end do
        rri = sqrt( ri(1) * ri(1) + ri(2) * ri(2) + ri(3) * ri(3) )
        rrj = sqrt( rj(1) * rj(1) + rj(2) * rj(2) + rj(3) * rj(3) )
!       Particles must be inside the largest region considered.
        if( rri > drmax ) return
        if( rrj > drmax ) return
        if( tcol <= drmax )then
            tc(icp)   = tcol
            lc(icp,1) = l
            lc(icp,2) = l1
        end if


        return
        end



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine intdis_h(l,l1,rsig)
!!      Calculates the interaction distance between particles l.
!!       and l1.
!!      It plays also the role of second range filter.
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        PARAMETER (KSZJ=80000)
        common/sa1_h/nsa,non1,ksa(kszj,5),psa(kszj,5),vsa(kszj,5)
        common/sa20_h/t0,cspipiKK,dep,ddt,edipi,epin,ecsnn,ekn,ecspsn,ecspsm, &
         rnt,rnp,ecsspn,ecsspm


        rsig  = 0D0
        kl    = ksa(l,2)
        kl1   = ksa(l1,2)
        idpl  = 0
        idpl1 = 0
        pio   = 3.141592653589793D0

!       rsig in fm, sigma_tot in fm^2.
!         rsig = SQRT( sigma_tot / pi )
!       ecsnn in fm, csnn in fm^2.
!        ecsnn = SQRT( csnn / pi ) (i.e. rsig for NN)
!       csnn in fm^2, sigma_NN in mb. Note 1 mb = 0.1 fm^2.
!        csnn = sigma_NN * 0.1, sigma_NN in mb.
!       UseS AQM cross sections. sigma_AQM is a function that returns in mb.
        sigma_NN = ecsnn**2 * pio * 10D0
        rsig = SQRT( sigma_AQM( kl, kl1, sigma_NN )*0.1D0 / pio )

        ! N
        if( abs(kl) == 2212 .or. abs(kl) == 2112 )      idpl = 1
        ! J/psi & psi'
        if( abs(kl) == 443  .or. abs(kl) == 100443 )    idpl = 2
        ! pi
        if( abs(kl) == 211  .or. abs(kl) == 111 )       idpl = 3
        ! K
        if( abs(kl) == 321  .or. abs(kl) == 311 )       idpl = 4
        ! Sigma, Lambda, Xi, Omega
        if(      abs(kl) == 3212 .or. abs(kl) == 3112 .or. abs(kl) == 3222 &
            .or. abs(kl) == 3122 .or. abs(kl) == 3312 .or. abs(kl) == 3322 &
            .or. abs(kl) == 3334 )                      idpl = 5
        ! rho
        if( abs(kl) == 213  .or. abs(kl) == 113 )       idpl = 6
        ! Delta
        if( kl == 1114 .or. kl == 2114 .or. kl == 2214 .or. kl == 2224 ) &
                                                        idpl = 7

        if( abs(kl1) == 2212 .or. abs(kl1) == 2112 )    idpl1 = 1
        if( abs(kl1) == 443  .or. abs(kl1) == 100443 )  idpl1 = 2
        if( abs(kl1) == 211  .or. abs(kl1) == 111 )     idpl1 = 3
        if( abs(kl1) == 321  .or. abs(kl1) == 311 )     idpl1 = 4
        if(      abs(kl1) == 3212 .or. abs(kl1) == 3112 .or. abs(kl1) == 3222 &
            .or. abs(kl1) == 3122 .or. abs(kl1) == 3312 .or. abs(kl1) == 3322 &
            .or. abs(kl1) == 3334 )                     idpl1 = 5
        if( abs(kl1) == 213  .or. abs(kl1) == 113 )     idpl1 = 6
        if( kl1 == 1114 .or. kl1 == 2114 .or. kl1 == 2214 .or. kl1 == 2224 ) &
                                                        idpl1 = 7

        if( idpl == 1 .and. idpl1 == 1 ) rsig = ecsnn

        if( idpl == 2 .and. idpl1 == 1 )then
            rsig = ecspsn
            if( kl == 100443 ) rsig = ecsspn
        end if
        if( idpl == 1 .and. idpl1 == 2)then
            rsig = ecspsn
            if( kl1 == 100443 ) rsig = ecsspn
        end if
        if( idpl == 2 .and. idpl1 == 3 )then
            rsig = ecspsm
            if( kl == 100443 ) rsig = ecsspm
        end if
        if( idpl == 3 .and. idpl1 == 2 )then
            rsig = ecspsm
            if( kl1 == 100443 ) rsig = ecsspm
        end if

        if( idpl == 2 .and. idpl1 == 6 )then
            rsig = ecspsm
            if( kl1 == 113 ) rsig = ecspsm * SQRT(2D0)
!       At the case of rho0, cross section enlarges a factor 2 to
!        consider the effect of omega.
            if( kl == 100443)then
                rsig = ecsspm
                if( kl1 == 113 ) rsig = ecsspm * SQRT(2D0)
            end if
        end if
        if( idpl == 6 .and. idpl1 == 2)then
            rsig = ecspsm
            if( kl == 113) rsig = ecspsm * SQRT(2D0)
            if( kl1 == 100443 )then
                rsig = ecsspm
                if( kl == 113 ) rsig = ecsspm * SQRT(2D0)
            end if
        end if

        if( idpl == 3 .and. idpl1 == 3 ) rsig = edipi
        if( idpl == 1 .and. idpl1 == 3 ) rsig = epin
        if( idpl == 3 .and. idpl1 == 1 ) rsig = epin
        if( idpl == 3 .and. idpl1 == 5 ) rsig = epin
        if( idpl == 5 .and. idpl1 == 3 ) rsig = epin
!       Assume the total cross section of (pion)Y,((pion)Cascade) and
!        ((pion)Omega) = (pion)N .
        if(  idpl == 4 .and. (idpl1 == 1  .or.  idpl1 == 5) ) rsig = ekn
        if( (idpl == 1 .or.   idpl  == 5) .and. idpl1 == 4 )  rsig = ekn
        if(  idpl == 4 .and.  idpl1 == 4 ) rsig = edipi

!       Assume the total cross section of KY (K Cascade) and (K Omega) = KN .
!       Assume the total cross section of KK=(pion)(pion) .
        if( idpl  == 1 .and. idpl1 == 6 ) rsig = epin
        if( idpl  == 1 .and. idpl1 == 7 ) rsig = ecsnn
        if( idpl  == 3 .and. idpl1 == 7 ) rsig = epin
        if( idpl1 == 1 .and. idpl  == 6 ) rsig = epin
        if( idpl1 == 1 .and. idpl  == 7 ) rsig = ecsnn
        if( idpl1 == 3 .and. idpl  == 7 ) rsig = epin
        if( idpl  == 1 .and. idpl1 == 5 ) rsig = ecsnn
        if( idpl  == 5 .and. idpl1 == 1 ) rsig = ecsnn
!       Assume the total cross section of (rho)N=(pion)N .
!       Assume the total cross section of N(Delta)=NN .
!       Assume the total cross section of (pion)(Delta)=(pion)N .
!       Assume the total cross section of NY, N(Cascade), and N(Omega)=NN .


        return
        end



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        function s0715(ss,ilo,i,the)
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        common/sa13/kjp20,non13,vjp20,vjp21,vjp22,vjp23
!       pion- + p to k+ + sigma-,channel 4
        ilo=0
        si=0.
        ii=i
        if(ii == 1) goto 30
        if(ss <= the) goto 100
30      ilo=1
        if(ss >= 1.9) goto 10
        IF(kjp20 == 0)then
        si=0.25*(1.-0.75*(ss-1.691))
        else
        si=vjp20
        endif
        goto 100
10      IF(kjp20 == 0)then
        si=309.1*exp(max(-40.,-3.77*ss))
        else
        si=vjp20
        endif
100     continue
        s0715=si
        return
        end



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        function s07122(ss,ilo,i,the)
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        common/sa13/kjp20,non13,vjp20,vjp21,vjp22,vjp23
!       pion- + p to k0 + lambda,channel 5
        ilo=0
        si=0.
        ii=i
        if(ii == 1) goto 30
        if(ss <= the) goto 100
30      ilo=1
        if(ss >= 1.684) goto 10
        IF(kjp20 == 0)THEN
        si=0.9/0.091*(ss-1.613)
        ELSE
        si=vjp20
        ENDIF
        goto 100
10      if(ss >= 2.1) goto 20
        IF(kjp20 == 0)THEN
        si=436.3*exp(max(-40.,-4.154*ss))
        ELSE
        si=vjp20
        ENDIF
        goto 100
20      IF(kjp20 == 0)THEN
        si=0.314*exp(max(-40.,-0.301*ss))
        ELSE
        si=vjp20
        ENDIF
100     s07122=si
        return
        end



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        function s07123(ss,ilo,i,the)
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        common/sa13/kjp20,non13,vjp20,vjp21,vjp22,vjp23
!       pion- + p to k0 + sigma0,channel 6
        ilo=0
        si=0.
        ii=i
        if(ii == 1) goto 30
        if(ss <= the) goto 100
30      ilo=1
        if(ss >= 1.722) goto 10
        IF(kjp20 == 0)THEN
        si=10.6*(ss-1.689)
        ELSE
        si=vjp20
        ENDIF
        goto 100
10      if(ss >= 3.) goto 20
        IF(kjp20 == 0)THEN
        si=13.7*exp(max(-40.,-1.92*ss))
        ELSE
        si=vjp20
        ENDIF
        goto 100
20      IF(kjp20 == 0)THEN
        si=0.188*exp(max(-40.,-0.611*ss))
        ELSE
        si=vjp20
        ENDIF
100     s07123=si
        return
        end



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        function s1724(ss,ilo,i,the)
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        common/sa13/kjp20,non13,vjp20,vjp21,vjp22,vjp23
!       pion+ + n to k+ + sigma0,channel 8
        ilo=0
        si=0.
        ii=i
        if(ii == 1) goto 30
        if(ss <= the) goto 100
30      ilo=1
        IF(kjp20 == 0)THEN
        si=0.25*(s0715(ss,ilo,ii,the)+s07123(ss,ilo,ii,the)+ &
         s1713(ss,ilo,ii,the))
        ELSE
        si=vjp20
        ENDIF
100     s1724=si
        return
        end



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        function s1727(ss,ilo,i,the)
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        common/sa13/kjp20,non13,vjp20,vjp21,vjp22,vjp23
!       pion+ + n to k+ + lambda,channel 9
        ilo=0
        si=0.
        ii=i
        if(ii == 1) goto 30
        if(ss <= the) goto 100
30      ilo=1
        if(ss >= 1.684) goto 10
        IF(kjp20 == 0)THEN
        si=0.9*(ss-1.613)/0.091
        ELSE
        si=vjp20
        ENDIF
        goto 100
10      if(ss >= 2.1) goto 20
        IF(kjp20 == 0)THEN
        si=436.3*exp(max(-40.,-4.154*ss))
        ELSE
        si=vjp20
        ENDIF
        goto 100
20      IF(kjp20 == 0)THEN
        si=0.314*exp(max(-40.,-0.301*ss))
        ELSE
        si=vjp20
        ENDIF
!*SA
!100    s1727=si*0.25
100     s1727=si
        return
        end



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        function s1713(ss,ilo,i,the)
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        common/sa13/kjp20,non13,vjp20,vjp21,vjp22,vjp23
!       pion+ + p to k+ + sigma+,channel 7
        ilo=0
        si=0.
        ii=i
        if(ii == 1) goto 30
        if(ss <= the) goto 100
30      ilo=1
        if(ss >= 1.934) goto 10
        IF(kjp20 == 0)THEN
        si=0.7*(ss-1.683)/0.218
        ELSE
        si=vjp20
        ENDIF
        goto 100
10      if(ss >= 3.) goto 20
        IF(kjp20 == 0)THEN
        si=60.26*exp(max(-40.,-2.31*ss))
        ELSE
        si=vjp20
        ENDIF
        goto 100
20      IF(kjp20 == 0)THEN
        si=0.36*exp(max(-40.,-0.605*ss))
        ELSE
        si=vjp20
        ENDIF
100     s1713=si
        return
        end



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        function s2325(ss,ilo,i,the)
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
!       pion0 + n to k+ + sigma-,channel 12
        ii=i
        s2325=s1724(ss,ilo,ii,the)
        return
        end



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        function s2314(ss,ilo,i,the)
!       pion0 + p to k+ + sigma0,channel 10
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        ii=i
        s2314=s1724(ss,ilo,ii,the)
        return
        end



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        function s2317(ss,ilo,i,the)
!       pion0 + p to k+ + lambda,channel 11
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        ii=i
        s2317=s1727(ss,ilo,ii,the)
        return
        end



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine ppdelta(ss,icp,ioo)
!       a part of 'prod' to deal with pp -> ...
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        INTEGER PYCOMP
        PARAMETER(NSIZE=750000)
!       calculate particle production weight and fill up lc(i,3-5),tw(i).
!       tw : the ratio of cross section of (special inela.)/tot
        COMMON/PYDAT2/KCHG(500,4),PMAS(500,4),PARF(2000),VCKM(4,4)
        COMMON/PYCIDAT2/KFMAXT,NONT2,PARAM(20),WEIGH(600)
        SAVE /PYCIDAT2/
        common/sa10_h/csnn,cspin,cskn,cspipi,cspsn,cspsm,rcsit,ifram, &
         iabsb,iabsm,i_sigma_AQM,csspn,csspm
        common/count_h/isinel(600)
        common/iloval/ilo1,ilo2,ilo3,ilo4,ilo5,ilo6,ilo7,ilo8,ilo9
        common/ctllist_t/ lc(nsize,5), tc(nsize), tw(nsize)


        if(isinel(174) == 0)then
        sigma1=0.
        goto 522
        endif
        the=pmas(pycomp(2224),1)+pmas(pycomp(2112),1)
        sigma1=snn(ss,ilo1,0,the)*WEIGH(46)
!       cross section of p + p to delta++ + n

522     if(isinel(173) == 0)then
        sigma2=0.
        goto 523
        endif
        the=pmas(pycomp(2214),1)+pmas(pycomp(2212),1)
        sigma2=snn(ss,ilo2,0,the)*WEIGH(45)
!       cross section of p+p to delta+ +p
523     if(ilo1 == 0.and.ilo2 == 0) goto 13
        if(sigma1 < 1.e-6 .and. sigma2 < 1.e-6) goto 13
        ik1=2224
        ik2=2112
        ic=174
!       p+p to delta++ +  n
        sigm12=sigma1+sigma2
        if(pyr(1) > sigma1/sigm12)then
        ik1=2214
        ik2=2212
        ic=173
!       cross section of p+p to delta+ +p
        endif
        lc(icp,3)=ik1
        lc(icp,4)=ik2
        lc(icp,5)=ic
        tw(icp)=sigm12/csnn/10.
        ioo=1


13      return
        end



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine pndelta(ss,icp,ioo)
!       a part of 'prod' to deal with pn -> ...
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        INTEGER PYCOMP
        PARAMETER(NSIZE=750000)
!       calculate particle production weight and fill up lc(i,3-5),tw(i).
!       tw : the ratio of cross section of (special inela.)/tot
        COMMON/PYDAT2/KCHG(500,4),PMAS(500,4),PARF(2000),VCKM(4,4)
        COMMON/PYCIDAT2/KFMAXT,NONT2,PARAM(20),WEIGH(600)
        SAVE /PYCIDAT2/
        common/sa10_h/csnn,cspin,cskn,cspipi,cspsn,cspsm,rcsit,ifram, &
         iabsb,iabsm,i_sigma_AQM,csspn,csspm
        common/count_h/isinel(600)
        common/iloval/ilo1,ilo2,ilo3,ilo4,ilo5,ilo6,ilo7,ilo8,ilo9
        common/ctllist_t/ lc(nsize,5), tc(nsize), tw(nsize)


        if(isinel(175) == 0)then
        sigma1=0.
        goto 524
        endif
        the=pmas(pycomp(2214),1)+pmas(pycomp(2112),1)
        sigma1=snn(ss,ilo1,0,the)*WEIGH(47)
!       cross section of p + n to delta+ + n

524     if(isinel(176) == 0)then
        sigma2=0.
        goto 525
        endif
        the=pmas(pycomp(2114),1)+pmas(pycomp(2212),1)
        sigma2=snn(ss,ilo2,0,the)*WEIGH(48)
!       cross section of p+n to delta0 +p
525     if(ilo1 == 0.and.ilo2 == 0) goto 13
        if(sigma1 < 1.e-6 .and. sigma2 < 1.e-6) goto 13
        ik1=2214
        ik2=2112
        ic=175
!       p+n to delta+ +  n
        sigm12=sigma1+sigma2
        if(pyr(1) > sigma1/sigm12)then
        ik1=2114
        ik2=2212
        ic=176
!       cross section of p+n to delta0 +p
        endif
        lc(icp,3)=ik1
        lc(icp,4)=ik2
        lc(icp,5)=ic
        tw(icp)=sigm12/csnn/10.
        ioo=1


13      return
        end



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine nndelta(ss,icp,ioo)
!       a part of 'prod' to deal with nn-> ...
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        INTEGER PYCOMP
        PARAMETER(NSIZE=750000)
!       calculate particle production weight and fill up lc(i,3-5),tw(i).
!       tw : the ratio of cross section of (special inela.)/tot
        COMMON/PYDAT2/KCHG(500,4),PMAS(500,4),PARF(2000),VCKM(4,4)
        COMMON/PYCIDAT2/KFMAXT,NONT2,PARAM(20),WEIGH(600)
        SAVE /PYCIDAT2/
        common/sa10_h/csnn,cspin,cskn,cspipi,cspsn,cspsm,rcsit,ifram, &
         iabsb,iabsm,i_sigma_AQM,csspn,csspm
        common/count_h/isinel(600)
        common/iloval/ilo1,ilo2,ilo3,ilo4,ilo5,ilo6,ilo7,ilo8,ilo9
        common/ctllist_t/ lc(nsize,5), tc(nsize), tw(nsize)


        if(isinel(177) == 0)then
        sigma1=0.
        goto 527
        endif
        the=pmas(pycomp(2114),1)+pmas(pycomp(2112),1)
        sigma1=snn(ss,ilo1,0,the)*WEIGH(49)
!       cross section of n + n to delta0 + n

527     if(isinel(178) == 0)then
        sigma2=0.
        goto 526

        endif
        the=pmas(pycomp(1114),1)+pmas(pycomp(2212),1)
        sigma2=snn(ss,ilo2,0,the)*WEIGH(50)
!       cross section of n+n to delta- + p
526     if(ilo1 == 0.and.ilo2 == 0) goto 13
        if(sigma1 < 1.e-6 .and. sigma2 < 1.e-6) goto 13
        ik1=2114
        ik2=2112
        ic=177
!       n+n to delta0 + n
        sigm12=sigma1+sigma2
        if(pyr(1) > sigma1/sigm12)then
        ik1=1114
        ik2=2212
        ic=178
!       cross section of n+n to delta- + p
        endif
        lc(icp,3)=ik1
        lc(icp,4)=ik2
        lc(icp,5)=ic
        tw(icp)=sigm12/csnn/10.
        ioo=1


13      return
        end



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine pip1(ss,icp,ioo)
!       a part of 'prod' to deal with pion- + p -> ...
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        INTEGER PYCOMP
        PARAMETER(NSIZE=750000)
!       calculate particle production weight and fill up lc(i,3-5),tw(i).
!       tw : the ratio of cross section of (special inela.)/tot
        COMMON/PYDAT2/KCHG(500,4),PMAS(500,4),PARF(2000),VCKM(4,4)
        COMMON/PYCIDAT2/KFMAXT,NONT2,PARAM(20),WEIGH(600)
        common/sa10_h/csnn,cspin,cskn,cspipi,cspsn,cspsm,rcsit,ifram, &
         iabsb,iabsm,i_sigma_AQM,csspn,csspm
        common/count_h/isinel(600)
        common/iloval/ilo1,ilo2,ilo3,ilo4,ilo5,ilo6,ilo7,ilo8,ilo9
        common/ctllist_t/ lc(nsize,5), tc(nsize), tw(nsize)


        if(isinel(11) == 0)then
        sigma1=0.
        goto 401
        endif
        the=pmas(pycomp(321),1)+pmas(pycomp(3112),1)
        sigma1=s0715(ss,ilo1,0,the)*WEIGH(15)
!       cross section of pion- + p to k+ + sigma-
!       cross section is here in the unit of mb
401     if(isinel(12) == 0)then
        sigma2=0.
        goto 402
        endif
        the=pmas(pycomp(311),1)+pmas(pycomp(3122),1)
        sigma2=s07122(ss,ilo2,0,the)*WEIGH(7)
!       cross section of pion- + p to k0 + lambda
402     if(isinel(13) == 0)then
        sigma3=0.
        goto 403
        endif
        the=pmas(pycomp(311),1)+pmas(pycomp(3212),1)

        sigma3=s07123(ss,ilo3,0,the)*WEIGH(16)
!       cross section of pion- + p to k0 + sigma0
403     if(isinel(147) == 0)then
        sigma4=0.
        goto 203
        endif
!       the: threshold energy of a reaction
        the=pmas(pycomp(1114),1)+pmas(pycomp(211),1)
        sigma4=sdelta(ss,ilo4,0,the)*WEIGH(57)
!       cross section of pion- + p to delta- + pi+
203     if(isinel(148) == 0)then
        sigma5=0.
        goto 204
        endif
        the=pmas(pycomp(113),1)+pmas(pycomp(2112),1)
        sigma5=srho(ss,ilo5,0,the)*WEIGH(63)
!       cross section of pion- + p to rho0 + n
204     if(isinel(149) == 0)then
        sigma6=0.
        goto 205
        endif
        the=pmas(pycomp(-213),1)+pmas(pycomp(2212),1)
        sigma6=srho(ss,ilo6,0,the)*WEIGH(64)
!       cross section of pion- + p to rho- + p
205     if(isinel(150) == 0)then
        sigma7=0.
        goto 310
        endif
        the=pmas(pycomp(-211),1)+pmas(pycomp(2214),1)
        sigma7=sdelta(ss,ilo7,0,the)*WEIGH(58)
!       cross section of pion- + p to delta+ + pion-
310     if(isinel(151) == 0)then
        sigma8=0.
        goto 805
        endif
        the=pmas(pycomp(111),1)+pmas(pycomp(2114),1)
        sigma8=sdelta(ss,ilo8,0,the)*WEIGH(58)
!       cross section of pion- + p to delta0 + pion0

805     if(ilo1 == 0 .and. ilo2 == 0 .and. ilo3 == 0 &
         .and.ilo4 == 0 .and. ilo5 == 0 .and. ilo6 == 0 &
         .and.ilo7 == 0.and.ilo8 == 0) goto 13
        if(sigma1 < 1.e-6 .and. sigma2 < 1.e-6 .and. sigma3 < 1.e-6 &
               .and.sigma4 < 1.e-6 .and. sigma5 < 1.e-6 .and. &
                sigma6 < 1.e-6.and. sigma7 < 1.e-6 &
               .and.sigma8 < 1.e-6) goto 13
        sigma12=sigma1+sigma2
        sigma13=sigma12+sigma3
        sigma14=sigma13+sigma4
        sigma15=sigma14+sigma5
        sigma16=sigma15+sigma6
        sigma17=sigma16+sigma7
        sigma18=sigma17+sigma8
        s1=sigma1/sigma18
        s2=sigma12/sigma18
        s3=sigma13/sigma18
        s4=sigma14/sigma18
        s5=sigma15/sigma18
        s6=sigma16/sigma18
        s7=sigma17/sigma18
        rlus=pyr(1)
        if(rlus <= s1)then
        ik1=321
        ik2=3112
!       3112 is the flavor code of sigma-
        ic=11
!       pion- + p to k+ + sigma-
        goto 416
        endif
        if(rlus > s1 .and. rlus <= s2)then
        ik1=311
        ik2=3122
!       3122 is the Kf code of lambda
        ic=12
!       pion- + p to k0 + lambda
        goto 416
        endif
        if(rlus > s2 .and. rlus <= s3)then
        ik1=311
        ik2=3212
!       3212 is the KF code of sigma0
        ic=13
!       pion- + p to k0 + sigma0
        goto 416
        endif
        if(rlus > s3 .and. rlus <= s4)then
        ik1=211
        ik2=1114
!       1114 is the Kf code of delta-
        ic=147
!       pion- + p to pi+ + delta-
        goto 416
        endif
        if(rlus > s4 .and. rlus <= s5)then
        ik1=113
        ik2=2112
!       113 is the Kf code of rho0
        ic=148
!       pion- + p to rho0 + n
        goto 416
        endif
        if(rlus > s5 .and. rlus <= s6)then
        ik1=-213
        ik2=2212
!       -213 is the Kf code of rho-
        ic=149
!       pion- + p to rho- + p
        goto 416
        endif
        if(rlus > s6 .and. rlus <= s7)then
        ik1=-211
        ik2=2214
!       2214 is the Kf code of delta+
        ic=150
!       pion- + p to delta++pion-
        goto 416
        endif
        ik1=111
        ik2=2114
!       2114 is the Kf code of delta0
        ic=151
!       pion- + p to delta0+pion0
        goto 416
416     lc(icp,3)=ik1
        lc(icp,4)=ik2
        lc(icp,5)=ic
        tw(icp)=sigma18/cspin/10.
        ioo=1


13      return
        end



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine pin3(ss,icp,ioo)
!       a part of 'prod' to deat with pion- + n -> ...
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        INTEGER PYCOMP
        PARAMETER(NSIZE=750000)
!       calculate particle production weight and fill up lc(i,3-5),tw(i).
!       tw : the ratio of cross section of (special inela.)/tot
        COMMON/PYDAT2/KCHG(500,4),PMAS(500,4),PARF(2000),VCKM(4,4)
        COMMON/PYCIDAT2/KFMAXT,NONT2,PARAM(20),WEIGH(600)
        SAVE /PYCIDAT2/
        common/sa10_h/csnn,cspin,cskn,cspipi,cspsn,cspsm,rcsit,ifram, &
         iabsb,iabsm,i_sigma_AQM,csspn,csspm
        common/count_h/isinel(600)
        common/iloval/ilo1,ilo2,ilo3,ilo4,ilo5,ilo6,ilo7,ilo8,ilo9
        common/ctllist_t/ lc(nsize,5), tc(nsize), tw(nsize)


        if(isinel(152) == 0)then
        sigma1=0.
        goto 222
        endif
        the=pmas(pycomp(111),1)+pmas(pycomp(1114),1)
        sigma1=sdelta(ss,ilo1,0,the)*WEIGH(61)
!       cross section of pion- + n to delta- + pi0

222     if(isinel(153) == 0)then
        sigma2=0.
        goto 223
        endif
        the=pmas(pycomp(-213),1)+pmas(pycomp(2112),1)
        sigma2=srho(ss,ilo2,0,the)*WEIGH(65)
!       cross section of pion- + n to rho- + n
223     if(isinel(154) == 0)then
        sigma3=0.
        goto 228
        endif
        the=pmas(pycomp(-211),1)+pmas(pycomp(2114),1)
        sigma3=sdelta(ss,ilo3,0,the)*WEIGH(62)
!       cross section of pion- + n to delta0 + pion-
228     if(isinel(14) == 0)then
        sigma4=0.
        goto 818
        endif
        ilo4=1
        the=pmas(pycomp(311),1)+pmas(pycomp(3112),1)

        sigma4=s1724(ss,ilo4,0,the)*WEIGH(24)

!       cross section of pion- + n to k0 + sigma-
818     if(ilo1 == 0.and.ilo2 == 0.and.ilo3 == 0 &
         .and.ilo4 == 0) goto 13
        if(sigma1 < 1.e-6 .and. sigma2 < 1.e-6.and. &
         sigma3 < 1.e-6.and. sigma4 < 1.e-6) goto 13
        sigm12=sigma1+sigma2
        sigm13=sigm12+sigma3
        sigm14=sigm13+sigma4
        s1=sigma1/sigm14
        s2=sigm12/sigm14
        s3=sigm13/sigm14
        rlu1=pyr(1)
        if(rlu1 <= s1)then
        ik1=111
        ik2=1114
        ic=152
!       pion- + n to delta-  +  pi0
        elseif(rlu1 > s1 .and. rlu1 <= s2)then
        ik1=-213
        ik2=2112
        ic=153
!       cross section of pion- + n to rho- + n
        elseif(rlu1 > s2 .and. rlu1 <= s3)then
        ik1=-211
        ik2=2114
!       cross section of pion- + n to delta0 + pion-
        ic=154
        else
        ik1=311
        ik2=3112
        ic=14
!       cross section of pion- + n to k0 + sigma-
        endif
        lc(icp,3)=ik1
        lc(icp,4)=ik2
        lc(icp,5)=ic
        tw(icp)=sigm14/cspin/10.
        ioo=1


13      return
        end



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine pin2(ss,icp,ioo)
!       a part of 'prod' to deal with pion0 + n to ...
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        INTEGER PYCOMP
        PARAMETER(NSIZE=750000)
!       calculate particle production weight and fill up lc(i,3-5),tw(i).
!       tw : the ratio of cross section of (special inela.)/tot
        COMMON/PYDAT2/KCHG(500,4),PMAS(500,4),PARF(2000),VCKM(4,4)
        COMMON/PYCIDAT2/KFMAXT,NONT2,PARAM(20),WEIGH(600)
        SAVE /PYCIDAT2/
        common/sa10_h/csnn,cspin,cskn,cspipi,cspsn,cspsm,rcsit,ifram, &
         iabsb,iabsm,i_sigma_AQM,csspn,csspm
        common/count_h/isinel(600)
        common/iloval/ilo1,ilo2,ilo3,ilo4,ilo5,ilo6,ilo7,ilo8,ilo9
        common/ctllist_t/ lc(nsize,5), tc(nsize), tw(nsize)


        if(isinel(172) == 0)then
        sigma1=0.
        goto 216
        endif
        ilo1=1
        the=pmas(pycomp(-213),1)+pmas(pycomp(2212),1)
        sigma1=srho(ss,ilo1,0,the)*WEIGH(72)
!       cross section of pion0 + n to rho- + p
216     if(isinel(168) == 0)then
        sigma2=0.
        goto 217
        endif
!       the--threshold energy of a reaction
        the=pmas(pycomp(-211),1)+pmas(pycomp(2214),1)
        sigma2=sdelta(ss,ilo2,0,the)*WEIGH(59)
!       cross section of pion0 + n to delta+ + pi-
217     if(isinel(171) == 0)then
        sigma3=0.
        goto 218
        endif
        the=pmas(pycomp(113),1)+pmas(pycomp(2112),1)
        sigma3=srho(ss,ilo3,0,the)*WEIGH(70)
!       cross section of pion0 + n to rho0 + n
218     if(isinel(169) == 0)then
        sigma4=0.
        goto 807
        endif
        the=pmas(pycomp(211),1)+pmas(pycomp(1114),1)
        sigma4=sdelta(ss,ilo4,0,the)*WEIGH(60)
!       cross section of pion0 + n to delta- + pion+
807     if(isinel(19) == 0)then
        sigma5=0.
        goto 811
        endif
        ilo5=1
        the=pmas(pycomp(311),1)+pmas(pycomp(3122),1)

        sigma5=s1724(ss,ilo5,0,the)*WEIGH(10)

!       cross section of pion0 + n to k0 + lambda is assumed
!       to be cross section of pion0 + p to k+ + lambda
811     if(isinel(18) == 0)then
        sigma6=0.
        goto 816
        endif
        ilo6=1
        the=pmas(pycomp(321),1)+pmas(pycomp(3112),1)

        sigma6=s2325(ss,ilo6,0,the)*WEIGH(22)

!       cross section of pion0 + n to k+ + sigma-
816     if(isinel(20) == 0)then
        sigma7=0.
        goto 817
        endif
        ilo7=1
        the=pmas(pycomp(311),1)+pmas(pycomp(3212),1)

        sigma7=s1724(ss,ilo7,0,the)*WEIGH(23)

!       cross section of pion0 + n to k0 + sigma0

817     if(isinel(170) == 0)then
        sigma8=0.
        goto 827
        endif
        the=pmas(pycomp(111),1)+pmas(pycomp(2114),1)
        sigma8=sdelta(ss,ilo8,0,the)*WEIGH(60)
!       cross section of pion0 + n to delta0 + pion0
827     if(ilo1 == 0 .and. ilo2 == 0 .and. ilo3 == 0.and. &
         ilo4 == 0.and. ilo5 == 0 &
         .and. ilo6 == 0 .and. ilo7 == 0 .and. ilo8 == 0) goto 13
        if(sigma1 < 1.e-6 .and. sigma2 < 1.e-6 .and. sigma3 < 1.e-6 &
                .and. sigma4 < 1.e-6.and. sigma5 < 1.e-6 &
                .and. sigma6 < 1.e-6.and. sigma7 < 1.e-6 &
                .and. sigma8 < 1.e-6) goto 13
        sigma12=sigma1+sigma2
        sigma13=sigma12+sigma3
        sigma14=sigma13+sigma4
        sigma15=sigma14+sigma5
        sigma16=sigma15+sigma6
        sigma17=sigma16+sigma7
        sigma18=sigma17+sigma8
        s1=sigma1/sigma18
        s2=sigma12/sigma18
        s3=sigma13/sigma18
        s4=sigma14/sigma18
        s5=sigma15/sigma18
        s6=sigma16/sigma18
        s7=sigma17/sigma18
        rlus=pyr(1)
        if(rlus <= s1)then
        ik1=-213
        ik2=2212
!       -213 is the flavor code of rho-
        ic=172
!       pion0 + n to p + rho-
        elseif(rlus > s1 .and. rlus <= s2)then
        ik1=-211
        ik2=2214
!       2214 is the flavor code of delta+
        ic=168
!       pion0 + n to delta+ + pi-
        elseif(rlus > s2 .and. rlus <= s3)then
        ik1=113
        ik2=2112
!       113 is the Kf code of rho
        ic=171
!       pion0 + n to rho  +  n
        elseif(rlus > s3 .and. rlus <= s4)then
        ik1=211
        ik2=1114
!       1114 is the Kf code of delta-
        ic=169
        elseif(rlus > s4 .and. rlus <= s5)then
        ik1=311
        ik2=3122
!       3122 is the Kf code of lambda
!       pion0 + n to k0 + lambda
        ic=19
        elseif(rlus > s5 .and. rlus <= s6)then
        ik1=321
        ik2=3112
!       3112 is the Kf code of sigma-
        ic=18
        elseif(rlus > s6 .and. rlus <= s7)then
!       pion0 + n to k+ + sigma-
        ik1=311
        ik2=3212
!       3212 is the Kf code of sigma0
        ic=20
!       pion0 + n to k0 + sigma0
        else
        ik1=111
        ik2=2114
!       2114 is the Kf code of delta0
        ic=170
!       pion0 + n to pi0 + delta0
        endif
        lc(icp,3)=ik1
        lc(icp,4)=ik2
        lc(icp,5)=ic
        tw(icp)=sigma18/cspin/10.
        ioo=1


13      return
        end



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine pin1(ss,icp,ioo)
!       a part of 'prod' to deal with pion+ + n to ...
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        INTEGER PYCOMP
        PARAMETER(NSIZE=750000)
!       calculate particle production weight and fill up lc(i,3-5),tw(i).
!       tw : the ratio of cross section of (special inela.)/tot
        COMMON/PYDAT2/KCHG(500,4),PMAS(500,4),PARF(2000),VCKM(4,4)
        COMMON/PYCIDAT2/KFMAXT,NONT2,PARAM(20),WEIGH(600)
        SAVE /PYCIDAT2/
        common/sa10_h/csnn,cspin,cskn,cspipi,cspsn,cspsm,rcsit,ifram, &
         iabsb,iabsm,i_sigma_AQM,csspn,csspm
        common/count_h/isinel(600)
        common/iloval/ilo1,ilo2,ilo3,ilo4,ilo5,ilo6,ilo7,ilo8,ilo9
        common/ctllist_t/ lc(nsize,5), tc(nsize), tw(nsize)


        if(isinel(8) == 0)then
        sigma1=0.
        goto 306
        endif
        the=pmas(pycomp(321),1)+pmas(pycomp(3212),1)
        sigma1=s1724(ss,ilo1,0,the)*WEIGH(18)
!       cross section of pion+ + n to k+ + sigma0
306     if(isinel(9) == 0)then
        sigma2=0.
        goto 207
        endif
        the=pmas(pycomp(321),1)+pmas(pycomp(3122),1)
        sigma2=s1727(ss,ilo2,0,the)*WEIGH(8)
!       cross section of pion+ + n to k+ + lambda
207     if(isinel(158) == 0)then
        sigma3=0.
        goto 208
        endif
!       the: threshold energy of a reaction
        the=pmas(pycomp(-211),1)+pmas(pycomp(2224),1)
        sigma3=sdelta(ss,ilo3,0,the)*WEIGH(55)
!       cross section of pion+ + n to delta++ + pi-
208     if(isinel(161) == 0)then
        sigma4=0.
        goto 209
        endif
        the=pmas(pycomp(113),1)+pmas(pycomp(2212),1)
        sigma4=srho(ss,ilo4,0,the)*WEIGH(66)
!       cross section of pion+ + n to rho0 + p
209     if(isinel(162) == 0)then
        sigma5=0.
        goto 210
        endif
        the=pmas(pycomp(213),1)+pmas(pycomp(2112),1)
        sigma5=srho(ss,ilo5,0,the)*WEIGH(67)
!       cross section of pion+ + n to rho+ + n
210     if(isinel(159) == 0)then
        sigma6=0.
        goto 806
        endif
        the=pmas(pycomp(211),1)+pmas(pycomp(2114),1)
        sigma6=sdelta(ss,ilo6,0,the)*WEIGH(56)
!       cross section of pion+ + n to delta0 + pion+
806     if(isinel(10) == 0)then
        sigma7=0.
        goto 814
        endif
        ilo7=1
        the=pmas(pycomp(311),1)+pmas(pycomp(3222),1)

        sigma7=s1724(ss,ilo7,0,the)*WEIGH(19)

!       cross section of pion+ + n to k0 + sigma+

814     if(isinel(160) == 0)then
        sigma8=0.
        goto 818
        endif
        ilo8=1
        the=pmas(pycomp(111),1)+pmas(pycomp(2214),1)

        sigma8=sdelta(ss,ilo8,0,the)*WEIGH(19)

!       cross section of pion+ + n to pi0 + delta+

818     if(ilo1 == 0 .and. ilo2 == 0 .and. ilo3 == 0 &
         .and.ilo4 == 0 .and. ilo5 == 0.and.ilo6 == 0 &
         .and.ilo7 == 0.and.ilo8 == 0) goto 13
        if(sigma1 < 1.e-6 .and. sigma2 < 1.e-6 .and. sigma3 < 1.e-6 &
               .and.sigma4 < 1.e-6 .and. sigma5 < 1.e-6.and.  &
                sigma6 < 1.e-6 &
               .and.sigma7 < 1.e-6.and.sigma8 < 1.e-6) goto 13
        sigma12=sigma1+sigma2
        sigma13=sigma12+sigma3
        sigma14=sigma13+sigma4
        sigma15=sigma14+sigma5
        sigma16=sigma15+sigma6
        sigma17=sigma16+sigma7
        sigma18=sigma17+sigma8
        s1=sigma1/sigma18
        s2=sigma12/sigma18
        s3=sigma13/sigma18
        s4=sigma14/sigma18
        s5=sigma15/sigma18
        s6=sigma16/sigma18
        s7=sigma17/sigma18
        rlus=pyr(1)
        if(rlus <= s1)then
        ik1=321
        ik2=3212
!       3212 is the flavor code of sigma0
        ic=8
!       pion+ + n to k+ + sigma0
        elseif(rlus > s1 .and. rlus <= s2)then
        ik1=321
        ik2=3122
!       3122 is the flavor code of lambda
        ic=9
!       pion+ + n to k+ + lambda
        elseif(rlus > s2 .and. rlus <= s3)then
        ik1=-211
        ik2=2224
!       2224 is the Kf code of delta++
        ic=158
!       pion+ + n to pi- + delta++
        elseif(rlus > s3 .and. rlus <= s4)then
        ik1=113
        ik2=2212
!       113 is the Kf code of rho0
        ic=161
!       pion+ + n to rho0 + p
        elseif(rlus > s4 .and. rlus <= s5)then
        ik1=213
        ik2=2112
!       213 is the Kf code of rho+
        ic=162
!       pion+ + n to rho+ + n
        elseif(rlus > s5 .and. rlus <= s6)then
        ik1=211
        ik2=2114
!       2114 is the Kf code of delta0
        ic=159
        elseif(rlus > s6 .and. rlus <= s7)then
        ik1=311
        ik2=3222
!       3222 is the Kf code of sigma+
        ic=10
!       pion+ + n to k0 + sigma+
        else
        ik1=111
        ik2=2214
!       2214 is the Kf code of delta+
        ic=160
!       pion+ + n to pi0 + delta+
        endif
        lc(icp,3)=ik1
        lc(icp,4)=ik2
        lc(icp,5)=ic
        tw(icp)=sigma18/cspin/10.
        ioo=1


13      return
        end



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine pip2(ss,icp,ioo)
!       a part of 'prod' to deal with pion+ + p to ...
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        INTEGER PYCOMP
        PARAMETER(NSIZE=750000)
!       calculate particle production weight and fill up lc(i,3-5),tw(i).
!       tw : the ratio of cross section of (special inela.)/tot
        COMMON/PYDAT2/KCHG(500,4),PMAS(500,4),PARF(2000),VCKM(4,4)
        COMMON/PYCIDAT2/KFMAXT,NONT2,PARAM(20),WEIGH(600)
        SAVE /PYCIDAT2/
        common/sa10_h/csnn,cspin,cskn,cspipi,cspsn,cspsm,rcsit,ifram, &
         iabsb,iabsm,i_sigma_AQM,csspn,csspm
        common/count_h/isinel(600)
        common/iloval/ilo1,ilo2,ilo3,ilo4,ilo5,ilo6,ilo7,ilo8,ilo9
        common/ctllist_t/ lc(nsize,5), tc(nsize), tw(nsize)


        if(isinel(7) == 0)then
        sigma1=0.
        goto 212
        endif
        the=pmas(pycomp(321),1)+pmas(pycomp(3222),1)
        sigma1=s1713(ss,ilo1,0,the)*WEIGH(17)
!       cross section of pion+ + p to k+ + sigma+
212     if(isinel(155) == 0)then
        sigma2=0.
        goto 213
        endif
!       the: threshold energy of a reaction
        the=pmas(pycomp(111),1)+pmas(pycomp(2224),1)
        sigma2=sdelta(ss,ilo2,0,the)*WEIGH(51)
!       cross section of pion+ + p to delta++ + pi0
213     if(isinel(157) == 0)then
        sigma3=0.
        goto 214
        endif
        the=pmas(pycomp(213),1)+pmas(pycomp(2212),1)
        sigma3=srho(ss,ilo3,0,the)*WEIGH(68)
!       cross section of pion+ + p to rho+ + p
214     if(isinel(156) == 0)then
        sigma4=0.
        goto 804
        endif
        the=pmas(pycomp(211),1)+pmas(pycomp(2214),1)
        sigma4=sdelta(ss,ilo4,0,the)*WEIGH(52)
!       cross section of pion+ + p to delta+ + pion+
804     if(ilo1 == 0 .and. ilo2 == 0 .and. ilo3 == 0 &
         .and. ilo4 == 0) goto 13
        if(sigma1 < 1.e-6 .and. sigma2 < 1.e-6 .and. sigma3 < 1.e-6 &
         .and. sigma4 < 1.e-6 ) goto 13
        sigma12=sigma1+sigma2
        sigma13=sigma12+sigma3
        sigma14=sigma13+sigma4
        s1=sigma1/sigma14
        s2=sigma12/sigma14
        s3=sigma13/sigma14
        rlus=pyr(1)
        if(rlus <= s1)then
        ik1=321
        ik2=3222
!       3222 is the flavor code of sigma+
        ic=7
!       pion+ + p to k+ + sigma+
        elseif(rlus > s1 .and. rlus <= s2)then
        ik1=111
        ik2=2224
!       2224 is the flavor code of delta++
        ic=155
!       pion+ + p to delta++  +  pi0
        elseif(rlus > s2 .and. rlus <= s3)then
        ik1=213
        ik2=2212
!       213 is the Kf code of rho+
        ic=157
!       pion+ + p to rho+ + p
        else
        ik1=211
        ik2=2214
!       2214 is the Kf code of delta+
        ic=156
!       pion+ + p to delta+ + pion+
        endif
        lc(icp,3)=ik1
        lc(icp,4)=ik2
        lc(icp,5)=ic
        tw(icp)=sigma14/cspin/10.
        ioo=1


13      return
        end



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        function sdelta(ss,ilo,i,the)
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        common/sa13/kjp20,non13,vjp20,vjp21,vjp22,vjp23
!       pion- + p to pi0 + delta0
        ilo=0
        si=0.
        ii=i
        if(ii == 1) goto 30
        if(ss <= the) goto 100
30      ilo=1
        if(ss >= 1.6941) goto 10
        IF(kjp20 == 0)THEN
        si=-61.127+42.9365*ss
        ELSE
        si=vjp21
        ENDIF
        goto 100
10      IF(kjp20 == 0)THEN
        sit=-0.0186959*ss**3+0.310359*ss**2-0.755106*ss+0.565481
        si=1.0/sit
        ELSE
        si=vjp21
        ENDIF
100     sdelta=si
        return
        end



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        function srho(ss,ilo,i,the)
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        common/sa13/kjp20,non13,vjp20,vjp21,vjp22,vjp23
!       pion- + p to pi0 + delta0
        ilo=0
        si=0.
        ii=i
        if(ii == 1) goto 30
        if(ss <= the) goto 100
30      ilo=1
        if(ss >= 1.8837) goto 10
        IF(kjp20 == 0)THEN
        si=-23.3607+13.9936*ss
        ELSE
        si=vjp22
        ENDIF
        goto 100
10      IF(kjp20 == 0)THEN
        sit=0.331583*ss**3-1.86123*ss**2+3.81364*ss-2.50068
        si=1.0/sit
        ELSE
        si=vjp22
        ENDIF
100     srho=si
        return
        end



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        function snn(ss,ilo,i,the)
!       n + n to n + delta
!       parameterized x-section is not correct for n-n
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        common/sa13/kjp20,non13,vjp20,vjp21,vjp22,vjp23
        ilo=0
        si=0.
        ii=i
        if(ii == 1) goto 30
        if(ss <= the) goto 100
30      ilo=1
        IF(kjp20 == 0)THEN
        si=20*(ss-2.015)**2/(0.015+(ss-2.015)**2)
        ELSE
        si=vjp23
        ENDIF
100     snn=si
        return
        end



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine pip3(ss,icp,ioo)
!       a part of 'prod' to deal with pion0 + p to ...
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        INTEGER PYCOMP
        PARAMETER(NSIZE=750000)
!       calculate particle production weight and fill up lc(i,3-5),tw(i).
!       tw : the ratio of cross section of (special inela.)/tot
        COMMON/PYDAT2/KCHG(500,4),PMAS(500,4),PARF(2000),VCKM(4,4)
        COMMON/PYCIDAT2/KFMAXT,NONT2,PARAM(20),WEIGH(600)
        SAVE /PYCIDAT2/
        common/sa10_h/csnn,cspin,cskn,cspipi,cspsn,cspsm,rcsit,ifram, &
         iabsb,iabsm,i_sigma_AQM,csspn,csspm
        common/count_h/isinel(600)
        common/iloval/ilo1,ilo2,ilo3,ilo4,ilo5,ilo6,ilo7,ilo8,ilo9
        common/ctllist_t/ lc(nsize,5), tc(nsize), tw(nsize)


        if(isinel(15) == 0)then
        sigma1=0.
        goto 308
        endif
        the=pmas(pycomp(321),1)+pmas(pycomp(3212),1)
        sigma1=s2314(ss,ilo1,0,the)*WEIGH(20)
!       cross section of pion0 + p to k+ + sigma0
308     if(isinel(16) == 0)then
        sigma2=0.
        goto 239
        endif
        the=pmas(pycomp(321),1)+pmas(pycomp(3122),1)
        sigma2=s2317(ss,ilo2,0,the)*WEIGH(9)
!       cross section of pion0 + p to k+ + lambda
239     if(isinel(163) == 0)then
        sigma3=0.
        goto 220
        endif
!       the: threshold energy of a reaction
        the=pmas(pycomp(211),1)+pmas(pycomp(2114),1)
        sigma3=sdelta(ss,ilo3,0,the)*WEIGH(53)
!       cross section of pion0 + p to delta0 + pi+
220     if(isinel(165) == 0)then
        sigma4=0.
        goto 240
        endif
        the=pmas(pycomp(213),1)+pmas(pycomp(2112),1)
        sigma4=srho(ss,ilo4,0,the)*WEIGH(69)
!       cross section of pion0 + p to rho+ + n
240     if(isinel(164) == 0)then
        sigma5=0.
        goto 808
        endif
        the=pmas(pycomp(-211),1)+pmas(pycomp(2224),1)
        sigma5=sdelta(ss,ilo5,0,the)*WEIGH(54)
!       cross section of pion0 + p to delta++ + pion-
808     if(isinel(17) == 0)then
        sigma6=0.
        goto 815
        endif
        ilo6=1
        the=pmas(pycomp(311),1)+pmas(pycomp(3222),1)

        sigma6=s1724(ss,ilo6,0,the)*WEIGH(21)

!       cross section of pion0 + P to k0 + sigma+

815     if(isinel(166) == 0)then
        sigma7=0.
        goto 835
        endif
        ilo7=1
        the=pmas(pycomp(113),1)+pmas(pycomp(2212),1)
        sigma7=srho(ss,ilo7,0,the)*WEIGH(71)
!       cross section of pion0 + P to rho0 + p
835     if(isinel(167) == 0)then
        sigma8=0.
        goto 845
        endif
        ilo8=1
        the=pmas(pycomp(111),1)+pmas(pycomp(2214),1)
        sigma8=sdelta(ss,ilo8,0,the)*WEIGH(54)
!       cross section of pion0 + P to pi0 + delta+
845     if(ilo1 == 0 .and. ilo2 == 0 .and. ilo3 == 0 &
         .and.ilo4 == 0.and.ilo5 == 0.and.ilo6 == 0 &
         .and.ilo7 == 0.and.ilo8 == 0) goto 13
        if(sigma1 < 1.e-6 .and. sigma2 < 1.e-6 .and. sigma3 < 1.e-6 &
         .and.sigma4 < 1.e-6.and.sigma5 < 1.e-6.and.sigma6 < 1.e-6 &
         .and.sigma7 < 1.e-6.and.sigma8 < 1.e-6 ) goto 13
        sigma12=sigma1+sigma2
        sigma13=sigma12+sigma3
        sigma14=sigma13+sigma4
        sigma15=sigma14+sigma5
        sigma16=sigma15+sigma6
        sigma17=sigma16+sigma7
        sigma18=sigma17+sigma8
        s1=sigma1/sigma18
        s2=sigma12/sigma18
        s3=sigma13/sigma18
        s4=sigma14/sigma18
        s5=sigma15/sigma18
        s6=sigma16/sigma18
        s7=sigma17/sigma18
        rlus=pyr(1)
        if(rlus <= s1)then
        ik1=321
        ik2=3212
        ic=15
!       3212 is the flavor code of sigma0
!       pion0 + p to k+ + sigma0
        elseif(rlus > s1 .and. rlus <= s2)then
        ik1=321
        ik2=3122
        ic=16
!       3122 is the flavor code of lambda
!       pion0 + p to k+ + lambda
        elseif(rlus > s2 .and. rlus <= s3)then
        ik1=211
        ik2=2114
!       2214 is the Kf code of delta+
        ic=163
!       pion0 + p to pi+ + delta0
        elseif(rlus > s3 .and. rlus <= s4)then
        ik1=213
        ik2=2112
!       213 is the Kf code of rho+
        ic=165
!       pion0 + p to rho+ + n
        elseif(rlus > s4 .and. rlus <= s5)then
        ik1=-211
        ik2=2224
!       2224 is the Kf code of delta++
        ic=164
        elseif(rlus > s5 .and. rlus <= s6)then
        ik1=311
        ik2=3222
!       3222 is the Kf code of sigma+
!       pion0 + p to k0 + sigma+
        ic=17
        elseif(rlus > s6 .and. rlus <= s7)then
        ik1=113
        ik2=2212
!       113 is the Kf code of rho0
        ic=166
!       pion0 + p to rho0 + p
        else
        ik1=111
        ik2=2214
!       2214 is the Kf code of delta+
!       pion0 + p to pi0 + delta+
        ic=167
        endif
        lc(icp,3)=ik1
        lc(icp,4)=ik2
        lc(icp,5)=ic
        tw(icp)=sigma18/cspin/10.
        ioo=1


13      return
        end



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine jpsin(ss,icp,ioo,ii)
!       a part of 'prod' to deal with J/psi (psi') + n to ...
!       ii = 1 for J/psi, ii = 2 for psi'
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        INTEGER PYCOMP
        PARAMETER(NSIZE=750000)
        COMMON/PYDAT2/KCHG(500,4),PMAS(500,4),PARF(2000),VCKM(4,4)
        COMMON/PYCIDAT2/KFMAXT,NONT2,PARAM(20),WEIGH(600)
        SAVE /PYCIDAT2/
        common/sa10_h/csnn,cspin,cskn,cspipi,cspsn,cspsm,rcsit,ifram, &
         iabsb,iabsm,i_sigma_AQM,csspn,csspm
        common/count_h/isinel(600)
        common/iloval/ilo1,ilo2,ilo3,ilo4,ilo5,ilo6,ilo7,ilo8,ilo9
        common/ctllist_t/ lc(nsize,5), tc(nsize), tw(nsize)


        if(isinel(179) == 0)then
        sigma1=0.
        goto 212
        endif
        the=pmas(pycomp(4122),1)+pmas(pycomp(411),1)
!       the--threshold energy of a reaction (squart root of s)
        if(ss <= the)then
        sigma1=0.
        goto 212
        endif
        sigma1=cspsn
!       cross section of J/psi + n (fm^2)
        if(ii == 2)sigma1=csspn
!       cross section of psi' + n (fm^2)
212     if(isinel(180) == 0)then
        sigma2=0.
        goto 213
        endif
        the=pmas(pycomp(4212),1)+pmas(pycomp(411),1)
        if(ss <= the)then
        sigma2=0.
        goto 213
        endif
        sigma2=cspsn
        if(ii == 2)sigma2=csspn
213     if(isinel(181) == 0)then
        sigma3=0.
        goto 214
        endif
        the=pmas(pycomp(4112),1)+pmas(pycomp(421),1)
        if(ss <= the)then
        sigma3=0.
        goto 214
        endif
        sigma3=cspsn
        if(ii == 2)sigma3=csspn
214     if(sigma1 < 1.e-16 .and. sigma2 < 1.e-16 .and. &
          sigma3 < 1.e-16) goto 13
        sigma12=sigma1+sigma2
        sigma13=sigma12+sigma3
        s1=sigma1/sigma13
        s2=sigma12/sigma13
        rlus=pyr(1)
        if(rlus <= s1)then
        ik1=4122
        ik2=-411
        ic=179
!       J/psi (psi') + n to lamdac + Dba
        elseif(rlus > s1 .and. rlus <= s2)then
        ik1=4212
        ik2=-411
        ic=180
!       J/psi (psi') + n to sigmac+  + Dba
        else
        ik1=4112
        ik2=-421
        ic=181
!       J/psi (psi') + n to sigmac0 + D0ba
        endif
        lc(icp,3)=ik1
        lc(icp,4)=ik2
        lc(icp,5)=ic
        tw(icp)=sigma13*rcsit/cspsn
        if(ii == 2)tw(icp)=sigma13*rcsit/csspn
        ioo=1


13      return
        end



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine jpsip(ss,icp,ioo,ii)
!       a part of 'prod' to deal with J/psi (psi') + p to ...
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        INTEGER PYCOMP
        PARAMETER(NSIZE=750000)
        COMMON/PYDAT2/KCHG(500,4),PMAS(500,4),PARF(2000),VCKM(4,4)
        COMMON/PYCIDAT2/KFMAXT,NONT2,PARAM(20),WEIGH(600)
        SAVE /PYCIDAT2/
        common/sa10_h/csnn,cspin,cskn,cspipi,cspsn,cspsm,rcsit,ifram, &
         iabsb,iabsm,i_sigma_AQM,csspn,csspm
        common/count_h/isinel(600)
        common/iloval/ilo1,ilo2,ilo3,ilo4,ilo5,ilo6,ilo7,ilo8,ilo9
        common/ctllist_t/ lc(nsize,5), tc(nsize), tw(nsize)


        if(isinel(182) == 0)then
        sigma1=0.
        goto 212
        endif
        the=pmas(pycomp(4122),1)+pmas(pycomp(421),1)
!       the: threshold energy of a reaction (squart root of s)
        if(ss <= the)then
        sigma1=0.
        goto 212
        endif
        sigma1=cspsn
        if(ii == 2)sigma1=csspn
212     if(isinel(183) == 0)then
        sigma2=0.
        goto 213
        endif
        the=pmas(pycomp(4212),1)+pmas(pycomp(421),1)
        if(ss <= the)then
        sigma2=0.
        goto 213
        endif
        sigma2=cspsn
        if(ii == 2)sigma2=csspn
213     if(isinel(184) == 0)then
        sigma3=0.
        goto 214
        endif
        the=pmas(pycomp(4222),1)+pmas(pycomp(411),1)
        if(ss <= the)then
        sigma3=0.
        goto 214
        endif
        sigma3=cspsn
        if(ii == 2)sigma3=csspn
214     if(sigma1 < 1.e-16 .and. sigma2 < 1.e-16 .and. &
          sigma3 < 1.e-16) goto 13
        sigma12=sigma1+sigma2
        sigma13=sigma12+sigma3
        s1=sigma1/sigma13
        s2=sigma12/sigma13
        rlus=pyr(1)
        if(rlus <= s1)then
        ik1=4122
        ik2=-421
        ic=182
!       J/psi (psi') + p to lamdac + D0ba
        elseif(rlus > s1 .and. rlus <= s2)then
        ik1=4212
        ik2=-421
        ic=183
!       J/psi (psi') + p to sigmac+  + D0ba
        else
        ik1=4222
        ik2=-411
        ic=184
!       J/psi (psi') + p to sigmac++ + Dba
        endif
        lc(icp,3)=ik1
        lc(icp,4)=ik2
        lc(icp,5)=ic
        tw(icp)=sigma13*rcsit/cspsn
        if(ii == 2)tw(icp)=sigma13*rcsit/csspn
        ioo=1


13      return
        end



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine jpsip1(ss,icp,ioo,ii)
!       a part of 'prod' to deal with J/psi (psi') + pion+ to ...
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        INTEGER PYCOMP
        PARAMETER(NSIZE=750000)
        COMMON/PYDAT2/KCHG(500,4),PMAS(500,4),PARF(2000),VCKM(4,4)
        COMMON/PYCIDAT2/KFMAXT,NONT2,PARAM(20),WEIGH(600)
        SAVE /PYCIDAT2/
        common/sa10_h/csnn,cspin,cskn,cspipi,cspsn,cspsm,rcsit,ifram, &
         iabsb,iabsm,i_sigma_AQM,csspn,csspm
        common/count_h/isinel(600)
        common/iloval/ilo1,ilo2,ilo3,ilo4,ilo5,ilo6,ilo7,ilo8,ilo9
        common/ctllist_t/ lc(nsize,5), tc(nsize), tw(nsize)


        if(isinel(185) == 0) goto 13
        the=pmas(pycomp(411),1)+pmas(pycomp(423),1)
        if(ss <= the) goto 13
!98/03/23   ww=cspsm
!       J/psi + pion+ to D + D*0ba
        lc(icp,3)=411
        lc(icp,4)=-423
        lc(icp,5)=185
        tw(icp)=rcsit
        ioo=1
        jj=ii


13      return
        end



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine jpsip0(ss,icp,ioo,ii)
!       a part of 'prod' to deal with J/psi (psi') + pion0 to ...
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        INTEGER PYCOMP
        PARAMETER(NSIZE=750000)
        COMMON/PYDAT2/KCHG(500,4),PMAS(500,4),PARF(2000),VCKM(4,4)
        COMMON/PYCIDAT2/KFMAXT,NONT2,PARAM(20),WEIGH(600)
        SAVE /PYCIDAT2/
        common/sa10_h/csnn,cspin,cskn,cspipi,cspsn,cspsm,rcsit,ifram, &
         iabsb,iabsm,i_sigma_AQM,csspn,csspm
        common/count_h/isinel(600)
        common/iloval/ilo1,ilo2,ilo3,ilo4,ilo5,ilo6,ilo7,ilo8,ilo9
        common/ctllist_t/ lc(nsize,5), tc(nsize), tw(nsize)


        if(isinel(186) == 0)then
        sigma1=0.
        goto 212
        endif
        the=pmas(pycomp(421),1)+pmas(pycomp(423),1)
!       the: threshold energy of a reaction (squart root of s)
        if(ss <= the)then
        sigma1=0.
        goto 212
        endif
        sigma1=cspsm
!       cross section of J/psi + m (fm^2)
        if(ii == 2)sigma1=csspm
!       cross section of psi' + m (fm^2)
212     if(isinel(187) == 0)then
        sigma2=0.
        goto 213
        endif
        the=pmas(pycomp(411),1)+pmas(pycomp(413),1)
        if(ss <= the)then
        sigma2=0.
        goto 213
        endif
        sigma2=cspsm
        if(ii == 2)sigma2=csspm
213     if(sigma1 < 1.e-16 .and. sigma2 < 1.e-16) goto 13
        sigma12=sigma1+sigma2
        s1=sigma1/sigma12
        rlus=pyr(1)
        if(rlus <= s1)then
        ik1=421
        ik2=-423
        ic=186
!       J/psi (psi') + pion0 to D0 + D*0ba
        else
        ik1=411
        ik2=-413
        ic=187
!       J/psi (psi') + pion0 to D + D*ba
        endif
        lc(icp,3)=ik1
        lc(icp,4)=ik2
        lc(icp,5)=ic
        tw(icp)=sigma12*rcsit/cspsm
        if(ii == 2)tw(icp)=sigma12*rcsit/csspm
        ioo=1


13      return
        end



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine jpsip2(ss,icp,ioo,ii)
!       a part of 'prod' to deal with J/psi (psi') + pion- to ...
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        INTEGER PYCOMP
        PARAMETER(NSIZE=750000)
        COMMON/PYDAT2/KCHG(500,4),PMAS(500,4),PARF(2000),VCKM(4,4)
        COMMON/PYCIDAT2/KFMAXT,NONT2,PARAM(20),WEIGH(600)
        SAVE /PYCIDAT2/
        common/sa10_h/csnn,cspin,cskn,cspipi,cspsn,cspsm,rcsit,ifram, &
         iabsb,iabsm,i_sigma_AQM,csspn,csspm
        common/count_h/isinel(600)
        common/iloval/ilo1,ilo2,ilo3,ilo4,ilo5,ilo6,ilo7,ilo8,ilo9
        common/ctllist_t/ lc(nsize,5), tc(nsize), tw(nsize)


        if(isinel(188) == 0) goto 13
        the=pmas(pycomp(421),1)+pmas(pycomp(413),1)
        if(ss <= the) goto 13
!98/03/23   ww=cspsm
!       J/psi (psi') + pion- to D0 + D*ba
        lc(icp,3)=421
        lc(icp,4)=-413
        lc(icp,5)=188
        tw(icp)=rcsit
        ioo=1
        jj=ii


13      return
        end



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine jpsir1(ss,icp,ioo,ii)
!       a part of 'prod' to deal with J/psi (psi') + rho+ to ...
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        INTEGER PYCOMP
        PARAMETER(NSIZE=750000)
        COMMON/PYDAT2/KCHG(500,4),PMAS(500,4),PARF(2000),VCKM(4,4)
        COMMON/PYCIDAT2/KFMAXT,NONT2,PARAM(20),WEIGH(600)
        SAVE /PYCIDAT2/
        common/sa10_h/csnn,cspin,cskn,cspipi,cspsn,cspsm,rcsit,ifram, &
         iabsb,iabsm,i_sigma_AQM,csspn,csspm
        common/count_h/isinel(600)
        common/iloval/ilo1,ilo2,ilo3,ilo4,ilo5,ilo6,ilo7,ilo8,ilo9
        common/ctllist_t/ lc(nsize,5), tc(nsize), tw(nsize)


        if(isinel(189) == 0) goto 13
        the=pmas(pycomp(411),1)+pmas(pycomp(421),1)
        if(ss <= the) goto 13
!98/03/23   ww=cspsm
!       J/psi (psi,) + rho+ to D + D0ba
        lc(icp,3)=411
        lc(icp,4)=-421
        lc(icp,5)=189
        tw(icp)=rcsit
        ioo=1
        jj=ii


13      return
        end



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine jpsir0(ss,icp,ioo,ii)
!       a part of 'prod' to deal with J/psi (psi') + rho0 to ...
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        INTEGER PYCOMP
        PARAMETER(NSIZE=750000)
        COMMON/PYDAT2/KCHG(500,4),PMAS(500,4),PARF(2000),VCKM(4,4)
        COMMON/PYCIDAT2/KFMAXT,NONT2,PARAM(20),WEIGH(600)
        SAVE /PYCIDAT2/
        common/sa10_h/csnn,cspin,cskn,cspipi,cspsn,cspsm,rcsit,ifram, &
         iabsb,iabsm,i_sigma_AQM,csspn,csspm
        common/count_h/isinel(600)
        common/iloval/ilo1,ilo2,ilo3,ilo4,ilo5,ilo6,ilo7,ilo8,ilo9
        common/ctllist_t/ lc(nsize,5), tc(nsize), tw(nsize)


        if(isinel(190) == 0)then
        sigma1=0.
        goto 212
        endif
        the=pmas(pycomp(421),1)*2.
!       the--threshold energy of a reaction (squart root of s)
        if(ss <= the)then
        sigma1=0.
        goto 212
        endif
        sigma1=cspsm
        if(ii == 2)sigma1=csspm
212     if(isinel(191) == 0)then
        sigma2=0.
        goto 213
        endif
        the=pmas(pycomp(411),1)*2.
        if(ss <= the)then
        sigma2=0.
        goto 213
        endif
        sigma2=cspsm
        if(ii == 2)sigma2=csspm
213     if(sigma1 < 1.e-16 .and. sigma2 < 1.e-16) goto 13
        sigma12=sigma1+sigma2
        s1=sigma1/sigma12
        rlus=pyr(1)
        if(rlus <= s1)then
        ik1=421
        ik2=-421
        ic=190
!       J/psi (psi') + rho0 to D0 + D0ba
        else
        ik1=411
        ik2=-411
        ic=191
!       J/psi (psi') + rho0 to D + Dba
        endif
        lc(icp,3)=ik1
        lc(icp,4)=ik2
        lc(icp,5)=ic
        tw(icp)=2.*sigma12*rcsit/cspsm
        if(ii == 2)tw(icp)=2.*sigma12*rcsit/csspm
!       at the case of rho0, cross section enlarges a factor 2 to
!        consider the effect of omega
        ioo=1


13      return
        end



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine jpsir2(ss,icp,ioo,ii)
!       a part of 'prod' to deal with J/psi (psi') + rho- to ...
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        INTEGER PYCOMP
        PARAMETER(NSIZE=750000)
        COMMON/PYDAT2/KCHG(500,4),PMAS(500,4),PARF(2000),VCKM(4,4)
        COMMON/PYCIDAT2/KFMAXT,NONT2,PARAM(20),WEIGH(600)
        SAVE /PYCIDAT2/
        common/sa10_h/csnn,cspin,cskn,cspipi,cspsn,cspsm,rcsit,ifram, &
         iabsb,iabsm,i_sigma_AQM,csspn,csspm
        common/count_h/isinel(600)
        common/iloval/ilo1,ilo2,ilo3,ilo4,ilo5,ilo6,ilo7,ilo8,ilo9
        common/ctllist_t/ lc(nsize,5), tc(nsize), tw(nsize)


        if(isinel(192) == 0) goto 13
        the=pmas(pycomp(421),1)+pmas(pycomp(411),1)
        if(ss <= the) goto 13
!98/03/23   ww=cspsm
!       J/psi + rho- to D0 + Dba
        lc(icp,3)=421
        lc(icp,4)=-411
        lc(icp,5)=192
        tw(icp)=rcsit
        ioo=1
        jj=ii


13      return
        end



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine DMeson_coll( KF_A,KF_B,E_AB,icp,ioo )
!!      A part of "prod" to deal with D + pion/rho --> ...
!       Find production and fill up lc(i,3-5), tw(i).
!       D + pion/rho channels are defined as 385-458.
!       D, M: D meson and pion/rho Meson.
!       KF_A, KF_B: the KF code.
!       E_AB: the total energy of two hadrons.
!       icp: the line number of the colliding pair in collision (time) list.
!       lc: collision list.
!       tc: collision time list.
!       tw: the ratio of cross section of (special inela.)/tot.
!       ioo: flag of elastic (0) and inelastic (1) treatment.
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        INTEGER PYCHGE
        PARAMETER(NSIZE=750000)
        COMMON/PYCIDAT2/KFMAXT,NONT2,PARAM(20),WEIGH(600)
        common/sa10_h/csnn,cspin,cskn,cspipi,cspsn,cspsm,rcsit,ifram, &
         iabsb,iabsm,i_sigma_AQM,csspn,csspm
        common/count_h/isinel(600)
        common/ctllist_t/ lc(nsize,5), tc(nsize), tw(nsize)


!       Make order as D + pion/rho.
        i_order = 1
        KF_D = KF_A   ! KF of D meson.
        KF1_D = KF_D / 100   ! 280224 Lei The 100-th position of KF_D.
        KF_M = KF_B   ! KF of pion/rho meson.

!       Exchange orders if first one is not D meson.
        if( KF1_D /= 4 )then
            KF   = KF_D
            KF_D = KF_M
            KF_M = KF
            i_order = -1
        end if
!       Charges of two hadrons.
        iCharge_D = PYCHGE( KF_D )
        iCharge_M = PYCHGE( KF_M )

!       Single production channels: 16 channels, 385-400.
        if( (iCharge_D > 0 .AND. iCharge_M > 0) .OR.  &  ! ++
            (iCharge_D < 0 .AND. iCharge_M < 0) .OR.  &  ! --
            (iCharge_D == 0 .AND. iCharge_M > 0 .AND. &   ! 0+, KF-+
                  KF_D < 0 .AND.      KF_M > 0) .OR.  &  ! 0+, KF-+
            (iCharge_D == 0 .AND. iCharge_M < 0 .AND. &   ! 0-, KF+-
                  KF_D > 0 .AND.      KF_M < 0)       &  ! 0-, KF+-
          )then
            call Dmeson1(KF_D,KF_M,E_AB,icp,ioo,i_order)
!       Dual production channels: 64 channels, 401-464,
!        else if( (iCharge_D == 0 .AND. iCharge_M == 0) .OR.  &   ! 00
!                 (iCharge_D > 0  .AND. iCharge_M < 0)  .OR.  &   ! +-
!                 (iCharge_D < 0  .AND. iCharge_M > 0)  .OR.  &   ! -+
!                 (iCharge_D > 0  .AND. iCharge_M == 0) .OR.  &   ! +0
!                 (iCharge_D < 0  .AND. iCharge_M == 0) .OR.  &   ! -0
!                 (iCharge_D == 0 .AND. iCharge_M > 0   .AND. &   ! 0+, KF++
!                       KF_D > 0  .AND.      KF_M > 0)  .OR.  &   ! 0+, KF++
!                 (iCharge_D == 0 .AND. iCharge_M < 0   .AND. &   ! 0-, KF--
!                       KF_D < 0  .AND.      KF_M < 0)            ! 0-, KF--
!               )then
        else
            call Dmeson2(KF_D,KF_M,E_AB,icp,ioo,i_order)
        end if


        return
        end



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine Dmeson1(KF_D,KF_M,E_AB,icp,ioo,i_od)
!       A part of "prod" to deal with D + pion/rho  --> ...
!       Single production channel: 385-400, total 16 reaction channels.
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        PARAMETER(NSIZE=750000)
        common/sa10_h/csnn,cspin,cskn,cspipi,cspsn,cspsm,rcsit,ifram, &
         iabsb,iabsm,i_sigma_AQM,csspn,csspm
        common/count_h/isinel(600)
        common/ctllist_t/ lc(nsize,5), tc(nsize), tw(nsize)


!       385. D+ + pi+  -->  D*+ + rho+
        if( KF_D == 411 .AND. KF_M == 211 )then
            i_channel = 385
            KF_out_A = 413
            KF_out_B = 213
!       386. D- + pi-  -->  D*- + rho-
        else if( KF_D == -411 .AND. KF_M == -211 )then
            i_channel = 386
            KF_out_A = -413
            KF_out_B = -213
!       387. D*+ + pi+  -->  D+ + rho+
        else if( KF_D == 413 .AND. KF_M == 211 )then
            i_channel = 387
            KF_out_A = 411
            KF_out_B = 213
!       388. D*- + pi-  -->  D- + rho-
        else if( KF_D == -413 .AND. KF_M == -211 )then
            i_channel = 388
            KF_out_A = -411
            KF_out_B = -213
!       389. D+ + rho+  -->  D*+ + pi+
        else if( KF_D == 411 .AND. KF_M == 213 )then
            i_channel = 389
            KF_out_A = 413
            KF_out_B = 211
!       390. D- + rho-  -->  D*- + pi-
        else if( KF_D == -411 .AND. KF_M == -213 )then
            i_channel = 390
            KF_out_A = -413
            KF_out_B = -211
!       391. D*+ + rho+  -->  D+ + pi+
        else if( KF_D == 413 .AND. KF_M == 213 )then
            i_channel = 391
            KF_out_A = 411
            KF_out_B = 211
!       392. D*- + rho-  -->  D- + pi-
        else if( KF_D == -413 .AND. KF_M == -213 )then
            i_channel = 392
            KF_out_A = -411
            KF_out_B = -211
!       393. D0 + pi-  -->  D*0 + rho-
        else if( KF_D == 421 .AND. KF_M == -211 )then
            i_channel = 393
            KF_out_A = 423
            KF_out_B = -213
!       394. Dbar0 + pi+  -->  Dbar*0 + rho+
        else if( KF_D == -421 .AND. KF_M == 211 )then
            i_channel = 394
            KF_out_A = -423
            KF_out_B = 213
!       395. D*0 + pi-  -->  D0 + rho-
        else if( KF_D == 423 .AND. KF_M == -211 )then
            i_channel = 395
            KF_out_A = 421
            KF_out_B = -213
!       396. D*bar0 + pi+  -->  Dbar0 + rho+
        else if( KF_D == -423 .AND. KF_M == 211 )then
            i_channel = 396
            KF_out_A = -421
            KF_out_B = 213
!       397. D0 + rho-  -->  D*0 + pi-
        else if( KF_D == 421 .AND. KF_M == -213 )then
            i_channel = 397
            KF_out_A = 423
            KF_out_B = -211
!       398. Dbar0 + rho+  -->  Dbar*0 + pi+
        else if( KF_D == -421 .AND. KF_M == 213 )then
            i_channel = 398
            KF_out_A = -423
            KF_out_B = 211
!       399. D*0 + rho-  -->  D0 + pi-
        else if( KF_D == 423 .AND. KF_M == -213 )then
            i_channel = 399
            KF_out_A = 421
            KF_out_B = -211
!       400. D*bar0 + rho+  -->  Dbar0 + pi+
        else if( KF_D == -423 .AND. KF_M == 213 )then
            i_channel = 400
            KF_out_A = -421
            KF_out_B = 211
!       Something is wrong.
        else
            write(9,*) "No inel. channel is found in hadcas "       // &
                       "Dmeson1(), treated as elastic collision, "  // &
                       "KF_D, KF_M =", KF_D, KF_M
            return
        end if

!       Inelastic channel switch.
        if( isinel( i_channel ) == 0 ) return
!       Threshold energy.
        E_thresh = PYMASS( KF_out_A ) + PYMASS( KF_out_B )
        if( E_AB <= E_thresh ) return

!       Order exchanging corresponding to that made in "DMeson_coll".
        if( i_od == -1 )then
            KF = KF_out_A
            KF_out_A = KF_out_B
            KF_out_B = KF
        end if

!       Feedback.
        lc(icp,3) = KF_out_A
        lc(icp,4) = KF_out_B
        lc(icp,5) = i_channel
        tw(icp)   = rcsit
        ioo = 1


        return
        end



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine Dmeson2(KF_D,KF_M,E_AB,icp,ioo,i_od)
!       A part of "prod" to deal with D + pion/rho  --> ...
!       Dual production channel: 401-464, total 64 reaction channels.
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        PARAMETER(NSIZE=750000)
        common/sa10_h/csnn,cspin,cskn,cspipi,cspsn,cspsm,rcsit,ifram, &
         iabsb,iabsm,i_sigma_AQM,csspn,csspm
        common/count_h/isinel(600)
        common/ctllist_t/ lc(nsize,5), tc(nsize), tw(nsize)


        if( KF_D == 421 .AND. KF_M == 111 )then
!       401. D0 + pi0  -->  D*0 + rho0
            i_channel_1 = 401
            KF_out_A_1 = 423
            KF_out_B_1 = 113
!       402. D0 + pi0  -->  D*+ + rho-
            i_channel_2 = 402
            KF_out_A_2 = 413
            KF_out_B_2 = -213
        else if( KF_D == -421 .AND. KF_M == 111 )then
!       403. Dbar0 + pi0  -->  Dbar*0 + rho0
            i_channel_1 = 403
            KF_out_A_1 = -423
            KF_out_B_1 = 113
!       404. Dbar0 + pi0  -->  D*- + rho+
            i_channel_2 = 404
            KF_out_A_2 = -413
            KF_out_B_2 = 213
        else if( KF_D == 423 .AND. KF_M == 111 )then
!       405. D*0 + pi0  -->  D0 + rho0
            i_channel_1 = 405
            KF_out_A_1 = 421
            KF_out_B_1 = 113
!       406. D*0 + pi0  -->  D+ + rho-
            i_channel_2 = 406
            KF_out_A_2 = 411
            KF_out_B_2 = -213
        else if( KF_D == -423 .AND. KF_M == 111 )then
!       407. D*bar0 + pi0  -->  Dbar0 + rho0
            i_channel_1 = 407
            KF_out_A_1 = -421
            KF_out_B_1 = 113
!       408. D*bar0 + pi0  -->  D- + rho+
            i_channel_2 = 408
            KF_out_A_2 = -411
            KF_out_B_2 = 213
        else if( KF_D == 421 .AND. KF_M == 113 )then
!       409. D0 + rho0  -->  D*0 + pi0
            i_channel_1 = 409
            KF_out_A_1 = 423
            KF_out_B_1 = 111
!       410. D0 + rho0  -->  D*+ + pi-
            i_channel_2 = 410
            KF_out_A_2 = 413
            KF_out_B_2 = -211
        else if( KF_D == -421 .AND. KF_M == 113 )then
!       411. Dbar0 + rho0  -->  Dbar*0 + pi0
            i_channel_1 = 411
            KF_out_A_1 = -423
            KF_out_B_1 = 111
!       412. Dbar0 + rho0  -->  D*- + pi+
            i_channel_2 = 412
            KF_out_A_2 = -413
            KF_out_B_2 = 211
        else if( KF_D == 423 .AND. KF_M == 113 )then
!       413. D*0 + rho0  -->  D0 + pi0
            i_channel_1 = 413
            KF_out_A_1 = 421
            KF_out_B_1 = 111
!       414. D*0 + rho0  -->  D+ + pi-
            i_channel_2 = 414
            KF_out_A_2 = 411
            KF_out_B_2 = -211
        else if( KF_D == -423 .AND. KF_M == 113 )then
!       415. D*bar0 + rho0  -->  Dbar0 + pi0
            i_channel_1 = 415
            KF_out_A_1 = -421
            KF_out_B_1 = 111
!       416. D*bar0 + rho0  -->  D- + pi+
            i_channel_2 = 416
            KF_out_A_2 = -411
            KF_out_B_2 = 211
        else if( KF_D == 411 .AND. KF_M == -211 )then
!       417. D+ + pi-  -->  D*+ + rho-
            i_channel_1 = 417
            KF_out_A_1 = 413
            KF_out_B_1 = -213
!       418. D+ + pi-  -->  D*0 + rho0
            i_channel_2 = 418
            KF_out_A_2 = 423
            KF_out_B_2 = 113
        else if( KF_D == -411 .AND. KF_M == 211 )then
!       419. D- + pi+  -->  D*- + rho+
            i_channel_1 = 419
            KF_out_A_1 = -413
            KF_out_B_1 = 213
!       420. D- + pi+  -->  D*bar0 + rho0
            i_channel_2 = 420
            KF_out_A_2 = -423
            KF_out_B_2 = 113
        else if( KF_D == 413 .AND. KF_M == -211 )then
!       421. D*+ + pi-  -->  D+ + rho-
            i_channel_1 = 421
            KF_out_A_1 = 411
            KF_out_B_1 = -213
!       422. D*+ + pi-  -->  D0 + rho0
            i_channel_2 = 422
            KF_out_A_2 = 421
            KF_out_B_2 = 113
        else if( KF_D == -413 .AND. KF_M == 211 )then
!       423. D*- + pi+  -->  D- + rho+
            i_channel_1 = 423
            KF_out_A_1 = -411
            KF_out_B_1 = 213
!       424. D*- + pi+  -->  Dbar0 + rho0
            i_channel_2 = 424
            KF_out_A_2 = -421
            KF_out_B_2 = 113
        else if( KF_D == 411 .AND. KF_M == -213 )then
!       425. D+ + rho-  -->  D*+ + pi-
            i_channel_1 = 425
            KF_out_A_1 = 413
            KF_out_B_1 = -211
!       426. D+ + rho-  -->  D*0 + pi0
            i_channel_2 = 426
            KF_out_A_2 = 423
            KF_out_B_2 = 111
        else if( KF_D == -411 .AND. KF_M == 213 )then
!       427. D- + rho+  -->  D*- + pi+
            i_channel_1 = 427
            KF_out_A_1 = -413
            KF_out_B_1 = 211
!       428. D- + rho+  -->  D*bar0 + pi0
            i_channel_2 = 428
            KF_out_A_2 = -423
            KF_out_B_2 = 111
        else if( KF_D == 413 .AND. KF_M == -213 )then
!       429. D*+ + rho-  -->  D+ + pi-
            i_channel_1 = 429
            KF_out_A_1 = 411
            KF_out_B_1 = -211
!       430. D*+ + rho-  -->  D0 + pi0
            i_channel_2 = 430
            KF_out_A_2 = 421
            KF_out_B_2 = 111
        else if( KF_D == -413 .AND. KF_M == 213 )then
!       431. D*- + rho+  -->  D- + pi+
            i_channel_1 = 431
            KF_out_A_1 = -411
            KF_out_B_1 = 211
!       432. D*- + rho+  -->  Dbar0 + pi0
            i_channel_2 = 432
            KF_out_A_2 = -421
            KF_out_B_2 = 111
        else if( KF_D == 411 .AND. KF_M == 111 )then
!       433. D+ + pi0  -->  D*+ + rho0
            i_channel_1 = 433
            KF_out_A_1 = 413
            KF_out_B_1 = 113
!       434. D+ + pi0  -->  D*0 + rho+
            i_channel_2 = 434
            KF_out_A_2 = 423
            KF_out_B_2 = 213
        else if( KF_D == -411 .AND. KF_M == 111 )then
!       435. D- + pi0  -->  D*- + rho0
            i_channel_1 = 435
            KF_out_A_1 = -413
            KF_out_B_1 = 113
!       436. D- + pi0  -->  D*bar0 + rho-
            i_channel_2 = 436
            KF_out_A_2 = -423
            KF_out_B_2 = -213
        else if( KF_D == 413 .AND. KF_M == 111 )then
!       437. D*+ + pi0  -->  D+ + rho0
            i_channel_1 = 437
            KF_out_A_1 = 411
            KF_out_B_1 = 113
!       438. D*+ + pi0  -->  D0 + rho+
            i_channel_2 = 438
            KF_out_A_2 = 421
            KF_out_B_2 = 213
        else if( KF_D == -413 .AND. KF_M == 111 )then
!       439. D*- + pi0  -->  D- + rho0
            i_channel_1 = 439
            KF_out_A_1 = -411
            KF_out_B_1 = 113
!       440. D*- + pi0  -->  Dbar0 + rho-
            i_channel_2 = 440
            KF_out_A_2 = -421
            KF_out_B_2 = -213
        else if( KF_D == 411 .AND. KF_M == 113 )then
!       441. D+ + rho0  -->  D*+ + pi0
            i_channel_1 = 441
            KF_out_A_1 = 413
            KF_out_B_1 = 111
!       442. D+ + rho0  -->  D*0 + pi+
            i_channel_2 = 442
            KF_out_A_2 = 423
            KF_out_B_2 = 211
        else if( KF_D == -411 .AND. KF_M == 113 )then
!       443. D- + rho0  -->  D*- + pi0
            i_channel_1 = 443
            KF_out_A_1 = -413
            KF_out_B_1 = 111
!       444. D- + rho0  -->  D*bar0 + pi-
            i_channel_2 = 444
            KF_out_A_2 = -423
            KF_out_B_2 = -211
        else if( KF_D == 413 .AND. KF_M == 113 )then
!       445. D*+ + rho0  -->  D+ + pi0
            i_channel_1 = 445
            KF_out_A_1 = 411
            KF_out_B_1 = 111
!       446. D*+ + rho0  -->  D0 + pi+
            i_channel_2 = 446
            KF_out_A_2 = 421
            KF_out_B_2 = 211
        else if( KF_D == -413 .AND. KF_M == 113 )then
!       447. D*- + rho0  -->  D- + pi0
            i_channel_1 = 447
            KF_out_A_1 = -411
            KF_out_B_1 = 111
!       448. D*- + rho0  -->  Dbar0 + pi-
            i_channel_2 = 448
            KF_out_A_2 = -421
            KF_out_B_2 = -211
        else if( KF_D == 421 .AND. KF_M == 211 )then
!       449. D0 + pi+  -->  D*0 + rho+
            i_channel_1 = 449
            KF_out_A_1 = 423
            KF_out_B_1 = 213
!       450. D0 + pi+  -->  D*+ + rho0
            i_channel_2 = 450
            KF_out_A_2 = 413
            KF_out_B_2 = 113
        else if( KF_D == -421 .AND. KF_M == -211 )then
!       451. Dbar0 + pi-  -->  Dbar*0 + rho-
            i_channel_1 = 451
            KF_out_A_1 = -423
            KF_out_B_1 = -213
!       452. Dbar0 + pi-  -->  D*- + rho0
            i_channel_2 = 452
            KF_out_A_2 = -413
            KF_out_B_2 = 113
        else if( KF_D == 423 .AND. KF_M == 211 )then
!       453. D*0 + pi+  -->  D0 + rho+
            i_channel_1 = 453
            KF_out_A_1 = 421
            KF_out_B_1 = 213
!       454. D*0 + pi+  -->  D+ + rho0
            i_channel_2 = 454
            KF_out_A_2 = 411
            KF_out_B_2 = 113
        else if( KF_D == -423 .AND. KF_M == -211 )then
!       455. D*bar0 + pi-  -->  Dbar0 + rho-
            i_channel_1 = 455
            KF_out_A_1 = -421
            KF_out_B_1 = -213
!       456. D*bar0 + pi-  -->  D- + rho0
            i_channel_2 = 456
            KF_out_A_2 = -411
            KF_out_B_2 = 113
        else if( KF_D == 421 .AND. KF_M == 213 )then
!       457. D0 + rho+  -->  D*0 + pi+
            i_channel_1 = 457
            KF_out_A_1 = 423
            KF_out_B_1 = 211
!       458. D0 + rho+  -->  D*+ + pi0
            i_channel_2 = 458
            KF_out_A_2 = 413
            KF_out_B_2 = 111
        else if( KF_D == -421 .AND. KF_M == -213 )then
!       459. Dbar0 + rho-  -->  Dbar*0 + pi-
            i_channel_1 = 459
            KF_out_A_1 = -423
            KF_out_B_1 = -211
!       460. Dbar0 + rho-  -->  D*- + pi0
            i_channel_2 = 460
            KF_out_A_2 = -413
            KF_out_B_2 = 111
        else if( KF_D == 423 .AND. KF_M == 213 )then
!       461. D*0 + rho+  -->  D0 + pi+
            i_channel_1 = 461
            KF_out_A_1 = 421
            KF_out_B_1 = 211
!       462. D*0 + rho+  -->  D+ + pi0
            i_channel_2 = 462
            KF_out_A_2 = 411
            KF_out_B_2 = 111
        else if( KF_D == -423 .AND. KF_M == -213 )then
!       463. D*bar0 + rho-  -->  Dbar0 + pi-
            i_channel_1 = 463
            KF_out_A_1 = -421
            KF_out_B_1 = -211
!       464. D*bar0 + rho-  -->  D- + pi0
            i_channel_2 = 464
            KF_out_A_2 = -411
            KF_out_B_2 = 111
!       Something is wrong.
        else
            write(9,*) "No inel. channel is found in hadcas "       // &
                       "Dmeson2(), treated as elastic collision, "  // &
                       "KF_D, KF_M =", KF_D, KF_M
            return
        end if

        ratio_sig1 = 0D0
        ratio_sig2 = 0D0

        E_thresh_1 = PYMASS( KF_out_A_1 ) + PYMASS( KF_out_B_1 )
        E_thresh_2 = PYMASS( KF_out_A_2 ) + PYMASS( KF_out_B_2 )
        if( isinel( i_channel_1 ) == 1 .AND. E_AB > E_thresh_1 ) &
            ratio_sig1 = 0.5D0
        if( isinel( i_channel_2 ) == 1 .AND. E_AB > E_thresh_2 ) &
            ratio_sig2 = 0.5D0

        if( ratio_sig1 < 1D-15 .AND.ratio_sig2 < 1D-15 ) return

        ratio_sig = ratio_sig1 + ratio_sig2
        ratio1 = ratio_sig1 / ratio_sig
        if( PYR(1) <= ratio1 )then
            KF_out_A  = KF_out_A_1
            KF_out_B  = KF_out_B_1
            i_channel = i_channel_1
        else
            KF_out_A  = KF_out_A_2
            KF_out_B  = KF_out_B_2
            i_channel = i_channel_2
        end if

!       Order exchanging corresponding to that made in "DMeson_coll".
        if( i_od == -1 )then
            KF = KF_out_A
            KF_out_A = KF_out_B
            KF_out_B = KF
        end if

!       Feedback.
        lc(icp,3) = KF_out_A
        lc(icp,4) = KF_out_B
        lc(icp,5) = i_channel
        tw(icp)   = rcsit
        ioo = 1


        return
        end



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine Dpn_coll(KF_D,KF_N,E_AB,icp,ioo)
!       A part of "prod" to deal with D + p/n  --> ...
!       Single production channel: 465-488, total 32 reaction channels.
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        PARAMETER(NSIZE=750000)
        common/sa10_h/csnn,cspin,cskn,cspipi,cspsn,cspsm,rcsit,ifram, &
         iabsb,iabsm,i_sigma_AQM,csspn,csspm
        common/count_h/isinel(600)
        common/ctllist_t/ lc(nsize,5), tc(nsize), tw(nsize)


!       Make order as D + N.
        i_order = 1
!       Exchange orders if first one is not D meson.
        if( KF_D > 999 )then
            KF   = KF_D
            KF_D = KF_N
            KF_N = KF
            i_order = -1
        end if

!       465. D+  + p+    -->  D*+ + p+
        if( KF_D == 411 .AND. KF_N == 2212 )then
            i_channel = 465
            KF_out_A = 413
            KF_out_B = 2212
!       466. D+  + p-    -->  D*+ + p-
        else if( KF_D == 411 .AND. KF_N == -2212 )then
            i_channel = 466
            KF_out_A = 413
            KF_out_B = -2212
!       467. D+  + n     -->  D*+ + n
        else if( KF_D == 411 .AND. KF_N == 2112 )then
            i_channel = 467
            KF_out_A = 413
            KF_out_B = 2112
!       468. D+  + nbar  -->  D*+ + nbar
        else if( KF_D == 411 .AND. KF_N == -2112 )then
            i_channel = 468
            KF_out_A = 413
            KF_out_B = -2112
!       469. D-  + p+    -->  D*- + p+
        else if( KF_D == -411 .AND. KF_N == 2212 )then
            i_channel = 469
            KF_out_A = -413
            KF_out_B = 2212
!       470. D-  + p-    -->  D*- + p-
        else if( KF_D == -411 .AND. KF_N == -2212 )then
            i_channel = 470
            KF_out_A = -413
            KF_out_B = -2212
!       471. D-  + n     -->  D*- + n
        else if( KF_D == -411 .AND. KF_N == 2112 )then
            i_channel = 471
            KF_out_A = -413
            KF_out_B = 2112
!       472. D-  + nbar  -->  D*- + nbar
        else if( KF_D == -411 .AND. KF_N == -2112 )then
            i_channel = 472
            KF_out_A = -413
            KF_out_B = -2112
!       473. D*+ + p+    -->  D+  + p+
        else if( KF_D == 413 .AND. KF_N == 2212 )then
            i_channel = 473
            KF_out_A = 411
            KF_out_B = 2212
!       474. D*+ + p-    -->  D+  + p-
        else if( KF_D == 413 .AND. KF_N == -2212 )then
            i_channel = 474
            KF_out_A = 411
            KF_out_B = -2212
!       475. D*+ + n     -->  D+  + n
        else if( KF_D == 413 .AND. KF_N == 2112 )then
            i_channel = 475
            KF_out_A = 411
            KF_out_B = 2112
!       476. D*+ + nbar  -->  D+  + nbar
        else if( KF_D == 413 .AND. KF_N == -2112 )then
            i_channel = 476
            KF_out_A = 411
            KF_out_B = -2112
!       477. D*- + p+    -->  D-  + p+
        else if( KF_D == -413 .AND. KF_N == 2212 )then
            i_channel = 477
            KF_out_A = -411
            KF_out_B = 2212
!       478. D*- + p-    -->  D-  + p-
        else if( KF_D == -413 .AND. KF_N == -2212 )then
            i_channel = 478
            KF_out_A = -411
            KF_out_B = -2212
!       479. D*- + n     -->  D-  + n
        else if( KF_D == -413 .AND. KF_N == 2112 )then
            i_channel = 479
            KF_out_A = -411
            KF_out_B = 2112
!       480. D*- + nbar  -->  D-  + nbar
        else if( KF_D == -413 .AND. KF_N == -2112 )then
            i_channel = 480
            KF_out_A = -411
            KF_out_B = -2112
!       481. D0     + p+    -->  D*0 + p+
        else if( KF_D == 421 .AND. KF_N == 2212 )then
            i_channel = 481
            KF_out_A = 423
            KF_out_B = 2212
!       482. D0     + p-    -->  D*0 + p-
        else if( KF_D == 421 .AND. KF_N == -2212 )then
            i_channel = 482
            KF_out_A = 423
            KF_out_B = -2212
!       483. D0     + n     -->  D*0 + n
        else if( KF_D == 421 .AND. KF_N == 2112 )then
            i_channel = 483
            KF_out_A = 423
            KF_out_B = 2112
!       484. D0     + nbar  -->  D*0 + nbar
        else if( KF_D == 421 .AND. KF_N == -2112 )then
            i_channel = 484
            KF_out_A = 423
            KF_out_B = -2112
!       485. Dbar0  + p+    -->  D*bar0 + p+
        else if( KF_D == -421 .AND. KF_N == 2212 )then
            i_channel = 485
            KF_out_A = -423
            KF_out_B = 2212
!       486. Dbar0  + p-    -->  D*bar0 + p-
        else if( KF_D == -421 .AND. KF_N == -2212 )then
            i_channel = 486
            KF_out_A = -423
            KF_out_B = -2212
!       487. Dbar0  + n     -->  D*bar0 + n
        else if( KF_D == -421 .AND. KF_N == 2112 )then
            i_channel = 487
            KF_out_A = -423
            KF_out_B = 2112
!       488. Dbar0  + nbar  -->  D*bar0 + nbar
        else if( KF_D == -421 .AND. KF_N == -2112 )then
            i_channel = 488
            KF_out_A = -423
            KF_out_B = -2112
!       489. D*0    + p+    -->  D0 + p+
        else if( KF_D == 423 .AND. KF_N == 2212 )then
            i_channel = 489
            KF_out_A = 421
            KF_out_B = 2212
!       490. D*0    + p-    -->  D0 + p-
        else if( KF_D == 423 .AND. KF_N == -2212 )then
            i_channel = 490
            KF_out_A = 421
            KF_out_B = -2212
!       491. D*0    + n     -->  D0 + n
        else if( KF_D == 423 .AND. KF_N == 2112 )then
            i_channel = 491
            KF_out_A = 421
            KF_out_B = 2112
!       492. D*0    + nbar  -->  D0 + nbar
        else if( KF_D == 423 .AND. KF_N == -2112 )then
            i_channel = 492
            KF_out_A = 421
            KF_out_B = -2112
!       493. D*bar0 + p+    -->  Dbar0 + p+
        else if( KF_D == -423 .AND. KF_N == 2212 )then
            i_channel = 493
            KF_out_A = -421
            KF_out_B = 2212
!       494. D*bar0 + p-    -->  Dbar0 + p-
        else if( KF_D == -423 .AND. KF_N == -2212 )then
            i_channel = 494
            KF_out_A = -421
            KF_out_B = -2212
!       495. D*bar0 + n     -->  Dbar0 + n
        else if( KF_D == -423 .AND. KF_N == 2112 )then
            i_channel = 495
            KF_out_A = -421
            KF_out_B = 2112
!       496. D*bar0 + nbar  -->  Dbar0 + nbar
        else if( KF_D == -423 .AND. KF_N == -2112 )then
            i_channel = 496
            KF_out_A = -421
            KF_out_B = -2112
!       Something is wrong.
        else
            write(9,*) "No inel. channel is found in hadcas "        // &
                       "Dpn_coll(), treated as elastic collision, "  // &
                       "KF_D, KF_M =", KF_D, KF_N
            return
        end if

!       Inelastic channel switch.
        if( isinel( i_channel ) == 0 ) return
!       Threshold energy.
        E_thresh = PYMASS( KF_out_A ) + PYMASS( KF_out_B )
        if( E_AB <= E_thresh ) return

!       Order exchanging corresponding to that made in "DMeson_coll".
        if( i_order == -1 )then
            KF = KF_out_A
            KF_out_A = KF_out_B
            KF_out_B = KF
        end if

!       Feedback.
        lc(icp,3) = KF_out_A
        lc(icp,4) = KF_out_B
        lc(icp,5) = i_channel
        tw(icp)   = rcsit
        ioo = 1


        return
        end



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine decay(time,deltt)
!!      Deals with particle decay in transport processes.
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        PARAMETER (KSZJ=80000,NSIZE=750000)
        common/sa1_h/nsa,non1,ksa(kszj,5),psa(kszj,5),vsa(kszj,5)
        common/sa8_h/tau(kszj),ishp(kszj)
        common/sa14/ipyth(2000),idec(2000),iwide(2000)
        common/sa17/nde,non17,kde(10,5),pde(10,5),vde(10,5)
        common/sa19_h/coor(3)
        common/ctllist_h/nctl,noinel(600),nctl0,noel
        common/ctllist_t/ lc(nsize,5), tc(nsize), tw(nsize)
        dimension pp(5),rr(4)
!       idec: stord the line number (in the particle list) of
!        particles after decay
!       iwide: stord the line number (in the particle list) of decaying
!        particles
!       the messages of decayed particles are
!        stored in the varibles and arrays in 'sa17'
!       pp, rr: momentum, position of decaying particle


        nde = 0
        kde = 0
        pde = 0D0
        vde = 0D0
        iwide =0
        ii = 0
        do i=1,nsa,1
            kf = ksa(i,2)
            if( ABS(kf) == 213 .OR. kf == 113 )then
                ii = ii + 1
                iwide(ii) = i
            end if
        end do
        if(ii == 0) return
        do 100 i=1,ii,1
!       deal with decay of decaying particles one by one
        jj = iwide(i)
        kf = ksa(jj,2)
        do i1=1,5,1
            pp(i1) = psa(jj,i1)
        end do
        do i1=1,3,1
            rr(i1) = vsa(jj,i1)
        end do
        rr(4) = time
        ee = pp(4)
        if(ee < 1D-15) ee = 1D-15
        damass = pp(5)
!       compo = damass * 0.151D0 * deltt / 0.197327D0 / ee
        tauj = tau(jj)
        tdel = deltt
        if( (time - deltt) < tauj ) tdel = time - tauj
        compo = damass * 0.151D0 * tdel / 0.197327D0 / ee
        if(compo < 1D-15) compo=1D-15
        prob = 1D0 - EXP( -compo )
!       0.151: full width in [GeV] of decaying particle (rho)
!       mass, full width, energy: in [GeV], time in [GeV^-1] originally
!       0.197327D0: conversion factor from [GeV^-1] to [fm/c]
        if( PYR(1) <= prob )then
!       perform a particle decay
            call decpr(pp,rr,kf)
!       it gives momentum and position to the particles after decay
            call upddep(jj,time)
!       update particle list after a particle decay
            call upddet(time)
!       update collision time list after a particle decay
!           if(nctl == 0)return
        end if
100     continue


        return
        end



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine decpr(pp,rr,kf)
!!      Gives momentum and position to decayed particles.
!       pp, rr: momentum, position of decaying particle
!       kf: flavour code of decaying particle
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        INTEGER PYCOMP
        PARAMETER (KSZJ=80000)
        COMMON/PYDAT2/KCHG(500,4),PMAS(500,4),PARF(2000),VCKM(4,4)
        common/sa8_h/tau(kszj),ishp(kszj)
        common/sa14/ipyth(2000),idec(2000),iwide(2000)
        common/sa17/nde,non17,kde(10,5),pde(10,5),vde(10,5)
        dimension pp(5),rr(4),pi(4),pj(4),rd(3),b(3)


        amas = pp(5)
        amas2 = amas * amas
        nde = 2
!       assume two body decay

        if(kf == 213)then
!       rho+ to pion+ + pion0
            am1 = pmas( pycomp(211), 1 )
            am2 = pmas( pycomp(111), 1 )
            kde(1,2) = 211
            kde(2,2) = 111
        end if
        if(kf == -213)then
!       rho- to pion- + pion0
        am1=pmas(pycomp(-211),1)
        am2=pmas(pycomp(111),1)
        kde(1,2)=-211
        kde(2,2)=111
        endif
        if(kf == 113)then
!       rho0 to pion0 + pion0
        am1=pmas(pycomp(111),1)
        am2=am1
        kde(1,2)=111
        kde(2,2)=111
        endif

        agum=(amas2-(am1+am2)**2)*(amas2-(am1-am2)**2)
        if(agum <= 0.)then
        write(9,*) 'kf,kf1,kf2=',kf,kde(1,2),kde(2,2)
        write(9,*) 'amas,am1,am2=',amas,am1,am2
        agum=1.e-15
        endif

        pab=sqrt(agum)/2./amas
        cita=2.*pyr(1)-1.
        if(cita > 0.9999)cita=0.9999
        sita=sqrt(1.-cita*cita)
        fi=2.*3.1416*pyr(1)
        px=pab*sita*cos(fi)
        py=pab*sita*sin(fi)
        pz=pab*cita
        pi(1)=px
        pi(2)=py
        pi(3)=pz
        pi(4)=sqrt(pab*pab+am1*am1)
        pj(1)=-px
        pj(2)=-py
        pj(3)=-pz
        pj(4)=sqrt(pab*pab+am2*am2)
        ee=pp(4)
        do i=1,3
        b(i)=pp(i)/ee
        enddo
        ilo=1
        call lorntz(ilo,b,pi,pj)
        do i=1,4
        pde(1,i)=pi(i)
        pde(2,i)=pj(i)
        enddo
        pde(1,5)=am1
        pde(2,5)=am2

        rrp=0.8
!       rrp: the radius (in fm) of decaying particle
        do i=1,nde
        cita=2.*pyr(1)-1.
        if(cita > 0.9999)cita=0.9999
        sita=sqrt(1.-cita*cita)
        fi=2.*3.1416*pyr(1)
        rd(1)=rrp*sita*cos(fi)
        rd(2)=rrp*sita*sin(fi)
        rd(3)=rrp*cita
!       arrange the particle i after decay on the surface
!        of sphere with radius rrp and centered at the position of parent
        do j=1,3
        vde(i,j)=rd(j)+rr(j)
        enddo
        vde(i,4)=rr(4)
        enddo


        return
        end



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine upddep(ij,time)
!!      update particle list after a particle decay and
!!       truncate collision list correspondingly.
!       ij: line number (in particle list 'sa1_h') of decaying particle
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        PARAMETER (KSZJ=80000,NSIZE=750000)
        common/sa1_h/nsa,non1,ksa(kszj,5),psa(kszj,5),vsa(kszj,5)
        common/sa8_h/tau(kszj),ishp(kszj)
        common/sa9_h/kfmax,kfaco(100),numb(100),non9,disbe(100,100)
        common/sa14/ipyth(2000),idec(2000),iwide(2000)
        common/sa17/nde,non17,kde(10,5),pde(10,5),vde(10,5)
        common/sa19_h/coor(3)
        common/sa20_h/t0,cspipiKK,dep,ddt,edipi,epin,ecsnn,ekn,ecspsn,ecspsm, &
         rnt,rnp,ecsspn,ecsspm
        common/ctllist_h/nctl,noinel(600),nctl0,noel
        common/ctllist_t/ lc(nsize,5), tc(nsize), tw(nsize)
!       idec: stord the order number (in the particle list) of
!        particles after decay


        do m=1,2000
        idec(m)=0
        enddo
!       put decayed particles (in 'sa17') into particle list 'sa1_h'
!        and turncate collision list correspondingly.
        ll=ij
        kd=ksa(ll,2)
        do 500 i=1,nde
        kf=kde(i,2)
        do 600 j=1,kfmax
        if(kf /= kfaco(j)) goto 600
        jj=numb(j)+1
!       update the particle list.
        do m=nsa,jj,-1
        mm=m+1
        ksa(mm,2)=ksa(m,2)
        ksa(mm,1)=1
        ksa(mm,3)=ksa(m,3)
        do m1=1,5
        psa(mm,m1)=psa(m,m1)
        vsa(mm,m1)=vsa(m,m1)
        enddo
        ishp(mm)=ishp(m)
        tau(mm)=tau(m)
        enddo
        if(ll >= jj)ll=ll+1
!       give proper values to particle jj.
        ksa(jj,2)=kf
        ksa(jj,1)=1
        ksa(jj,3)=0
        do m=1,4
        psa(jj,m)=pde(i,m)
        vsa(jj,m)=vde(i,m)
        enddo
        psa(jj,5)=pde(i,5)   ! 220924 Lei
        ishp(jj)=1
        tau(jj)=time
!       the values of 'ishp' and 'tau' of decayed particles
!        are given here, its formation time are assumed to be zero
        do m=j,kfmax
        numb(m)=numb(m)+1
        enddo
        idec(i)=jj
        do m=1,2000
        ipym=idec(m)
        iwi=iwide(m)
        if(ipym > jj)idec(m)=ipym+1
        if(iwi > jj)iwide(m)=iwi+1
        enddo
!       update the values of lc(m,1-2) with value >= jj
        do m=1,nctl
        lc1=lc(m,1)
        if(lc1 >= jj)lc(m,1)=lc1+1
        lc2=lc(m,2)
        if(lc2 >= jj)lc(m,2)=lc2+1
        enddo
        goto 200
600     continue
200     continue
        nsa=nsa+1
500     continue

!       Throws away the colli. pairs with partner of ll.
!        (decaying particle)
        j=0
        do i=1,nctl
        i1=lc(i,1)
        j1=lc(i,2)
        if(i1 == ll) goto 401
        if(j1 == ll) goto 401
        ! 061123 Lei Avoid the particle collideing with itself.
        if(i1 == j1) goto 401
        if((tc(i)-time) <= ddt) goto 401
!       Throws away the pairs with tc <= time.
        j=j+1
        tc(j)=tc(i)
        tw(j)=tw(i)
        do m=1,5
        lc(j,m)=lc(i,m)
        enddo
401     continue
        enddo
        do i=j+1,nctl+1
        tc(i)=0.0
        tw(i)=0.0
        do m=1,5
        lc(i,m)=0
        enddo
        enddo
        nctl=j

!       remove decaying particle (ll) from particle list and
!        truncate the collision list correspondingly.
        kf=ksa(ll,2)
!       kf=kf1
!       ll=l
        if(ll == nsa)then
        kfd=numb(kfmax)-numb(kfmax-1)

!       following statements are added by sa on 07/March/97
        i2 = 0
        if(kfd == 0)then
        numbm=numb(kfmax)
        do i1=1,kfmax
        if(numb(i1) /= numbm) goto 3000
        i2=i1
        goto 3001
3000    enddo
3001    continue
        if(i2 > 0)then
        do i1=i2,kfmax
        numb(i1)=numb(i1)-1
        enddo
        end if
        goto 400
        endif
        numb(kfmax)=numb(kfmax)-1
400     goto 100
        endif
        do j=ll+1,nsa
        jj=j-1
        ksa(jj,2)=ksa(j,2)
        ksa(jj,1)=1
        ksa(jj,3)=ksa(j,3)
        do m=1,5
        psa(jj,m)=psa(j,m)
        vsa(jj,m)=vsa(j,m)
        enddo
        ishp(jj)=ishp(j)
        tau(jj)=tau(j)
        enddo
        do m=1,nctl
        lc1=lc(m,1)
        lc2=lc(m,2)
        if(lc1 > ll)lc(m,1)=lc1-1
        if(lc2 > ll)lc(m,2)=lc2-1
        enddo
        do 800 j=1,kfmax
        if(kf /= kfaco(j)) goto 800
        do m=j,kfmax
        numb(m)=numb(m)-1
        enddo
        goto 100
800     continue
100     continue
        nsa=nsa-1
        do m=1,2000
        ipym=idec(m)
        iwi=iwide(m)
        if(ipym > ll)idec(m)=ipym-1
        if(iwi > ll)iwide(m)=iwi-1
        enddo


        return
        end



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine upddet(time)
!!      Updates collision list after a particle decay.
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        PARAMETER (NSIZE=750000)
        common/sa9_h/kfmax,kfaco(100),numb(100),non9,disbe(100,100)
        common/sa10_h/csnn,cspin,cskn,cspipi,cspsn,cspsm,rcsit,ifram, &
         iabsb,iabsm,i_sigma_AQM,csspn,csspm
        common/sa14/ipyth(2000),idec(2000),iwide(2000)
        common/sa17/nde,non17,kde(10,5),pde(10,5),vde(10,5)
        common/sa24/adj1(40),nnstop,non24,zstop
        common/ctllist_h/nctl,noinel(600),nctl0,noel
        common/ctllist_t/ lc(nsize,5), tc(nsize), tw(nsize)
!       idec: store line number (in the particle list) of decayed particles


        dddt = adj1(11)

        nctl=nctl+1
!       loop over particle list

!       Considers channels that are not well defined and treats them as elastic.
        m_upp = numb(kfmax)

        j21 = idec(1)
        j22 = idec(2)
        do j11=1,nde,1
            j1 = idec(j11)
            if( j1 > m_upp ) cycle
            do i=1,m_upp,1
!       decayed particles could not collide with each other immediately
                if(j11 == 1 .and. i == j22) cycle
                if(j11 == 2 .and. i == j21) cycle
!       Avoids the particle collideing with itself.
                if( i == j11 ) cycle
                i1 = i
                call rsfilt_h(j1,i1,iflag)
                if(iflag == 0) cycle
                tc(nctl) = 0D0
                tw(nctl) = 0D0
                call tcolij_h(i1,j1,time,nctl)
!               if( tc(nctl) > 1D-7 ) nctl = nctl + 1
                tci = tc(nctl)
                if(tci > 1D-7)then
                    do jj = 1, nctl-1, 1
                        tcj = tc(jj)
                        if( ABS(tcj - tci) < dddt )then
                            nctl = nctl + 1
                            exit
                        end if
                    end do
                    nctl = nctl + 1
                end if
            end do
        end do
        nctl = nctl - 1


        return
        end



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        function sigma_AQM( KF1_IN, KF2_IN, sigma_NN_IN )
!!      Calculates hh cross-section via AQM model. Unit in mb.
!!      Effective valence quark number:
!!          n_eff = n_d + n_u + 0.6*n_s + 0.2*n_c + 0.07*n_b.
!!      Cross section of hh collision:
!!          sigma_AQM = sigma_NN * n_eff_h1/3 * n_eff_h2/3.
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        real*8 n_eff(2)
        dimension n_q(6)


        KF1 = KF1_IN
        KF2 = KF2_IN
        sigma_NN = sigma_NN_IN
!       Temporary treatment.
        if( KF1 == 100443 ) KF1 = 443
        if( KF1 == 100553 ) KF1 = 553
        if( KF1 == 10441 )  KF1 = 441
        if( KF1 == 10551 )  KF1 = 551
        if( KF1 == 20443 )  KF1 = 443
        if( KF1 == 20553 )  KF1 = 553
        !
        if( KF2 == 100443 ) KF2 = 443
        if( KF2 == 100553 ) KF2 = 553
        if( KF2 == 10441 )  KF2 = 441
        if( KF2 == 10551 )  KF2 = 551
        if( KF2 == 20443 )  KF2 = 443
        if( KF2 == 20553 )  KF2 = 553
        ! For h1.
        KFABS = ABS( KF1 )
        KF = KF1
        do i_KF = 1,2,1
            ! Number of valence quarks with different flavors of h1 and h2.
            n_q = 0
!       Normal meson i_KF.
            if( KFABS > 99 .AND. KFABS < 1000 )then
                KF_q1 = ABS( KFABS / 100 )
                ! Error.
                if( KF_q1 < 1 .OR. KF_q1 > 6 ) goto 100
                n_q( KF_q1 ) = 1
                KF_q2 = ABS( ( KFABS - KF_q1*100 ) / 10 )
                ! Error.
                if( KF_q2 < 1 .OR. KF_q2 > 6 ) goto 100
                n_q( KF_q2 ) = n_q( KF_q2 ) + 1
!       Normal baryon i_KF.
            else if( KFABS > 999 .AND. KFABS < 10000 )then
                KF_q1 = ABS( KFABS / 1000 )
                ! Error.
                if( KF_q1 < 1 .OR. KF_q1 > 6 ) goto 200
                n_q( KF_q1 ) = 1
                KF_q2 = ABS( ( KFABS - KF_q1*1000 ) / 100 )
                ! Error.
                if( KF_q2 < 1 .OR. KF_q2 > 6 ) goto 200
                n_q( KF_q2 ) = n_q( KF_q2 ) + 1
                KF_q3 = ABS( ( KFABS - KF_q1*1000 - KF_q2*100 ) / 10 )
                n_q( KF_q3 ) = n_q( KF_q3 ) + 1
            end if
!       Effective valence quark number.
            n_eff( i_KF ) =   1D0*n_q(1) + 1D0*n_q(2) + &
                            0.6D0*n_q(3) + &
                            0.2D0*n_q(4) + 0.07D0*n_q(5)
            ! For h2.
            KFABS = ABS( KF2 )
            KF = KF2
        end do

!       hh cross section.
        sigma_AQM = sigma_NN * n_eff(1)/3D0 * n_eff(2)/3D0

        return

100     write(22,*) "Abort! Abnormal meson in xsection_AQM, i_KF, KF1:", i_KF,KF
        STOP
200     write(22,*) "Abort! Abnormal baryon in xsection_AQM, i_KF, KF1:",i_KF,KF
        STOP


        return
        end



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        SUBROUTINE PYCIDATA
!!      Loads some default values for hadronic scatterings.
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        COMMON/PYCIDAT1/KFACOT(100),DISDET(100),ISINELT(600)
        COMMON/PYCIDAT2/KFMAXT,NONCI2,PARAM(20),WEIGH(600)
        COMMON/SA13/KJP20,NON13,VJP20,VJP21,VJP22,VJP23
        SAVE /PYCIDAT1/,/PYCIDAT2/,/SA13/
!       In order of
!               1            2          3          4            5             6             7             8         9           10
!    0     proton,     neutron,     pbar-,      nbar,         pi+,          pi-,          pi0,           K-,    Kbar0,      Sigma0,
!    1     Sigma-,      Sigma+, Sigmabar0, Sigmabar+,   Sigmabar-,      Lambda0,   Lambdabar0,           K0,       K+,         Xi-,
!    2     Xibar+,         Xi0,    Xibar0,    Omega-,   Omegabar+,       Delta-,       Delta0,       Delta+,  Delta++,        rho+,
!    3       rho-,        rho0,     J/psi,      psi',   Deltabar+,    Deltabar0,    Deltabar-,   Deltabar--,       D+,          D-,
!    4         D0,       Dbar0,       D*+,       D*-,         D*0,       D*0bar,    Lambda_c+, Lambda_cbar-,     D_s+,        D_s-,
!    5      D*_s+,       D*_s-,       K*+,       K*-,         K*0,       K*bar0,      Upsilon,     Upsilon',   chi_0c,      chi_1c,
!    6     chi_2c,    Sigma_c0,  Sigma_c+, Sigma_c++, Sigma_cbar0,  Sigma_cbar+, Sigma_cbar++,        omega,       B0,       B0bar,
!    7         B+,          B-,      B_s0,   B_sbar0,        B_c+,         B_c-,          B*0,       B*bar0,      B*+,         B*-,
!    8      B*_s0,    B*_sbar0,     B*_c+,     B*_c-,   Lambda_b0, Lambda_bbar0,     Sigma_b0,  Sigma_bbar0, Sigma_b-, Sigma_bbar+,
!    9   Sigma_b+, Sigma_bbar-,        8*0
!       (92 particles)
        INTEGER :: KFACOT0(100)
        DATA KFACOT0 / &
!           1      2      3       4      5      6      7       8      9     10
         2212,  2112, -2212,  -2112,   211,  -211,   111,   -321,  -311,  3212, &   ! 10
         3112,  3222, -3212,  -3112, -3222,  3122, -3122,    311,   321,  3312, &   ! 20
        -3312,  3322, -3322,   3334, -3334,  1114,  2114,   2214,  2224,   213, &   ! 30
         -213,   113,   443, 100443, -1114, -2114, -2214,  -2224,   411,  -411, &   ! 40
          421,  -421,   413,   -413,   423,  -423,  4122,  -4122,   431,  -431, &   ! 50
          433,  -433,   323,   -323,   313,  -313,   553, 100553, 10441, 20443, &   ! 60
          445,  4112,  4212,   4222, -4112, -4212, -4222,    223,   511,  -511, &   ! 70
          521,  -521,   531,   -531,   541,  -541,   513,   -513,   523,  -523, &   ! 80
          533,  -533,   543,   -543,  5122, -5122,  5212,  -5212,  5112, -5112, &   ! 90
         5222, -5222,   8*0       /                                                 ! 100
        ! Nucleon & Delta = 0.5 fm.
        REAL(KIND=8) :: DISDET0(100)
        DATA DISDET0/ 0.5D0, 0.5D0, 0.5D0, 0.5D0, 21*0D0, &
                      0.5D0, 0.5D0, 0.5D0, 0.5D0, 5*0D0,  &
                      0.5D0, 0.5D0, 0.5D0, 0.5D0, 62*0D0 /
        ! 112*1 for D meson.
        INTEGER :: ISINELT0(600)
        DATA ISINELT0 / 384*1, 112*1, 96*0, 8*1 /
        DATA KFMAXT0 / 92 /
        REAL(KIND=8) :: PARAM0(20)
        DATA PARAM0 /  40D0, 25D0,   21D0,   10D0, 2D0, 0.85D0,  1D0, 0.02D0, &
                     0.1D0,  4D0, 0.16D0, 0.04D0, 6D0,    3D0, 12D0, 6D0, 4*0D0/
        REAL(KIND=8) :: WEIGH0(600)
        DATA WEIGH0 / 600*1D0 /
        DATA KJP020, VJP020,VJP021,VJP022,VJP023 / 1, 0.3D0, 4D0, 1.5D0, 8D0 /


        KFACOT  = KFACOT0
        DISDET  = DISDET0
        ISINELT = ISINELT0
        KFMAXT  = KFMAXT0
        PARAM   = PARAM0
        WEIGH   = WEIGH0
        KJP20   = KJP020
        VJP20   = VJP020
        VJP21   = VJP021
        VJP22   = VJP022
        VJP23   = VJP023


        RETURN
        END
!***********************************************************************
!...............   Switches and Parameters from LUCIAE   ...............
! KFACOT: flavor order of considered particles
! DISDET: allowable minimum distance between two
!          particles,=0.5 between two necleons (Delta),=0 otherwise
! ISINELT: switch for i-th inelastic channel
!           =0 closed, =1 opened
! KFMAXT (D=92): KFMAXT kinds of particles are involved in rescattering
! PARAM(1) (D=40.0mb) total cross-section of nucleon-nucleon
! PARAM(2) (D=25.0mb) total cross-section of pi-nucleon
! PARAM(3) (D=21.0mb) total cross-section of K-nucleon
! PARAM(4) (D=10.0mb) total cross-section of pi-pi
! PARAM(5) (D=2.0mb) cross-section of pi+pi --> K Kbar
! PARAM(6) (D=0.85) ratio of inelastic cross-section to total x-section
! PARAM(7) (D=1.0fm) formation time at rest-frame of particle
! PARAM(8) (D=0.02fm) time accuracy used in hadron cascade
! PARAM(9) (D=0.1) accuracy of four-momentum conservation
! PARAM(10) (D=4.0) size of effective rescattering region is product of
!            PARAM(10) and radius of target, origin is set on (0,0,0)
! PARAM(11) (D=0.16fm^-3) nucleon density of nucleus
! PARAM(12) (D=0.04 GeV^2/c^2) The <Pt^2> for the Gaussian distribution
!            of spectator, not used anymore
! PARAM(13) (D=6.0mb) total cross-section of J/psi + n
! PARAM(14) (D=3.0mb) total cross-section of J/psi + meson
! PARAM(15) (D=12.0mb) total cross-section of psi' + n
! PARAM(16) (D=6.0mb) total cross-section of psi' + meson
!       kjp20 = 0 : energy dependent cross section
!             = 1 : constant cross section
!       vjp20 : constant cross section of strangeness production
!       vjp21 : cross section of pion + p to pion + Delta
!       vjp22 : cross section of pion + p to rho + p
!       vjp23 : cross section of n + n to n + Delta
!***********************************************************************



!ccccccccccccccccccccccccccccccccccccc end ccccccccccccccccccccccccccccccccccccc