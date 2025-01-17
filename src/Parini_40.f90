!! parini_40.f90 is a part of the PACIAE event generator.
!! Copyright (C) 2024 PACIAE Group.
!! PACIAE is licensed under the GNU GPL v2 or later, see LICENSE for details.
!! Open source: https://github.com/ArcsaberHep/PACIAE4
!! Author: Ben-Hao Sa, December 2003 - January 2025.

!> This is the initialization program to generate the initial partonic state
!!  for C-framework or the intermediate hadronic state for A- and B-framework.

!!                                             By Ben-Hao at CIAE on 04/12/2003
!!                                  Last updated by An-Ke at UiO  on 17/01/2025


        subroutine parini( time_ini, ijk )
!!      Generates the initial partonic state in the relativistic
!!       NA, AN, AA, lN, & lA collision based on 'PYTHIA' for C-framework
!!       (i_mode=3/6). Generates the intermediate hadronic state before the
!!       hadronic rescattering for A-, B- and D-frameworks (i_mode=1, 2/5, 4,
!!       respectively).
!
!       It was composed by Ben-Hao Sa on 04/12/2003.
!       Its input messages are in 'PYJETS'.
!       Its intermediate working arraies are in 'sa2'.
!       'saf' is consistent with 'sa2'.
!       'saf' to 'PYJETS' after calling 'scat'.
!       Its output messages are in 'PYJETS' (partons) and 'sbh' (hadrons)
!        for the case of mstptj=0 (C- and D-framework), but is in 'sbh' for
!        mstptj=1 (A- and B-frameworks).
!
!       The hh, ll & lh collisions would be executed in "main.f90" directly.
!       However, some historical statements has been reserved here, which give
!        the possibility to perform all kinds of collisions (NN, NA, AN, ll, lN,
!        and lA).
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        PARAMETER (KSZJ=80000)
        PARAMETER (NSIZE=750000)
        COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
        COMMON/PYDAT3/MDCY(500,3),MDME(8000,2),BRAT(8000),KFDP(8000,5)
        COMMON/PYJETS/N,NPAD,K(KSZJ,5),P(KSZJ,5),V(KSZJ,5)
!       Those variables in above common blocks are defined in "JETSET".
        COMMON/PYSUBS/MSEL,MSUB(500),KFIN(2,-40:40),NON,CKIN(200)
        COMMON/PYPARS/MSTP(200),PARP(200),MSTI(200),PARI(200)
!       Those variables in above common block are defined in "PYTHIA".
        COMMON/PYCIDAT2/KFMAXT,NONCI2,PARAM(20),WEIGH(600)
        common/sa1/kjp21,non1,bp,iii,neve,nout,nosc
        common/sa2/nsa,non2,ksa(kszj,5),psa(kszj,5),vsa(kszj,5)
        common/sa4/tau(kszj),tlco(kszj,4)
        common/sa5/kfmax,kfaco(100),numb(100),numbs(100),non5, &
         disbe(100,100)
        common/sa6/kfmaxi,nwhole
        common/sa10/csnn,cspin,cskn,cspipi,cspsn,cspsm,rcsit,ifram, &
         iabsb,iabsm,i_sigma_AQM,ajpsi,csspn,csspm,csen
        common/sa14/ipyth(2000),idec(2000),iwide(2000)
        common/sa21/pincl(5),pscal(5),pinch(5),vnu,fq2,w2l,yyl,zl,xb, &
         pph,vnlep
        common/sa24/adj1(40),nnstop,non24,zstop
        common/sa26/ndiq(kszj),npt(kszj),ifcom(kszj),idi,idio
        common/sa27/itime,kjp22,gtime,astr,akapa(6),parj1,parj2,parj3, &
         parj21,parj4,adiv,gpmax,nnc
        common/sa30/vneump,vneumt,mstptj
        common/sbe/nbe,nonbe,kbe(kszj,5),pbe(kszj,5),vbe(kszj,5)
        common/saf/naf,nonaf,kaf(kszj,5),paf(kszj,5),vaf(kszj,5)
        common/sbh/nbh,nonbh,kbh(kszj,5),pbh(kszj,5),vbh(kszj,5)
        common/sgam/ngam,nongam,kgam(kszj,5),pgam(kszj,5),vgam(kszj,5)
        common/wz/c17(500,3),ishp(kszj),tp(500),coor(3),p17(500,4)
        common/count/isinel(600)
        common/papr/t0,cspipiKK,dep,ddt,edipi,epin,ecsnn,ekn,ecspsn,ecspsm, &
         rnt,rnp,rao,rou0,vneu,vneum,ecsspn,ecsspm,ecsen
        common/syspar/ipden,itden,suppm,suptm,suppc,suptc,r0p,r0t, &
         nap,nat,nzp,nzt,pio
        common/ctllist/nctl,noinel(600),nctl0,nctlm
        common/ctllist_t/ lc(nsize,5), tc(nsize), tw(nsize)
        common/sa12/ppsa(5),nchan,nsjp,sjp,taup,taujp
        common/sa15/nps,npsi,pps(5000,5),ppsi(5000,5)
        common/sa23/kpar,knn,kpp,knp,kep,NON_p,skpar,sknn,skpp,sknp,skep
        common/sa33/smadel,ecce,secce,parecc,iparres
        common/schuds/schun,schudn,schudsn,sfra
!       For the gamma related statistics.
        common/anly_gamma1/ egam,  egam1,  egam2,  egam3, &
                           segam, segam1, segam2, segam3
!       For the parini related statistics.
        common/anly_parini1/ woun, swoun, snpctl0, snpctlm, snpar
!       For the simulation control.
        COMMON/SA1_PY8/ i_mode, i_tune, KF_woDecay(100), &
               KF_proj, KF_targ, win, energy_B, psno, b_min, b_max
        double precision bst(4),bzp,bzt
!       iii : number of current event
!       csen: e+p total x section in fm^2
!       neve : total number of events
!       bp : impact parameter
!       'sbe': store initial parton confiquration (with diquark)
!       'saf': store parton configuration after parton rescattering
!              (w/o diquark)
!       c17(i,1-3) : three position of i-th nucleon
!       tp(i) : time of i-th nucleon counted since collision of two nuclei
!       p17(i,1-4) : four momentum of i-th nucleon
!       ishp(i)=1 if i-th particle inside the simulated volume
!              =0 if i-th particle outside the simulated volume
!       tau(i): formation time of particle i.
!       ecsen: largest collision distance between lepton and p
!       note: incident lepton collides with nucleon in nucleus once
!        only, because of very low total x-section. that collision is the one
!        with lowest minimum approaching distance.
!       cspipiKK (fm^2): cross section of pion + pion to kaon + kaon
!       edipi: largest interaction distance between two pions.
!       epin: largest interaction distance between pion and nucleon.
!       ekn: largest interaction distance between kaon and nucleon.
!       ecsnn: largest interaction distance between two nucleons.
!       t0 : average proper formation time at rest.
!       dep : the accuracy in four momentum conservation
!       ddt : time accuracy used in parton initiation
!       time accuracy used in parton cascade is dddt
!       rou0 : normal nucleon density.
!       rao : enlarged factor for the radius of simulated volume.
!       nap and nzp (nat and nzt) are the mass and charge numbers of
!        projectile (target) nucleus
!       r0p=rnp : the standard radius of projectile
!       r0t=rnt : the standard radius of target
!       nctl: number of collision pairs in current collision list
!       nctl0: number of collision pairs in initial collision list
!       nctlm: maxmimum number of collision pairs
!       noinel(1): statistics of nn elas. colli.;
!       noinel(592): statistics of nn colli. calling PYTHIA
!       noinel(593): statistics of nn colli. not calling PYTHIA
!       suppm: the upper bound in sampling the radius of projectile nucleon
!       suptm: the upper bound in sampling the radius of target nucleon
!       suppc: the maximum radius in sample for projectile
!       suptc: the maximum radius in sample for target
!       coor: CM position of collision system
!       kfmax: the maximum # of particle KF code considered
!       kfaco(i): KF code of i-th particle
!       numb(i): order number of last particle among particles with same flavor
!        code of kfaco(i) in particle list


!-------------------------------------------------------------------------------
!-----------------------   Local Variable Initializing   -----------------------
!       Counters of sub-collisions pairs (zero-pT, nn, pp, np/pn and lp/ln).
        kpar = 0
        knn  = 0
        kpp  = 0
        knp  = 0
        kep  = 0
!       To avoide infinite loop in parini.
        ijk = 0
        tp   = 0D0
        N    = 0
        nbe  = 0
        naf  = 0
        nsa  = 0
        nbh  = 0
        idi  = 0
        idio = 0
        ndiq = 0
        npt  = 0
        ifcom = 0
        nstr0 = 0
        nstr1 = 0
        nctl = 0
!-----------------------   Local Variable Initializing   -----------------------
!-------------------------------------------------------------------------------


!       Initiates the pp (pA,Ap,AA,lp & lA) collision system.

!       Creats the initial particle (nucleon) list.

        call initialize_position

        call initialize_momentum( win, energy_B )

!       Calculates the velocity of the CM of collision system in LAB or
!        in nucleon-nucleon CM system.
        bst(1) = p17(1,1)*nap + p17(nap+1,1)*nat
        bst(2) = p17(1,2)*nap + p17(nap+1,2)*nat
        bst(3) = p17(1,3)*nap + p17(nap+1,3)*nat
        bst(4) = p17(1,4)*nap + p17(nap+1,4)*nat
        bst(1) = -bst(1)/bst(4)
        bst(2) = -bst(2)/bst(4)
        bst(3) = -bst(3)/bst(4)


!-------------------------------------------------------------------------------
!-----------------------   Particle Properties Giving   ------------------------
!      '1 -> |nzp|' are projectile protons or lepton, '|nzp|+1 -> nap'
!        are projectile neutrons; 'nap+1 -> nap+nzt' are targer protons,
!        the rest are target nuctrons in 'PYJETS' after nuclear initiation above
        napt = nap + nat
        n = napt
        do i=1,N,1
            K(i,1) = 1
            K(i,2) = 2112
            P(i,5) = PYMASS( 2212 )
!       For NN, NA(AN) and AB.
            if( (i <= abs(nzp) .and. ipden < 2) &
                 .OR. (i > nap .and. i <= nap+nzt) )then
                K(i,2) = 2212
                P(i,5) = PYMASS( 2212 ) 
!       For l+N & lbar + N
            else if( i <= nap .AND. (ipden >= 11 .AND. ipden <= 16 &
                    .AND. ABS(nzp) == 1) )then
                K(i,2) = SIGN( ipden, -nzp )
                P(i,5) = PYMASS( ipden )
            end if
            do j=1,3,1
                p(i,j) = p17(i,j)
                v(i,j) = c17(i,j)
            end do
            p(i,4) = p17(i,4)
            v(i,4) = tp(i)
        end do
!       v, vbh and vsa arraies are the position four vector
!       note: for v etc., we do not take care of their fifth component now.
!        for array k, we take care of only first three components now.
!-----------------------   Particle Properties Giving   ------------------------
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
!---------------------   CMS Boost & Lorentz Contraction   ---------------------
!       boost PYJETS into cms of initial nucleus-nucleus collision system
!        from lab or initial nucleon-nucleon cms system.
!       call pyrobo( 1, n, 0D0, 0D0, bst(1), bst(2), bst(3) )
!       Lorentz contraction
        bzp3 = 0D0
        bzp4 = 0D0
        bzt3 = 0D0
        bzt4 = 0D0
        do i=1,nap,1
            bzp3 = bzp3 + p(i,3)
            bzp4 = bzp4 + p(i,4)
        end do
        do i = nap+1, napt, 1
            bzt3 = bzt3 + p(i,3)
            bzt4 = bzt4 + p(i,4)
        end do
        bzp = bzp3 / bzp4
        bzt = bzt3 / bzt4
        gamp = 1.d0 / dsqrt( dmax1(1.d-20, (1.0d0 - bzp*bzp) ) )
!       no Lorentz contraction for incident lepton
        if(ipden >= 2) gamp = 1.
        gamt = 1.d0 / dsqrt( dmax1(1.d-20, (1.0d0 - bzt*bzt) ) )
!       try no lorentz contract for target
!       gamt = 1.
        do i=1,nap
            c17(i,3) = c17(i,3) / gamp
            v(i,3)   = v(i,3)   / gamp
        enddo
        do i = nap+1, napt
            c17(i,3) = c17(i,3) / gamt
            v(i,3)   = v(i,3)   / gamt
        enddo

!       Positions two particles at ( b/2, 0, -20 ) and ( -b/2, 0, 20 ).
!       20 fm is just a large enough distance in the z-direction.
        do i=1,nap,1
            V(i,3) = V(i,3) - 20D0
        end do
        do i = nap+1, napt, 1
            V(i,3) = V(i,3) + 20D0
        end do
!---------------------   CMS Boost & Lorentz Contraction   ---------------------
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
!----------------------   Particle Filtering & Ordering   ----------------------
!       filter out those kind of particles wanted to study and make
!        the order of proton, neutron, ... (cf. 'filt')
        call filt
!       since lepton was moved to last position after calling filt, one has to
!        remove it to the fist position
        if(ipden >= 2) call ltof(n)
!       'PYJETS' to 'sa2'
        nsa = n
        do j=1,5
            do i=1,n
                ksa(i,j) = k(i,j)
                psa(i,j) = p(i,j)
                vsa(i,j) = v(i,j)
            enddo
        enddo
        do i=1,n
            ishp(i) = 1
            tau(i)  = 0D0
        enddo
        numb = numbs
!----------------------   Particle Filtering & Ordering   ----------------------
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
!       note: particle list is composed of the arraies in common block
!        'sa2', the array 'ishp' in common block 'wz', the array 'tau' in
!        common block 'sa4', and the array 'numb' in common block 'sa5'
        time = time_ini
!       full_events_history of OSC1999A
        call oscar(0)
!       calculate the position for the center of mass of the
!        non-freeze-out system. The distance of a particle, when checking
!        is it freezing out or not, is measured with respect to this center
        call copl(time)
!       creat the initial collision list, note: be sure that the initial
!       collision list must not be empty
        call ctlcre
!       Resamples a 'b' and positions now.
        if( nctl <= 0 )then
            iii = iii - 1
            ijk = 1
            return
        end if
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
!------------------------   System Time Initializing   -------------------------
!       time origin is set at the time of first NN/lN collision (1D-5 fm/c)
!       find out colli. pair with least colli. time
        call find(icp,tcp,0)
        if(icp == 0) stop 'initial collision list is empty'
        time = tcp

!       perform classical Newton motion in Lab. system for all particles, i.e.
!        bring the two nuclei (or the nucleon/leton and a nucleus) into contact.
        call his(time,istop)
        do ij=1,nsa
            vsa(ij,4) = 1D-5
        enddo

!       move origin of time to collision time of first nucleon-nucleon collision
        do ij=1,nctl
            tc(ij) = tc(ij) - time + 1D-5
        enddo
        time = time_ini
        call copl(time)
!------------------------   System Time Initializing   -------------------------
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
!-------------------------   Collision Implementing   --------------------------
!       loop over implementing NN (hh), ll & lN collision, updating hadron list,
!        and updating collision time list untill it is empty.
!       It is equivalent to implementing a nucleus-nucleus collision.
        call scat(time,ijk,ipau)

        if(ijk == 1) return
        time_ini = time   ! 081010

!       'saf' to 'PYJETS'
        if(mstptj == 0) call tran_saf
        naf = 0
!       'sa2' to 'sbh'
!       'sbh' stores spectators, hadrons & leptons from the diffractive & the
!         special sub-processes, e.g. NRQCD onia, or hadron beam remnants, or
!         B-framework.
        nbh = 0
        if(nsa > 0)then
            nbh = nsa
            do i2=1,5
                do i1=1,nsa
                    kbh(i1,i2) = ksa(i1,i2)
                    pbh(i1,i2) = psa(i1,i2)
                    vbh(i1,i2) = vsa(i1,i2)
                enddo
            enddo
            nsa = 0
        endif

!       boost PYJETS back to lab or nucleon-nucleon cms system.
!       P(N,5) = SQRT( MAX( -P(N,1)**2 -P(N,2)**2 -P(N,3)**2 +P(N,4)**2, 0D0 ) )
!       P(N-1,5) = SQRT( MAX( - P(N-1,1)**2 - P(N-1,2)**2 - P(N-1,3)**2 &
!                + P(N-1,4)**2,0.0))
!       call PYBORO( 1, N, 0D0, 0D0, -bst(1), -bst(2), -bst(3) )

!-----------------------   A, B & D-framework Treating   -----------------------
        if( mstptj == 1 )then
!       PYTHIA like simulation for NA (AN) & AB, and low energy simulation.
!       'sbh' to 'PYJETS'
            n=nbh
            if(n >= 1)then
                do i2=1,5
                    do i1=1,n
                        k(i1,i2) = kbh(i1,i2)
                        p(i1,i2) = pbh(i1,i2)
                        v(i1,i2) = vbh(i1,i2)
                    enddo
                enddo
            endif
        endif
        ! Removes photons from "PYJETS" to "sgam".
        if( i_mode == 1 ) call remo_gam(22)
        ! Appends "sbh" to "PYJETS".
        if( i_mode == 4 )then
            do l=1,nbh,1
                l1 = N + l
                do m=1,5,1
                    K(l1,m) = kbh(l,m)
                    P(l1,m) = pbh(l,m)
                    V(l1,m) = vbh(l,m)
                enddo
            enddo
            N = N + nbh
        end if
!-----------------------   A, B & D-framework Treating   -----------------------

!-------------------------   Collision Implementing   --------------------------
!-------------------------------------------------------------------------------


!       egam1: energy of gamma after partonic initiation
        if( mstptj == 0 .AND. i_mode /= 4 )then
            call prt_sgam(ngam,egam1,1)
        end if
!       egam3: gamma energy after hadronization
        if( (mstptj == 1 .AND. i_mode /= 1) .OR. i_mode == 4 )then
            call prt_sgam(ngam,egam3,3)
        end if


!-------------------------------------------------------------------------------
!------------------------   Wounded Nucleons Counting   ------------------------
        if( i_mode /= 8 .AND. i_mode /= 9 )then
!       woun: # of wounded nucleons
!       unwoun: # of unwounded nucleons
            unwoun = 0D0
            do i=1,nbh,1
                KF = kbh(i,2)
                if( KF == 2212 .OR. KF == 2112 )then
                    pT2_N = pbh(i,1)**2 + pbh(i,2)**2
                    if( pT2_N <= 1D-12 ) unwoun = unwoun + 1D0
                end if
            end do
            woun = (nap + nat)*1D0 - unwoun
        end if
        if( mstptj == 1 .OR. i_mode == 4 ) nbh = 0
!------------------------   Wounded Nucleons Counting   ------------------------
!-------------------------------------------------------------------------------


        return
        end



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine sysini(win)
!!      Gives the initial values to quantities needed in calculation.
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        COMMON/PYCIDAT1/KFACOT(100),DISDET(100),ISINELT(600)
        common/sa5/kfmax,kfaco(100),numb(100),numbs(100),non5, &
         disbe(100,100)
        common/count/isinel(600)
        COMMON/PYCIDAT2/KFMAXT,NONCI2,PARAM(20),WEIGH(600)
        common/sa6/kfmaxi,nwhole
        common/papr/t0,cspipiKK,dep,ddt,edipi,epin,ecsnn,ekn,ecspsn,ecspsm, &
         rnt,rnp,rao,rou0,vneu,vneum,ecsspn,ecsspm,ecsen
        common/syspar/ipden,itden,suppm,suptm,suppc,suptc,r0p,r0t, &
         nap,nat,nzp,nzt,pio
        common/sa10/csnn,cspin,cskn,cspipi,cspsn,cspsm,rcsit,ifram, &
         iabsb,iabsm,i_sigma_AQM,ajpsi,csspn,csspm,csen
        common/sa25/i_inel_proc,i_time_shower,para1_1,para1_2
!       For cross sections of hh collisions etc.
        common/para_h1/ para(20)


        anat = nat
        anap = nap

!       rou0 = PARAM(11)

!       considering the nucleus as a sphere with radii rnt for target
!        and rnp for projectile.
!       rnt = ( 3D0*anat / ( 4D0 * pio * rou0 ) )**(0.33333)
!       rnp = ( 3D0*anap / ( 4D0 * pio * rou0 ) )**(0.33333)
        ! 1.05 to 1.12
        rp00 = 1.12D0
        rt00 = 1.12D0
        ! rp00 = 1.122D0 (nat=208)
        ! rt00 = 1.12D0  (nat=197)
!       if(nap > 16) rp00 = 1.16D0 * ( 1D0 - 1.16D0 * anap**(-0.666666) )
!       if(nat > 16) rt00 = 1.16D0 * ( 1D0 - 1.16D0 * anat**(-0.666666) )

!       Uses PDG RPP2024 charge radius of proton 0.8409 +- 0.0004 fm.
!           (PDG RPP2024 magnetic radius of neutron 0.864 +0.009 -0.008 fm)
!       if(itden == 0) rnt = rt00 * anat**(0.33333)
        if(itden == 0) rnt = 0.841D0
        ! +0.54
        if(itden == 1) rnt = rt00 * anat**(0.33333)
        ! 2.60 2.095 1.54, deuteron
        if(nat == 2 .and. nzt == 1) rnt = 4.0D0
        ! lepton
        if(ipden >= 2) rnt = 0.5D0
!       if(ipden == 0) rnp = rp00*anap**(0.33333)
        if(ipden == 0) rnp = 0.841D0
        ! +0.54
        if(ipden == 1) rnp = rp00*anap**(0.33333)
        ! lepton
        if(ipden >= 2) rnp = 0.5D0
        ! 2.60 2.095 1.54
        if(nap == 2 .and. nzp == 1) rnp = 4.0D0
        rou0 = 3D0 / 4D0 / pio * anat / rnt**3
        r0p = rnp
        r0t = rnt

!       Sets initial values to some quantities
!       in the program the x-sections are given in a unit of fm^2.
        PARAM(1) = para1_1
        PARAM(5) = para(5)
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
!       1 mb = 0.1 fm^2
!       cspin : total cross section of pion + nucleon
!       cskn  : total cross section of Kaon + nucleon
!       cspipi: total cross section of pion + pion
!       cspsn : total cross section of J/psi (psi') + nucleon
!       cspsm : total cross section of J/psi (psi') + meson
!       iabsb = 0 : without J/psi (psi') + baryon (nucleon)
!             = 1 : with J/Psi (psi')    + baryon (nucleon)
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
        if( i_sigma_AQM == 3 .OR. i_sigma_AQM == 4 )then
            ! sigma_NN_tot = para1_1
            sigma_NN_tot = para1_2
            cspin  = sigma_AQM( 211, 2212, sigma_NN_tot ) * 0.1D0
            cskn   = sigma_AQM( 321, 2212, sigma_NN_tot ) * 0.1D0
            cspipi = sigma_AQM( 211,  211, sigma_NN_tot ) * 0.1D0
            cspsn  = sigma_AQM( 443, 2212, sigma_NN_tot ) * 0.1D0
            cspsm  = sigma_AQM( 443,  221, sigma_NN_tot ) * 0.1D0
            csspn  = cspsn
            csspm  = cspsm
        end if

        if( ipden >= 2 )then
            if(ifram == 0)then
                ept=sqrt(win*win+0.938*0.938)
                rots=sqrt((ept+0.938)*(ept+0.938)-win*win)
            endif
            if(ifram == 1) rots=win
            ! temporary using e^-p total x-section
            call crosep(rots,csen)
            ! e^-p total x-section
!           if(nzp < 0) call crosep(rots,csen)
            ! e^+p total x-section
!           if(nzp >= 0) call crosepp(rots,csen)
            csen = csen * 0.1D0
        endif

!       largest collision distance between two colliding particles.
        edipi  = dsqrt(cspipi / pio)
        epin   = dsqrt(cspin  / pio)
        ekn    = dsqrt(cskn   / pio)
        ecsnn  = dsqrt(csnn   / pio)
        ecspsn = dsqrt(cspsn  / pio)
        ecspsm = dsqrt(cspsm  / pio)
        ecsspn = dsqrt(csspn  / pio)
        ecsspm = dsqrt(csspm  / pio)
        ecsen  = dsqrt(csen   / pio)

        anp = nap**0.3333D0
        ant = nat**0.3333D0
        do ia=1,2
            if(ia == 1) napt = nap
            if(ia == 2) napt = nat
            ! In parini.f90.
            alpt = get_nuclear_skin_depth( napt )
            if(ia == 1) alp = alpt
            if(ia == 2) alt = alpt
        end do

!       Sets max R=10 fm (max b=20 fm) directly, i.e. bmax = 20 fm.
        ! suppc = rp00*anp + 2D0*alp
        ! suptc = rt00*ant + 2D0*alt
        suppc = 10D0
        suptc = 10D0
        suppm = 1D0 / ( 1D0 + EXP(0D0 - r0p/alp) )
        suptm = 1D0 / ( 1D0 + EXP(0D0 - r0t/alt) )

        rcsit = PARAM(6)
        t0    = PARAM(7)
!       t0 = 0D0
!       proper formation time of particle from 'PYTHIA'
        dep = PARAM(9)
        ddt = PARAM(8)
        rao = PARAM(10)
        kfmax  = KFMAXT
        kfaco  = KFACOT
        isinel = ISINELT
        disbe  = 0D0
        do j=1,kfmax
!       something might be missing here ?
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
        subroutine initialize_position
!!      Initializes positions of particles, i.e. distributes nucleons inside
!!       nuclei (and other particles).
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        PARAMETER (KSZJ=80000)
        common/sa1/kjp21,non1,bp,iii,neve,nout,nosc
        common/sa24/adj1(40),nnstop,non24,zstop
        common/sa30/vneump,vneumt,mstptj
        common/sa33/smadel,ecce,secce,parecc,iparres
        common/syspar/ipden,itden,suppm,suptm,suppc,suptc,r0p,r0t, &
         nap,nat,nzp,nzt,pio
         common/wz/c17(500,3),ishp(kszj),tp(500),coor(3),p17(500,4)


!       initialization of x, y,xy, x^2, y^2 and sump (statistics of the
!         number of nucleons in overlap region, the initial geometry)
        sumx  = 0D0
        sumy  = 0D0
        sumxy = 0D0
        sumx2 = 0D0
        sumy2 = 0D0
        sump  = 0D0
        c17   = 0D0


!-------------------------------------------------------------------------------
!--------------------------   Position Initializing   --------------------------
!       in position phase space
        iadj130 = INT(adj1(30))

!*****************************   A+B Collisions   ******************************
!       A+B (nucleus-nucleus)
        if( ipden == 1 .and. itden == 1 )then
!       distribute projectile nucleons by Woods-Saxon
            napt = nap
            ! In parini.f90.
            alpt = get_nuclear_skin_depth( napt )
            alp = alpt
            r0  = r0p
            ! upper bound in sampling the radius of projectile nucleon
            am = suppm
            ! maximum radius for projectile
            ac = suppc
            ! ratio of projectile participant nucleons to total
            ratps = vneump / nap
            do i1=1,nap
                if( iadj130 == 1 )then
                    rann = pyr(1)
                    if(rann < ratps)then
!       sample position of projectile nucleon in overlap region of nuclei
                        call arrove(i1,1,alp,r0,am,ac)
                    else
!       sample position of projectile nucleon according to Woods-Saxon
!        distribution
                        call woodsax_samp(i1,1,alp,r0,am,ac,1)
!       last argument here is 'iway', iway=1: particle i1 must be outside the
!        overlap region of colliding nuclei, iway=0: no more requirement
                    endif
                elseif( iadj130 == 0 )then
                    call woodsax_samp(i1,1,alp,r0,am,ac,0)
                endif
            enddo
!       distribute target nucleons by Woods-Saxon
            napt = nat
            alpt = get_nuclear_skin_depth( napt )
            alp = alpt
            r0  = r0t
            ! upper bound in sampling the radius of target
            am  = suptm
            ! maximum radius for target
            ac  = suptc
            ! ratio of target participant nucleons to total
            ratps = vneumt / nat
            do i1=1,nat
                i2 = i1 + nap
                if( iadj130 == 1 )then
                    rann = pyr(1)
                    if(rann < ratps)then
!       sample position of target nucleon in overlap region of colliding nuclei
                        call arrove(i2,0,alp,r0,am,ac)
                    else
!       sample position of target nucleon according to Woods-Saxon
!        distribution
                        call woodsax_samp(i2,0,alp,r0,am,ac,1)
                    endif
                elseif( iadj130 == 0 )then
                    call woodsax_samp(i2,0,alp,r0,am,ac,0)
                endif
            enddo

!-------------   Impact-Parameter & Initial Geometry Calculating   -------------
            do i=1,nap
                ! c17(i,1)=c17(i,1)+bp
                ! move x-component of origin to 0.5*bp
                c17(i,1) = c17(i,1) + 0.5*bp
!       Calculates eccentricity correctly for both adj1(30)=0 and 1.
                x = c17(i,1)
                y = c17(i,2)
                z = c17(i,3)
!       Relative distance between the projectile nucleon i and the
!        target center (-bp/2., 0, 0)
                rel_dist = SQRT( (x+bp/2.)**2 + y**2 + z**2 )
!       The projectile nucleon i is inside the target, i.e. inside the
!        overlap region.
                if( rel_dist  <=  r0t )then
                    sumx  = sumx  + x
                    sumy  = sumy  + y
                    sumxy = sumxy + x*y
                    sumx2 = sumx2 + x**2
                    sumy2 = sumy2 + y**2
                    sump  = sump  + 1D0
                end if
            enddo
            do i = nap+1, nap+nat
                c17(i,1) = c17(i,1) - 0.5*bp
!       Calculates eccentricity correctly for both adj1(30)=0 and 1.
                x = c17(i,1)
                y = c17(i,2)
                z = c17(i,3)
!       Relative distance between the target nucleon i and the
!        projectile center (+bp/2., 0, 0)
                rel_dist = SQRT( (x-bp/2.)**2 + y**2 + z**2 )
!       The target nucleon i is inside the projectile, i.e. inside the
!        overlap region.
                if( rel_dist  <=  r0p )then
                    sumx  = sumx  + x
                    sumy  = sumy  + y
                    sumxy = sumxy + x*y
                    sumx2 = sumx2 + x**2
                    sumy2 = sumy2 + y**2
                    sump  = sump  + 1D0
                end if
            enddo
!-------------   Impact-Parameter & Initial Geometry Calculating   -------------

!*****************************   A+B Collisions   ******************************

!***************************   NA & lA Collisions   ****************************
!       N+A or lepton+A
        elseif( (ipden == 0.or.ipden > 1) .and. itden == 1 )then
!       distribute projectile proton
            do i=1,3
                c17(1,i) = 0D0
                if(i == 1) c17(1,i) = c17(1,i) + 0.5*bp
            enddo
!       distribute target nucleons by Woods-Saxon
            napt = nat
            alpt = get_nuclear_skin_depth( napt )
            alp = alpt
            r0  = r0t
            ! upper bound in sampling the radius of target
            am  = suptm
            ! maximum radius for target
            ac  = suptc
            do i1=1,napt
                i2 = i1 + nap
                call woodsax_samp(i2,0,alp,r0,am,ac,0)
            enddo
            do i = nap+1, nap+nat
                c17(i,1) = c17(i,1) - 0.5*bp
            enddo
!***************************   NA & lA Collisions   ****************************

!******************************   AN Collisions   ******************************
!       A+N
        elseif( ipden == 1 .and. itden == 0 )then
!       distribute projectile nucleons by Woods-Saxon
            napt = nap
            alpt = get_nuclear_skin_depth( napt )
            alp = alpt
            r0  = r0p
            ! upper bound in sampling the radius of projectile nucleon
            am  = suppm
            ! maximum radius for projectile
            ac  = suppc
            do i1=1,napt
                call woodsax_samp(i1,1,alp,r0,am,ac,0)
            enddo
            do i=1,napt
                c17(i,1) = c17(i,1) + 0.5*bp
            enddo
            do i=1,3
                c17(nap+1,i) = 0D0
                if(i == 1) c17(nap+1,i) = -0.5*bp
            enddo
!******************************   AN Collisions   ******************************

!***************************   NN & lN Collisions   ****************************
!       N+N or lepton+N
!       Now NN and e+e- collisions have been processed in main.f90 directly.
        elseif( (ipden == 0  .and. itden == 0) .or. &
                (ipden >= 11 .and. itden == 0) )then
            do i=1,3
                c17(1,i) = 0D0
                c17(2,i) = 0D0
            enddo
        endif
!***************************   NN & lN Collisions   ****************************

        r0pt = r0p + r0t

!**********************   Initial Geometry Calculating   ***********************
        if(sump > 0D0)then   !!!
            asumx  = sumx / sump
            sigmx2 = sumx2/sump - asumx*asumx
            asumy  = sumy / sump
            sigmy2 = sumy2/sump - asumy*asumy
            asumxy = sumxy / sump
            sigmxy = asumxy - asumx*asumy
            sigmsu = sigmy2 + sigmx2
            sigmde = sigmy2 - sigmx2
            argu   = sigmde*sigmde + 4D0*sigmxy*sigmxy
!       participant eccentricity of participant nucleons
            if(argu > 0. .and. sigmsu > 0.) &
                ecce = sqrt(argu) / sigmsu
!       calculate transverse overlap area
            argu1 = sigmx2*sigmy2 - sigmxy*sigmxy
            if(argu1 > 0.) secce = 3.1416*sqrt(argu1)

!---------------------   Transverse-Momentum Deformation   ---------------------
!       assuming ecce = geometric eccentricity of ellipsoid
!           sqrt( 1-b^2/a^2 )
!        with half major axis
!           b = pt * ( 1 + smadel )
!        and half minor axis
!           a = pt * ( 1 - smadel ),
!        the resulted
!           smadel = -ecce*ecce/4,
!        neglecting the samll term of
!           ecce*ecce*( -2*smadel + smadel*smadel ).
            ecc2 = ecce*ecce
            ! approximated deformation parameter
            smadel_a = parecc*ecc2/4.
            delta1 = (2. - ecc2 + 2.*(1. - ecc2)**0.5 ) / ecc2
            delta2 = (2. - ecc2 - 2.*(1. - ecc2)**0.5 ) / ecc2
            if(delta1 <= 1.)then
                smadel = parecc*delta1  ! exact deformation parameter
            elseif(delta2 <= 1.)then
                smadel = parecc*delta2  ! exact deformation parameter
            endif
!       here a sign change is introduced because of asymmetry of initial
!        spatial space is oppsed to the final momentum space
!---------------------   Transverse-Momentum Deformation   ---------------------

        endif   !!!
!**********************   Initial Geometry Calculating   ***********************

!       the beam direction is identified as the z axis
!       the origin in position space is set on (0, 0, 0)
!       and the origin of time is set at the moment of
!       first nn colission assumed to be 1.e-5
!--------------------------   Position Initializing   --------------------------
!-------------------------------------------------------------------------------


        return
        end



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        function get_nuclear_skin_depth( nA )
!!      Gets the nuclear skin depth. nA: atomic mass number.
        integer :: nA
        real(kind=8) :: depth, get_nuclear_skin_depth


        if( nA < 27 )then
            depth = 0.47D0
        elseif( nA >= 27 .AND. nA <= 108 )then
            depth = 0.488D0
        else
            depth = 0.54D0
        endif

        if( nA == 27 )then
            depth = 0.478D0
        else if( nA == 28 )then
            depth = 0.48D0
        else if( nA == 32 )then
            depth = 0.49D0
        else if( nA == 56 )then
            depth = 0.49D0
        else if( nA == 64 )then
            depth = 0.49D0
        else if( nA == 108 )then
            depth = 0.495D0
        else if( nA == 184 )then
            depth = 0.53D0
        else if( nA == 197 )then
            depth = 0.54D0
        else if( nA == 207 )then
            depth = 0.545D0
        else if( nA == 238 )then
            depth = 0.55D0
        end if

        get_nuclear_skin_depth = depth


        return
        end



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine arrove(ii,jj,alp,r0,am,ac)
!!      Arranges randomly particle ii in overlap region of colliding nuclei
!        jj=0 and 1 for target and projectile, respectively
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        PARAMETER (KSZJ=80000)
        common/sa1/kjp21,non1,bp,iii,neve,nout,nosc
        common/wz/c17(500,3),ishp(kszj),tp(500),coor(3),p17(500,4)
        common/syspar/ipden,itden,suppm,suptm,suppc,suptc,r0p,r0t, &
         nap,nat,nzp,nzt,pio


        b = bp
        iiii = 0
54      iiii = iiii + 1
        if( iiii == 10000 )then
            write(9,*) "Warning, difficult to arrange produced nucleons in" &
                    // "subroutine arrove, infinite loop may occur"
            return
        end if
!       ii in target (-b/2.)
        if(jj == 0)then
!       sample a point according to woodsax distribution
            call woodsax_samp(ii,jj,alp,r0,am,ac,0)
            x = c17(ii,1)
            y = c17(ii,2)
            z = c17(ii,3)
!       relative to projectile center, they are b-x, y, and z, respectively
!       adjudge does (x-b,y,z) is in the sphere of projectile
            r1 = sqrt( (b-x)*(b-x) + y*y + z*z )
            if(r1 > r0p) goto 54
        endif
!       ii in projectile (+b/2.)
        if(jj == 1)then
!       sample a point according to woodsax distribution
            call woodsax_samp(ii,jj,alp,r0,am,ac,0)
            x = c17(ii,1)
            y = c17(ii,2)
            z = c17(ii,3)
!       relative to target center, they are x+b, y, and z, respectively
!       adjudge does (x+b,y,z) is in the sphere of target
            r1 = sqrt( (x+b)*(x+b) + y*y + z*z )
            if(r1 > r0t) goto 54
        endif


        return
        end



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine woodsax_samp(ii,jj,alp,r0,am,ac,iway)
!!      Samples position of nucleon ii in nucleus according to
!!       Woods-Saxon distribution
!       jj=0 and 1 for target and projectile, respectively
!       alp: diffusion length
!       r0: radius of nucleus
!       am: upper bound in sampling the radius
!       ac: maximum radius
!       iway=1: ii must be outside overlap region of colliding nuclei
!       iway=0: no more requirement
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        PARAMETER (KSZJ=80000)
        common/sa1/kjp21,non1,bp,iii,neve,nout,nosc
        common/wz/c17(500,3),ishp(kszj),tp(500),coor(3),p17(500,4)
        common/syspar/ipden,itden,suppm,suptm,suppc,suptc,r0p,r0t, &
         nap,nat,nzp,nzt,pio


        b = bp
!       selecting sample method
!       Woods-Sax distribution
        iiii = 0
100     iiii = iiii + 1
        if( iiii == 100000 )then
            write(9,*) "Warning, difficult to arrange produced nucleons in" &
                    // "subroutine woodsax, infinite loop may happen"
            return
        end if

        a1 = pyr(1)
        xf = ac*(a1)**( 1./3. )
        b1 = pyr(1)
        deno2 = 1. + exp( (xf-r0) / alp )
!       if(deno2 == 0.) deno2 = 1.e-10
        yf = 1. / deno2
!       Gaussian distribution
!       yf = exp( -xf*xf / 2. /r0 )
        if(b1 > yf/am) goto 100
        call samp(xf,ii)
!       subroutine 'samp' is sampling the direction according to isotropic
!        distribution
        x = c17(ii,1)
        y = c17(ii,2)
        z = c17(ii,3)

!       ii must be outside overlap region of colliding nuclei
        if(iway == 0) return

        ! ii in target (-b/2.)
        if(jj == 0)then
!       relative to projectile center, above x, y, and z are b-x, y, and z,
!        respectively
!       (b-x,y,z) is inside or not inside the sphere of projectile
            r1 = sqrt( (b-x)*(b-x) + y*y + z*z )
            if(r1 < r0p) goto 100
        endif
        ! ii in projectile (+b/2.)
        if(jj == 1)then
!       relative to target center, they are x+b, y, and z, respectively
!       (x+b,y,z) is inside or not inside the sphere of target
            r1 = sqrt( (x+b)*(x+b) + y*y + z*z )
            if(r1 < r0t) goto 100
        endif


        return
        end



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine samp(xf,i)
!!      Arranges i-th particle on the surface of sphere with radius xf
!!       sampling on the surface of a sphere with radius xf
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        PARAMETER (KSZJ=80000)
        common/wz/c17(500,3),ishp(kszj),tp(500),coor(3),p17(500,4)
        common/papr/t0,cspipiKK,dep,ddt,edipi,epin,ecsnn,ekn,ecspsn,ecspsm, &
         rnt,rnp,rao,rou0,vneu,vneum,ecsspn,ecsspm,ecsen   ! 060813
        common/syspar/ipden,itden,suppm,suptm,suppc,suptc,r0p,r0t, &
         nap,nat,nzp,nzt,pio


        cita = 2D0 * PYR(1) - 1D0
        fi   = 2D0 * pio * PYR(1)
        sita = SQRT( 1D0 - cita**2 )
        c17(i,1) = xf * sita * COS(fi)
        c17(i,2) = xf * sita * SIN(fi)
        c17(i,3) = xf * cita


        return
        end



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine sampi(xf,i)
!!      Sampling in a sphere with radius xf.
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        PARAMETER (KSZJ=80000)
        common/wz/c17(500,3),ishp(kszj),tp(500),coor(3),p17(500,4)


100     continue
        x = 2D0 * PYR(1) - 1D0
        y = 2D0 * PYR(1) - 1D0
        z = 2D0 * PYR(1) - 1D0
        rr = x*x + y*y + z*z
        if( rr > 1D0 ) goto 100
        c17(i,1) = xf * x
        c17(i,2) = xf * y
        c17(i,3) = xf * z


        return
        end



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine initialize_momentum( win, energy_B )
!!      Initializes momenta of particles (nucleons inside nuclei and
!!       other particles).
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        PARAMETER (KSZJ=80000)
        common/sa1/kjp21,non1,bp,iii,neve,nout,nosc
        common/sa10/csnn,cspin,cskn,cspipi,cspsn,cspsm,rcsit,ifram, &
         iabsb,iabsm,i_sigma_AQM,ajpsi,csspn,csspm,csen
        common/sa21/pincl(5),pscal(5),pinch(5),vnu,fq2,w2l,yyl,zl,xb, &
         pph,vnlep
        common/wz/c17(500,3),ishp(kszj),tp(500),coor(3),p17(500,4)
        common/syspar/ipden,itden,suppm,suptm,suppc,suptc,r0p,r0t, &
         nap,nat,nzp,nzt,pio
!       For the charge and 4-momentum check.
        common/cp_check/ p_init(4), p_final(4)


!-------------------------------------------------------------------------------
!--------------------------   Momentum Initializing   --------------------------
!       In momentum phase space.
        ep1 = 0D0
        ep2 = 0D0
        et1 = 0D0
        et2 = 0D0
        pp1 = 0D0
        pp2 = 0D0
        pt1 = 0D0
        pt2 = 0D0
        p17 = 0D0
        ! Sets the default frame as CMS.
        if( ifram > 2 ) ifram = 1

!       Fixed target in lab.
        if( ifram == 0 )then
            pp1 = win     ! momentum of projetile particle (if proton)
            pt1 = 1D-20   ! momentum of target    particle (if proton)
            pp2 = win     ! momentum of projetile particle (if neutron)
            pt2 = 1D-20   ! momentum of target    particle (if neutron)
            pm2 = PYMASS(2212)**2         ! square mass of proton
            ep1 = SQRT( pp1*pp1 + pm2 )   ! energy of projetile (proton)
            et1 = SQRT( pt1*pt1 + pm2 )   ! energy of target    (proton)
            pm2 = PYMASS(2112)**2         ! square mass of neutron
            ep2 = SQRT( pp2*pp2 + pm2 )   ! energy of projetile (neutron)
            et2 = SQRT( pt2*pt2 + pm2 )   ! energy of target    (neutron)
!       Sets four momentum and mass for incident lepton.
            if( ipden >= 11 .AND. ipden <= 16 )then
                pincl(1) = 0D0
                pincl(2) = 0D0
                pincl(3) = win
                pincl(5) = PYMASS(ipden)
                pincl4   = pincl(3)*pincl(3) + pincl(5)*pincl(5)
                pincl4   = dmax1( pincl4, 1D-20 )
                pincl(4) = dsqrt(pincl4)
                pinch(1) = 0D0
                pinch(2) = 0D0
                pinch(3) = 0D0
                pinch(5) = PYMASS(2212)
                pinch(4) = pinch(5)
            end if
        end if

!       In cms.
        if( ifram == 1 )then
            ep1 = 0.5D0*win   ! energy of projetile particle (if proton)
            et1 = ep1         ! energy of target    particle (if proton)
            ep2 = 0.5D0*win   ! energy of projetile particle (if neutron)
            et2 = ep2         ! energy of target    particle (if neutron)
            pm2 = PYMASS(2212)**2          ! square mass of proton
            pp1 =  SQRT( ep1*ep1 - pm2 )   ! momentum of projetile (proton)
            pt1 = -SQRT( et1*et1 - pm2 )   ! momentum of target    (proton)
            pm2 = PYMASS(2112)**2          ! square mass of nucleon
            pp2 =  SQRT( ep2*ep2 - pm2 )   ! momentum of projetile (neutron)
            pt2 = -SQRT( et2*et2 - pm2 )   ! momentum of target    (neutron)
!       Sets four momentum and mass for incident lepton.
            if( ipden >= 11 .AND. ipden <= 16 )then
                pincl(1) = 0D0
                pincl(2) = 0D0
                pincl(4) = 0.5D0*win
                pincl(5) = PYMASS(ipden)
                pincl3   = pincl(4)*pincl(4) - pincl(5)*pincl(5)
                pincl3   = dmax1( pincl3, 1D-20)
                pincl(3) = SQRT(pincl3)
                pinch(1) = 0D0
                pinch(2) = 0D0
                pinch(4) = 0.5D0*win
                pinch(5) = PYMASS(2212)
                pinch3   = pinch(4)*pinch(4) - pinch(5)*pinch(5)
                pinch3   = dmax1( pinch3, 1D-20 )
                pinch(3) = SQRT(pinch3)
            end if
        end if

!       Back-to-back in lab.
        if( ifram == 2 )then
            eA = win
            eB = energy_B
            ep1 = eA
            et1 = eB
            ep2 = eA
            et2 = eB
            if( ep1 < PYMASS(2212) )then
                if( iii == 1 )then
                    write(22,*)
                    write(22,*) "Warning! win(energy_A) < m_pronton " // &
                                "in the back-to-back collision frame. " // &
                                "Re-evaluated the energy " // &
                                "according to m_proton."
                    write(22,*)
                end if
                ep1 = PYMASS(2212)
            end if
            if( et1 < PYMASS(2212) )then
                if( iii == 1 )then
                    write(22,*)
                    write(22,*) "Warning! energy_B < m_pronton " // &
                                "in the back-to-back collision frame. " // &
                                "Re-evaluated the energy " // &
                                "according to m_proton."
                    write(22,*)
                end if
                et1 = PYMASS(2212)
            end if
            if( ep2 < PYMASS(2112) )then
                if( iii == 1 )then
                    write(22,*)
                    write(22,*) "Warning! win(energy_A) < m_neutron " // &
                                "in the back-to-back collision frame. " // &
                                "Re-evaluated the energy " // &
                                "according to m_neutron."
                    write(22,*)
                end if
                ep2 = PYMASS(2112)
            end if
            if( et2 < PYMASS(2112) )then
                if( iii == 1 )then
                    write(22,*)
                    write(22,*) "Warning! energy_B < m_neutron " // &
                                "in the back-to-back collision frame. " // &
                                "Re-evaluated the energy " // &
                                "according to m_neutron."
                    write(22,*)
                end if
                et2 = PYMASS(2112)
            end if
            pm2 = PYMASS(2212)**2
            pp1 =  SQRT( ep1*ep1 - pm2 )
            pt1 = -SQRT( et1*et1 - pm2 )
            pm2 = PYMASS(2112)**2
            pp2 =  SQRT( ep2*ep2 - pm2 )
            pt2 = -SQRT( et2*et2 - pm2 )
            if(       ( ABS(pp1) < 1D-10 .AND. ABS(pp2) < 1D-10 ) &
                .AND. ( ABS(pp1) < 1D-10 .AND. ABS(pt2) < 1D-10 ) &
                .AND. ( ABS(pt1) < 1D-10 .AND. ABS(pp2) < 1D-10 ) &
                .AND. ( ABS(pt1) < 1D-10 .AND. ABS(pt2) < 1D-10 ) )then
                write(22,*)
                write(22,*) "Abort! Too low win(energy_A) and energy_B!"
                write(22,*)
                stop
            end if
            if( ipden >= 11 .AND. ipden <= 16 )then
                pincl(1) = 0D0
                pincl(2) = 0D0
                pincl(4) = eA
                pincl(5) = PYMASS(ipden)
                pincl3   = pincl(4)*pincl(4) - pincl(5)*pincl(5)
                pincl3   = MAX( pincl3, 1D-20 )
                pincl(3) = SQRT( pincl3 )
                pinch(1) = 0D0
                pinch(2) = 0D0
                pinch(4) = eB
                pinch(5) = PYMASS(2212)
                pinch3   = pinch(4)*pinch(4) - pinch(5)*pinch(5)
                pinch3   = MAX( pinch3, 1D-20 )
                pinch(3) = SQRT( pinch3 )
            end if
        end if

        inzp = ABS(nzp)
        inzt = ABS(nzt)
        do i=1,nap,1
            p17(i,1) = 0D0
            p17(i,2) = 0D0
            if( i <= inzp )then
                p17(i,3) = pp1
                p17(i,4) = ep1
            else
                p17(i,3) = pp2
                p17(i,4) = ep2
            end if
            ! The initial total 4-momentum.
            p_init(3) = p_init(3) + p17(i,3)
            p_init(4) = p_init(4) + p17(i,4)
        enddo
        napt = nap + nat
        do i = nap+1, napt, 1
            p17(i,1) = 0D0
            p17(i,2) = 0D0
            if( i <= nap+inzt )then
                p17(i,3) = pt1
                p17(i,4) = et1
            else
                p17(i,3) = pt2
                p17(i,4) = et2
            end if
            ! The initial total 4-momentum.
            p_init(3) = p_init(3) + p17(i,3)
            p_init(4) = p_init(4) + p17(i,4)
        end do
!--------------------------   Momentum Initializing   --------------------------
!-------------------------------------------------------------------------------


        return
        end



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine filt
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
!       In case of lepton+A, one images lepton as an initial projectile proton.
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        PARAMETER (KSZJ=80000)
        common/sa5/kfmax,kfaco(100),numb(100),numbs(100),non5, &
         disbe(100,100)
        common/sa6/kfmaxi,nwhole
        COMMON/PYJETS/N,NPAD,K(KSZJ,5),P(KSZJ,5),V(KSZJ,5)
        COMMON/PYCIDAT2/KFMAXT,NONCI2,PARAM(20),WEIGH(600)
        SAVE/PYCIDAT2/


        iiii = 0
        jjjj = 0
        do i=1,kfmax,1
            KF = kfaco(i)
            do j = iiii+1, N, 1
                call ord( jjjj, j, KF )
            end do
            iiii = jjjj
            numbs(i) = jjjj
        end do


        return
        end



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine ord(ipi,j,kf)
!!      Orders particles according to flavor code.
!       j: the particle needed to order
!       ipi: j-th particle should order after ipi
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        PARAMETER (KSZJ=80000)
        COMMON/PYJETS/N,NPAD,K(KSZJ,5),P(KSZJ,5),V(KSZJ,5)
        dimension pp(5),vv(5),kk(5)


        ik=k(j,2)
        if(ik == kf)then
            ipi=ipi+1
            do jj=1,5
                kk(jj)=k(ipi,jj)
                pp(jj)=p(ipi,jj)
                vv(jj)=v(ipi,jj)
            enddo
            do jj=1,5
                k(ipi,jj)=k(j,jj)
                p(ipi,jj)=p(j,jj)
                v(ipi,jj)=v(j,jj)
            enddo
            do jj=1,5
                k(j,jj)=kk(jj)
                p(j,jj)=pp(jj)
                v(j,jj)=vv(jj)
            enddo
        endif


        return
        end



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine prt_pyj(nn,cc)
!!      Prints particle list and sum of momentum and energy.
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        INTEGER PYCHGE
        PARAMETER (KSZJ=80000)
        COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
        common/PYJETS/nsa,nonsa,ksa(kszj,5),psa(kszj,5),vsa(kszj,5)
        dimension peo(4)


        do i=1,nn
            write(mstu(11),*)i,ksa(i,2),(psa(i,j),j=1,4)
        enddo
        call psum(psa,1,nsa,peo)
        ich1=0
        do i1=1,nn
            kf=ksa(i1,2)
            ich1=ich1+pychge(kf)
        enddo
        cc=ich1/3D0
        write(22,*) 'pyj nn=',nn
        write(mstu(11),*) 'c & p sum=',cc,peo


        return
        end



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine prt_delt(nn,cc)
!!      Prints particle list and sum of momentum and energy.
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        INTEGER PYCHGE
        PARAMETER (KSZJ=80000)
        COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
        common/delt/nsa,nonsa,ksa(kszj,5),psa(kszj,5),vsa(kszj,5)
        dimension peo(4)


        do i=1,nn
            write(mstu(11),*)i,ksa(i,2),(psa(i,j),j=1,4)
        enddo
        call psum(psa,1,nsa,peo)
        ich1=0
        do i1=1,nn
            kf=ksa(i1,2)
            ich1=ich1+pychge(kf)
        enddo
        cc=ich1/3D0
        write(22,*) 'delt nn=',nn
        write(mstu(11),*) 'c & p sum=',cc,peo


        return
        end



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine prt_sbh(nn,cc)
!!      Prints particle list and sum of momentum and energy.
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        INTEGER PYCHGE
        PARAMETER (KSZJ=80000)
        COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
        common/sbh/nbh,nonbh,kbh(kszj,5),pbh(kszj,5),vbh(kszj,5)
        dimension peo(4)


        do i=1,nn
            write(mstu(11),*)i,kbh(i,2),(pbh(i,j),j=1,4)
        enddo
        call psum(pbh,1,nbh,peo)
        ich1=0
        do i1=1,nn
            kf=kbh(i1,2)
            ich1=ich1+pychge(kf)
        enddo
        cc=ich1/3D0
        write(22,*) 'sbh nn=',nn
        write(mstu(11),*) 'c & p sum=',cc,peo


        return
        end



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine prt_sa2(nn,cc)
!!      Prints particle list and sum of momentum and energy.
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        INTEGER PYCHGE
        PARAMETER (KSZJ=80000)
        COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
        common/sa2/nbh,nonbh,kbh(kszj,5),pbh(kszj,5),vbh(kszj,5)
        dimension peo(4)


        do i=1,nn
            write(22,*)i,kbh(i,2),(pbh(i,j),j=1,4)
        enddo
        call psum(pbh,1,nbh,peo)
        ich1=0
        do i1=1,nn
            kf=kbh(i1,2)
            ich1=ich1+pychge(kf)
        enddo
        cc=ich1/3D0
        write(22,*) 'sa2 nn=',nn
        write(22,*) 'c & p sum=',cc,peo


        return
        end



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine prt_sbe(nn,cc)   ! 220110
!!      Prints particle list and sum of momentum and energy.
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        INTEGER PYCHGE
        PARAMETER (KSZJ=80000)
        COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
        common/sbe/nbh,nonbh,kbh(kszj,5),pbh(kszj,5),vbh(kszj,5)
        dimension peo(4)


        call psum(pbh,1,nbh,peo)
        ich1=0
        do i1=1,nn
            kf=kbh(i1,2)
            ich1=ich1+pychge(kf)
        enddo
        cc=ich1/3D0
        ! write(22,*)'sbe nn=',nn
        write(22,*)'c & p sum=',cc,peo   !
        do i=1,nn
            write(22,*)i,kbh(i,1),kbh(i,2),(pbh(i,j),j=1,4)
            enddo
            write(22,*) "--------------------------------------------"// &
                        "--------------------------------------------"// &
                        "--------------------------------------------"// &
                        "------------"


        return
        end



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine prt_saf(nn,cc)   ! 220110
!!      Prints particle list and sum of momentum and energy.
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        INTEGER PYCHGE
        PARAMETER (KSZJ=80000)
        COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
        common/saf/nbh,nonbh,kbh(kszj,5),pbh(kszj,5),vbh(kszj,5)
        dimension peo(4)


        do i=1,nn
            write(22,*)i,kbh(i,2),(pbh(i,j),j=1,4)
        enddo
        call psum(pbh,1,nbh,peo)
        ich1=0
        do i1=1,nn
            kf=kbh(i1,2)
            ich1=ich1+pychge(kf)
        enddo
        cc=ich1/3D0
        write(22,*) 'saf nn=',nn
        write(22,*) 'c & p sum=',cc,peo


        return
        end



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine psum(pei,il,ih,peo)
!!      Calculates sum of momentum and energy.
!       pei: two dimension array of input momentum and energy
!       il and ih: lower and upper limits of sum
!       peo : one dimension array of output momentum and energy
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        PARAMETER (KSZJ=80000)
        dimension pei(kszj,5),peo(4)


        peo = 0D0
        do i=il,ih
            do j=1,4
                peo(j)=peo(j)+pei(i,j)
            end do
        end do


        return
        end



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine scat(time,ijk,ipau)
!!      Loops over implementing NN (hh) collision, updating hadron list, and
!!       updating collision time list untill the collision time list is empty.
!!      It is equivalent to implementing a nucleus-nucleus collision.
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        LOGICAL IS_PYTHIA8,IS_EXIST,IS_DIQUARK
        PARAMETER (KSZJ=80000)
        PARAMETER (NSIZE=750000)
        COMMON/PYDAT3/MDCY(500,3),MDME(8000,2),BRAT(8000),KFDP(8000,5)
        COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
        COMMON/PYSUBS/MSEL,MSUB(500),KFIN(2,-40:40),NON,CKIN(200)
        COMMON/PYPARS/MSTP(200),PARP(200),MSTI(200),PARI(200)
        COMMON/PYJETS/N,NPAD,K(KSZJ,5),P(KSZJ,5),V(KSZJ,5)
        common/papr/t0,cspipiKK,dep,ddt,edipi,epin,ecsnn,ekn,ecspsn,ecspsm, &
         rnt,rnp,rao,rou0,vneu,vneum,ecsspn,ecsspm,ecsen
        common/syspar/ipden,itden,suppm,suptm,suppc,suptc,r0p,r0t, &
         nap,nat,nzp,nzt,pio
        common/wz/c17(500,3),ishp(kszj),tp(500),coor(3),p17(500,4)
        common/sa1/kjp21,non1,bp,iii,neve,nout,nosc
        common/sa2/nsa,non2,ksa(kszj,5),psa(kszj,5),vsa(kszj,5)
        common/sa4/tau(kszj),tlco(kszj,4)
        common/sa5/kfmax,kfaco(100),numb(100),numbs(100),non5, &
         disbe(100,100)
        common/sa6/kfmaxi,nwhole
        common/sa7/ispmax,isdmax,iflmax,ispkf(100),n_bin_hist,asd(10), &
         afl(100,10,2),i_y_or_eta
        common/sa10/csnn,cspin,cskn,cspipi,cspsn,cspsm,rcsit,ifram, &
         iabsb,iabsm,i_sigma_AQM,ajpsi,csspn,csspm,csen
        common/sa12/ppsa(5),nchan,nsjp,sjp,taup,taujp
        common/sa13/kjp20,non13,vjp20,vjp21,vjp22,vjp23
        common/sa15/nps,npsi,pps(5000,5),ppsi(5000,5)
        common/sa16/x_ratio,decpro,dni(10),dpi(10),edi(10),bmin,bmax, &
         bar(10),abar(10),barf(10),abarf(10),emin(10),eminf(10), &
         eplu(10),epluf(10),nmax
        common/sa18/i_deex,n_deex_step,i_pT_coal,i_pT_endpoint,a_FF,aPS_c,aPS_b
        common/sa21/pincl(5),pscal(5),pinch(5),vnu,fq2,w2l,yyl,zl,xb, &
          pph,vnlep
        common/sa23/kpar,knn,kpp,knp,kep,NON_p,skpar,sknn,skpp,sknp,skep
        common/sa24/adj1(40),nnstop,non24,zstop
        common/sa25/i_inel_proc,i_time_shower,para1_1,para1_2
        common/sa26/ndiq(kszj),npt(kszj),ifcom(kszj),idi,idio
        common/sa27/itime,kjp22,gtime,astr,akapa(6),parj1,parj2,parj3, &
         parj21,parj4,adiv,gpmax,nnc
        common/sa28/nstr,nstra(kszj),nstrv(kszj),nstr0, &
         nstr1,nstr1a(kszj),nstr1v(kszj)
        common/sa29/i_color_reconnection,NCR,parp2
        common/sa30/vneump,vneumt,mstptj
        common/sa33/smadel,ecce,secce,parecc,iparres
        common/sa34/itorw,iikk,cp0,cr0,kkii
        common/sbe/nbe,nonbe,kbe(kszj,5),pbe(kszj,5),vbe(kszj,5)
        common/saf/naf,nonaf,kaf(kszj,5),paf(kszj,5),vaf(kszj,5)
        common/sbh/nbh,nonbh,kbh(kszj,5),pbh(kszj,5),vbh(kszj,5)
        common/sgam/ngam,nongam,kgam(kszj,5),pgam(kszj,5),vgam(kszj,5)
        common/ctllist/nctl,noinel(600),nctl0,nctlm
        common/ctllist_t/ lc(nsize,5), tc(nsize), tw(nsize)
        common/coal1/bmrat,i_mm
!       For the simulation control.
        COMMON/SA1_PY8/ i_mode, i_tune, KF_woDecay(100), &
               KF_proj, KF_targ, win, energy_B, psno, b_min, b_max
        dimension pi(4),pj(4),pint(4)
        dimension b(3),pl(100,5)
        logical flag_hh, flag_lh, INEL_coll
!       arraies in 'PYJETS' are given after calling 'PYTHIA'
!       arraies in 'sa2' are used in the collision processes
!       arraies in 'sbh' are used to store hadron after calling 'PYTHIA'
!       numbs(i) is given in 'filt', updated with transport processes, and
!        numbs(i)->numb(i) in the initiation of nucleus-nucleus collisin only
!        numb(i) is updated with transport processes
!        i_mm: only the hadrons below numb(i_mm) (cf. 'filt' in parini.f90)
!        join in updating the hh collision time list (cf. parini.f90).
!        Originally i_mm=2, i.e. only p (proton) and n (neutron) join in
!        updating the collision time list in parton initialization;
!        if i_mm=6, p, n, pbar, nbar, pi+ and pi- will join in updating
!        the collision time list in parton initialization.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!       0. is the hard distance between two pions
!       0.5 is the hard distance between two nucleons
!       0. is the hard distance between pion and nucleon
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!       lc(i,1) and lc(i,2) are the line # of colliding particles 1 and 2 of
!        i-th collision pair in particle list, respectively
!       lc(i,3) and lc(i,4) are the flavor codes of scattered particles 3 and 4
!        of i-th collision, respectively
!       lc(i,5) identifies the different inelastic processes,
!       lc(i,5)=592 refers to the process calling 'PYTHIA'
!       tc(i) is the collision time of i-th colli.
!       tw(i) is the cross section ratio of (i-th inelas.)/tot
!       array 'sbe' stores cumulatively parton (q,qq,g and their anti-particle)
!        configuration before breaking the diquarks
!       array 'saf' stores cumulatively parton (q,g and their anti-particle)
!        configuration after breaking the diquarks
!       idi: counts cumunatively the number of diquark (anti-diquark)
!       idio: value of idi after last nn collision
!       ndiq(j): = 0 if j is quark (antiquark)
!                = idi if j is diquark (anti-diquark)
!       note: j is line number in 'sbe' ('saf')
!       ifcom(idi): line number of first component of idi-th diquark
!       npt(idi): line number of second component of idi-th diquark
!        (anti-diquark) in 'sbe' ('saf')
!       nstr1: statitics of number of strings in a nucleus-nucleus collision
!        when fragmentation string-by-string
!       nstr1a(i): line number of first component of i-th string
!       nstr1v(i): line number of last component of i-th string
!       nstr0: number of strings after call break
!       iparres: = 0/1, with diquarks breaking-up; 2/3, without.


!-------------------------------------------------------------------------------
!-----------------------   Local Variable Initializing   -----------------------
        ijk = 0
        m1 = numb(1)
        m2 = numb(2)
        m3 = numb(3)
        m4 = numb(4)
        m6 = numb(6)
        m7 = numb(7)
        nctl0 = nctl
        nctlm = nctl0
!       The hadronization model.
        i_had_model = INT( adj1(12) )
!-----------------------   Local Variable Initializing   -----------------------
!-------------------------------------------------------------------------------



!*******************************************************************************
!---------------------------   Sub-Collisions Loop   ---------------------------
!       loop over hadron-hadron (NN), ll, lN sub-collisions.
        iiii = 1
10      if(iiii > 1) call copl(time)
!       find out the binary colli. with minimum collsion time
        call find(icp,tcp,1)
        if(icp == 0) goto 100
!       icp=0 means the collision list is empty
        naf00  = naf
        ngam00 = ngam
        if( i_mode == 4 )then
            idi   = 0
            idio  = 0
            ndiq  = 0
            npt   = 0
            ifcom = 0
            nstr1 = 0
            nbe   = 0
            naf   = 0
            nbh   = 0
        end if
        l  = lc(icp,1)
        l1 = lc(icp,2)
!       l and l1: line numbers of current hh colliding pair in 'sa2'
!       Prompts the initialization information of an NN collision.
        if(iii == 1 .and. iiii == 1)then
            MSTP(122) = 1
        else
            MSTP(122) = 0
        end if
        time0 = time
        kfa   = ksa(l,2)
        kfb   = ksa(l1,2)
        ikfa  = iabs(kfa)
        ikfb  = iabs(kfb)
        kfaab = iabs(kfa)
        kfbab = iabs(kfb)
!       record this collision time
        time  = tcp

!       tlco(l,4)  = tcp
!       tlco(l1,4) = tcp


!-------------------------------------------------------------------------------
!------------------------   Kinematics Pre-Estimating   ------------------------
        ilo = 0
        pi(4) = psa(l,4)
        pj(4) = psa(l1,4)
!       pi (pj): four momentum of current colliding hadron of l (l1)
        if(pi(4) < 1.e-20) pi(4)=1.e-20   ! 041204
        if(pj(4) < 1.e-20) pj(4)=1.e-20   ! 041204
        do i=1,3
            pi(i) = psa(l,i)
            pj(i) = psa(l1,i)
            b(i)  = ( pi(i) + pj(i) ) / ( pi(4) + pj(4) )
        end do
        pti = dsqrt( pi(1)**2 + pi(2)**2 )
        ptj = dsqrt( pj(1)**2 + pj(2)**2 )

!       Recalculates energies to avoid potential errors from the hadrons
!        with bad energies due to last "PYTHIA" execution.
        pi4_R = SQRT( pi(1)**2 + pi(2)**2 + pi(3)**2 + PYMASS(kfa)**2 )
        pj4_R = SQRT( pj(1)**2 + pj(2)**2 + pj(3)**2 + PYMASS(kfb)**2 )
        ss    = SQRT( ( pi4_R + pj4_R )**2 &
                    - ( pi(1) + pj(1) )**2 &
                    - ( pi(2) + pj(2) )**2 &
                    - ( pi(3) + pj(3) )**2 )

!       boost to CMS frame of colliding pair from CMS or Lab. frame of heavy
!        ion collision system
        call lorntz(ilo,b,pi,pj)
!       ss = pi(4)+pj(4)
!       if(ss < 1.e-18 ) ss = 1D-18

!       calculate the angular 'theta' of the momenta pi and pj
        ctai  = PYANGL( pi(3), SQRT( pi(1)**2 + pi(2)**2 ) )
        ctaj  = PYANGL( pj(3), SQRT( pj(1)**2 + pj(2)**2 ) )
        cctai = COS(ctai)
        cctaj = COS(ctaj)
        if( cctai > 0D0 )then
!       calculate the 'orentation' of the vector pi
            call codi(pi,cfi1,sfi1,ccta1,scta1)
        else
            call codi(pj,cfi1,sfi1,ccta1,scta1)
        end if
!------------------------   Kinematics Pre-Estimating   ------------------------
!-------------------------------------------------------------------------------


!       Performs classical Newton motion.
        call his( time, istop )
        if(istop == 1) goto 100
!       istop = 1 means all particles have get out of considered volume

!       Updated numb(i).
        ! Upto p.
        m1 = numb(1)
        ! Upto p and n.
        m2 = numb(2)
        ! Upto p, n and pbar.
        m3 = numb(3)
        ! Upto p, ..., nbar.
        m4 = numb(4)
        ! Upto p, ..., pi+.
        m5 = numb(5)
        ! Upto p, ..., pi-.
        m6 = numb(6)
        ! Upto p, ..., pi0.
        m7 = numb(7)


!-------------------------------------------------------------------------------
!---------------------------   B, C & D-framework   ----------------------------
!       High energy frameworks.
        if( i_mode /= 1 )then
            ! Enough energy to call PYTHIA (inelastic).
            INEL_coll = ss >= parp2
            ! Hadrons to be considered for a hh collision.
            flag_hh = ipden < 2 .AND. ( l <= numb(i_mm) .AND. l1 <= numb(i_mm) )
            ! Leptons/hadrons to be considered  for a lh collision.
            flag_lh = ipden > 2 &
                    .AND. ( ( kfaab >= 11 .AND. kfaab <= 16 .AND. l1 <= m2 ) &
                    .OR.    ( kfbab >= 11 .AND. kfbab <= 16 .AND. l  <= m2 ) )

!**************************   Inelastic Collisions   ***************************
            if( INEL_coll .AND. ( flag_hh .OR. flag_lh ) )then

!-----------------------------   Event Executing   -----------------------------
!       Executes collision event.
                call xevent( l, l1, iiii, i_error )
!       Hard to generate. Gives up this collision pair.
                if( i_error == 1 )then
                    do j1= icp+1, nctl, 1
                        j = j1 - 1
                        tc(j) = tc(j1)
                        tw(j) = tw(j1)
                        do m=1,5,1
                            lc(j,m) = lc(j1,m)
                        end do
                    end do
                    nctl = nctl - 1
                    naf  = naf00
                    ngam = ngam00
                    goto 10
                end if
!       Gives four position to the particles after calling PAEVNT
!        for particles in 'PYJETS'
                call ptcre(l,l1,time)
!       592-th scattering process is referred to calling 'PYTHIA'
                noinel(592) = noinel(592) + 1
!-----------------------------   Event Executing   -----------------------------

!-----------------------------   Lepton Counting   -----------------------------
!       statistics of number of leptons studied, identify scattered lepton,
!        and fill up pscal(5)
                if( ipden >= 11.and.ipden <= 16 )then
!       identify the studied leptons
                    kfl = ipden
                    if(nzp > 0) kfl = -ipden
                    nlep = 0
                    do j=1,N,1
                        ikl = k(j,2)
                        if(ikl == kfl)then
                            nlep = nlep + 1
                            pl(nlep,1) = p(j,1)
                            pl(nlep,2) = p(j,2)
                            pl(nlep,3) = p(j,3)
                            pl(nlep,4) = p(j,4)
                            pl(nlep,5) = p(j,5)
                        end if
                    end do
!       find the scattered lepton (with largest energy among studied leptons)
                    if( nlep > 1 )then
                        vnlep = vnlep + nlep
                        elep = 1D0
                        jj = 0
                        pscal = 0D0
                        do j1=1,nlep,1
                            plj14 = pl(j1,4)
                            if(plj14 >= elep)then
                                elep = plj14
                                jj = j1
                            end if
                        end do
                        if(jj > 0)then
                            do j2=1,5
                                pscal(j2) = pl(jj,j2)
                            end do
                        end if
                    else if(nlep == 1)then
                        vnlep = vnlep + nlep
                        do j2=1,5
                            pscal(j2) = pl(nlep,j2)
                        end do
                    end if
!       calculate kinematic variables relevant to incident and scattered
!        lepton only, in cms
                    pdotk = pinch(4)*pincl(4) - pinch(1)*pincl(1) &
                          - pinch(2)*pincl(2) - pinch(3)*pincl(3)   ! P.k
                    q11   = pincl(1) - pscal(1)
                    q22   = pincl(2) - pscal(2)
                    q33   = pincl(3) - pscal(3)
                    q44   = pincl(4) - pscal(4)
                    q112  = q11 * q11
                    q222  = q22 * q22
                    q332  = q33 * q33
                    q442  = q44 * q44
                    pdotq = pinch(4)*q44 - pinch(1)*q11 - pinch(2)*q22 &
                          - pinch(3)*q33   ! P.q
                    vnu   = pdotq / pinch(5)   ! \nu
                    fq2   = -( q442 - q112 - q222 - q332 )   ! Q^2=-q^2
                    w2l   = ( pinch(4) + q44)**2 - ( pinch(1) + q11 )**2 &
                          - ( pinch(2) + q22)**2 - ( pinch(3) + q33 )**2   ! W^2
                    pdotk = dmax1( pdotk, 1D-20 )
                    yyl   = pdotq / pdotk   ! y
                    pdotq = dmax1( pdotq, 1D-20 )
                    xb    = fq2 / 2D0 / pdotq   ! x_b
                end if
!-----------------------------   Lepton Counting   -----------------------------

!------------------   B-framework Gamma 22 & Hadron Removing  ------------------
                if( mstptj == 1 )then
                    ! Removes gamma from "PYJETS" to "sgam".
                    call remo_gam(22)
                    ! "PYJETS" to "sbh".
                    do im=1,5,1
                        do li=1,N,1
                            kbh(li,im) = K(li,im)
                            pbh(li,im) = P(li,im)
                            vbh(li,im) = V(li,im)
                        end do
                    end do
                    nbh = N
                    N = 0
                    ! For B-framework.
                    goto 997
                end if
!------------------   B-framework Gamma 22 & Hadron Removing  ------------------

!----------------------   C-framework Gamma 44 Removing  -----------------------
                ! Do not remove any entries for CR if C-framework via sfm.
                if( i_had_model /= 0 .OR. ( i_had_model == 0 .AND. kjp22 < 5 ) &
                    .OR. i_mode == 4 )then
                    ! Removes gamma (labeled as "44") from "PYJETS" to "sgam".
                    n44=0
                    do j=1,N,1
                        kf = K(j,2)
                        if(kf == 22)then
                            K(j,2) = 44
                            n44 = n44 + 1
                        end if
                    end do
                    if(n44 > 0) call remo_gam(44)
!----------------------   C-framework Gamma 44 Removing  -----------------------

!-----------------------------   debug   ------------------------------
                    if( i_had_model ==1 .AND. i_deex > 99 )then
                        ! In coales.f90
                        call do_debug( i_deex )
!-----------------------------   debug   ------------------------------

!-----------------------   C-framework Hadron Removing  ------------------------
!       Removes hadrons/leptons from "PYJETS" to "sbh".
                    else
                        call remo
                    end if
                end if
!-----------------------   C-framework Hadron Removing  ------------------------

!----------------------------   Diquark Locating   -----------------------------
!       Prepares for the parton rescattering.
                if( i_had_model == 0 .OR. i_had_model == 3 )then
                    i3 = nbe
                    do i1=1,N,1
                        i3 = i1 + nbe
                        KS = K(i1,1)
                        KF = K(i1,2)
                        ! Identifies diquarks.
                        ndiq( i1 + naf ) = 0
                        if( IS_DIQUARK(KF) .AND. IS_EXIST(KS,i_mode) )then
                            idi = idi + 1
                            ndiq( i1 + naf ) = idi
                        end if
                        ! 'PYJETS' to 'sbe'
                        do i2=1,5
                            kbe(i3,i2) = K(i1,i2)
                            pbe(i3,i2) = P(i1,i2)
                            vbe(i3,i2) = V(i1,i2)
                        end do
                    end do
                    nbeo = nbe
                    nbe = i3
                end if
!----------------------------   Diquark Locating   -----------------------------

!---------------------------   Diquark Breaking-up   ---------------------------
!       Breaks up diquark and give four momentum and four position
!        to the broken quarks (working in 'PYJETS').
                if( iparres == 0 .OR. iparres == 1 )  call break
!---------------------------   Diquark Breaking-up   ---------------------------

!-----------------------------   String Locating   -----------------------------
!       Finds number of strings and line number of first and last components
!        of each string
!#TODO(Lei20240218): not work well for PYTHIA 8, need improvement.
                if( .NOT.IS_PYTHIA8(i_mode) .AND. i_had_model /= 1 )then
                    jb = 0
10000               continue
                    do i1 = jb+1, N, 1
                        ! i1 is 'A'
                        if( K(i1,1) == 2 )then
                            do i2 = i1+1, N, 1
                                ! i2 is 'V'
                                if( K(i2,1) == 1 )then
                                    nstr1 = nstr1 + 1
                                    ! line numbers of nstr1-th string
                                    nstr1a(nstr1) = i1 + naf
                                    nstr1v(nstr1) = i2 + naf
                                    jb = i2
                                    goto 10000
                                end if
                            end do
                        end if
                    end do
                end if
                nstr0 = nstr1
                ! nstr1: number of strings after call break
!-----------------------------   String Locating   -----------------------------

!-----------------------------   Pauli Blocking   ------------------------------
                ipau = 0
                ! Without Pauli blocking
                goto 777
!       Pauli effect (working in 'PYJETS' (current), in 'saf' (past))
                tpaul = 1D0
!       tpaul: product of the unoccupation probabilities
                ! Current
                do i1=1,N,1
                    kfp = K(i1,2)
                    kfp = ABS(kfp)
                    ppaul = 1D0
                    if( kfp == 1 .OR. kfp == 2 .OR. kfp == 3 )then
                        call pauli( i1, ppaul )
!       ppaul: the unoccupation probability of particle i1
                        ! over occupation, should be blocked
                        if( ppaul  <  0D0 )then
                            tpaul = 0D0
                            goto 666
                        end if
                    end if
                    tpaul = tpaul * ppaul
                end do
!       Blocked
666             if( PYR(1) >= tpaul )then
                    ipau = 1
!       Removes "current part of 'sbe'" from 'sbe' and truncates 'ndiq' and
!        'npt' correspondingly.
                    do i1 = nbeo+1, nbe, 1
                        do i2=1,5
                            kbe(i1,i2) = 0
                            pbe(i1,i2) = 0D0
                            vbe(i1,i2) = 0D0
                        end do
                        ndiq( i1 + naf ) = 0
                    end do
                    do i1=idio,idi,1
                        npt(i1) = 0
                        ifcom(i1) = 0   ! 220110
                    end do
                    nbe = nbeo
                    idi = idio
!       Removes the current NN collision pair from the collision list.
                    do j1 = icp+1, nctl, 1
                        j = j1 - 1
                        tc(j) = tc(j1)
                        tw(j) = tw(j1)
                        do m=1,5
                            lc(j,m) = lc(j1,m)
                        end do
                    end do
                    nctl = nctl - 1
                    goto 10
                end if
!       Unblocked
777             continue
                idio = idi
!-----------------------------   Pauli Blocking   ------------------------------

!--------------------------   CME Charge Separation   --------------------------
!       add CME charge separation for u d s, c
                icme = INT( adj1(23) )
                if( (nap == nat) .AND. (nzp == nzt) .AND. (icme == 1) ) &
                    call chargecme
!--------------------------   CME Charge Separation   --------------------------

!-------------------------   B & C-framework Treating  -------------------------
!       Appends 'PYJETS' to 'saf'. etc.
                if( i_mode /= 4 )then
                    do i1=1,N,1
                        naf = naf + 1
                        do i2=1,5,1
                            kaf(naf,i2) = K(i1,i2)
                            paf(naf,i2) = P(i1,i2)
                            vaf(naf,i2) = V(i1,i2)
                        end do
                    end do
                end if
!-------------------------   B & C-framework Treating  -------------------------

997             continue

!-------------------------------   D-framework   -------------------------------
                if( i_mode == 4 )then
                    call D_framework( iiii, naf00, l, l1, i_error )
!       Throws away the current NN collision pair from collision list
!        due to errors occurred in "D_framework".
                    if( i_error == 1 )then
                        do j1= icp+1, nctl, 1
                            j = j1 - 1
                            tc(j) = tc(j1)
                            tw(j) = tw(j1)
                            do m=1,5,1
                                lc(j,m) = lc(j1,m)
                            end do
                        end do
                        nctl = nctl - 1
                        naf  = naf00
                        ngam = ngam00
                        N = 0
                        nbe = 0
                        goto 10
                    end if
                end if
!-------------------------------   D-framework   -------------------------------

!       Counts cross sections.
                call PASTAT(-1,iiii)

!---------------------   Particle & Time Lists Updating   ----------------------
!       Updates hadron list 'sa2' after calling PYTHIA ('sbh' to 'sa2'),
!        removes collision pair composed of l and/or l1, removes l (l1)
!        from 'sa2'.
                call updpip(l,l1,time)
!       Updates collision list after calling 'PYTHIA'
                call updtlp(time,i_mode)   ! 250423
!       Resets "nbh=0" to avoid potential errors.
                nbh = 0
!---------------------   Particle & Time Lists Updating   ----------------------

!**************************   Inelastic Collisions   ***************************

!***************************   Elastic Collisions   ****************************
            else
!       If ss is not enough to call PYTHIA or current hadron-hadron
!        collision pair is not in the plan of calling PYTHIA,
!        the current hh collision is impremented as ela. scattering.
!       noinel(593): counter of # of hadron-hadron coll. not calling PYTHIA.
                noinel(593) = noinel(593) + 1
                call coelas(l,l1,ss,pi,pj)
!       Updates the particle list for elastic scattering, pi and pj have been
!        boosted back to Lab frame or cms of nucleus-nucleus collision.
                call updple(l,l1,b,pi,pj)
!       Statistics of the number of ela. hadron-hadron collisions for
!        high energy channel
                noinel(1) = noinel(1) + 1
!       Updates the collision list after ela. scattering.
                call updatl(l,l1,time,i_mode)
!       Recovers for D-framework.
                if( i_mode == 4 )then
                    naf  = naf00
                    ngam = ngam00
                end if
            end if
!***************************   Elastic Collisions   ****************************

        end if
!       High energy frameworks ended.
!---------------------------   B, C & D-framework   ----------------------------
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
!-------------------------------   A-framework   -------------------------------
!       Low energy framework.
        if( i_mode == 1 )then
            call simulate_hh_of_AB_in_framework_A( l, l1, pi, pj, &
                                time, b, ss, inorex, jorn )

!**********************   Inelastic Scattering Treating   **********************
            if( inorex /= 1 )then
!----------------------------   Rotation & Boost   -----------------------------
!       Originally, the orientation of scattered particle is relative to the
!        scattering particle, so before boost back to Lab. or cms of
!        nucleus-nucleus collision system one has first rotating to being
!        relative to the cms of incident channel
                do j=1,N,1
                    do j1=1,4
                        pint(j1) = P(j,j1)
                    end do
!       perform the rotate for inelas. scattered particles
                    call rosa(cfi1,sfi1,ccta1,scta1,cfis,sfis,cctas,sctas,pint)
                    do j1=1,4,1
                        P(j,j1) = pint(j1)
                    end do
                end do
!       boost back to Lab. or nucleus-nucleus collider
                ilo = 1
                do j=1,N,2
                    ! Even N.
                    do j1=1,4
                        pi(j1) = P(j,j1)
                        pj(j1) = P(j+1,j1)
                    end do
                    ! Odd N.
                    if( j == N ) pj = pi
                    call lorntz( ilo, b, pi, pj )
                    do j1=1,4,1
                        P(j,j1)   = pi(j1)
                        P(j+1,j1) = pj(j1)
                    end do
                end do
!----------------------------   Rotation & Boost   -----------------------------

!       Gives four position to the particles after after inelas. scattering
!        in 'PYJETS'.
                call ptcre(l,l1,time)
!       'PYJETS' to 'sbh'
                nbh=0
                nbh = N
                do mm1=1,5,1
                    do li=1,nbh,1
                        kbh(li,mm1) = K(li,mm1)
                        pbh(li,mm1) = P(li,mm1)
                        vbh(li,mm1) = V(li,mm1)
                    enddo
                enddo
!       update hadron list 'sa2' after inela. scattering  ('sbh' to 'sa2')
                call updpip_nn(l,l1)
!       update collision time list after inela. scattering
                call updatl_nn(l,l1,time,jorn,nsa0,i_mode)
            end if
!**********************   Inelastic Scattering Treating   **********************

        end if
!       Low energy framework ended.
!-------------------------------   A-framework   -------------------------------
!-------------------------------------------------------------------------------


!       Next sub-collision.
        if( nctl > 0 )then
            iiii = iiii + 1
            if( iiii > 100*(nctl0) )then
                write(22,*) "Warning, infinite loop may have happened in" &
                         // "subroutine scat iii=", iii
                iii = iii - 1
                ijk = 1
                return
            end if
            if( nctl > nctlm ) nctlm = nctl
            goto 10
        end if
100     continue
        call copl(time)
!---------------------------   Sub-Collisions Loop   ---------------------------
!*******************************************************************************


        return
        end



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine xevent( i_a, i_b, i_NN_pair, i_error )
!!      Executes a sub-collisions (lN, NN) event via PYTHIA.
!       Uses "3MOM" frame in PYTHIA 6 covering "CMS" and "FIXT".
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        INTEGER PYCHGE
        PARAMETER (KSZJ=80000,KSZJ_PY8=300000)
        LOGICAL IS_PYTHIA8,IS_EXIST
        COMMON/PYINT1/MINT(400),VINT(400)
        COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
        COMMON/PYPARS/MSTP(200),PARP(200),MSTI(200),PARI(200)
        COMMON/PYJETS/N,NPAD,K(KSZJ,5),P(KSZJ,5),V(KSZJ,5)
        common/syspar/ipden,itden,suppm,suptm,suppc,suptc,r0p,r0t, &
         nap,nat,nzp,nzt,pio
        common/sa1/kjp21,non1,bp,iii,neve,nout,nosc
        common/sa2/nsa,non2,ksa(kszj,5),psa(kszj,5),vsa(kszj,5)
        common/sa6_p/ithroq_p,ithrob_p,ich_p,non6_p,throe_p(4)
        common/sa18/i_deex,n_deex_step,i_pT_coal,i_pT_endpoint,a_FF,aPS_c,aPS_b
        common/sa23/kpar,knn,kpp,knp,kep,NON_p,skpar,sknn,skpp,sknp,skep
        common/sa24/adj1(40),nnstop,non24,zstop
        common/sa30/vneump,vneumt,mstptj
        common/sa33/smadel,ecce,secce,parecc,iparres
        common/sa34/itorw,iikk,cp0,cr0,kkii
        COMMON/SA1_PY8/ i_mode, i_tune, KF_woDecay(100), &
               KF_proj, KF_targ, win, energy_B, psno, b_min, b_max
        COMMON/PYJETS_PY8/ N_PY8, NPAD_PY8, &
               K_PY8(KSZJ_PY8,8), P_PY8(KSZJ_PY8,7), V_PY8(KSZJ_PY8,5)
        COMMON/PYJETS_AA/ N_AA(KSZJ_PY8), N_TOT_AA, N_COLL_NN_AA, &
               K_AA(KSZJ_PY8,8), P_AA(KSZJ_PY8,7), V_AA(KSZJ_PY8,5)
        dimension ps0(6), ps1(6), dp0(6)
        data iii_0 / 0 /
        data i_NN_pair_0 / 0 /
        data dp0 / 6*0D0 /


!       Hadronization model.
        i_had_model = INT( adj1(12) )
!       Backup.
        MSTJ21  = MSTJ(21)
        MSTP111 = MSTP(111)
!       Charge and 4-momentum before the collision.
        ps0 = 0D0
!       Charge and 4-momentum after the collision.
        ps1 = 0D0
!       Regeneration counter.
        i_xevent = 0
        if( iii > iii_0 )then
            i_NN_pair_0 = 0
            dp0 = 0D0
        end if
        iii_0 = iii
        i_error = 0


!-------------------------------------------------------------------------------
!----------------------   Colllision System Identifying   ----------------------
!       Gets names of particles a and b.
        kf_a = ksa(i_a,2)
        kf_b = ksa(i_b,2)
!       Mainly for collisions in PYTHIA 6 mode, especially when eCM < 10 GeV.
!       Sums of incident px, py, pz, E, inv. m, and charge.
        do i=1,3,1
            ps0(i) = psa(i_a,i) + psa(i_b,i)
        end do
        ps0(4) = SQRT( psa(i_a,1)**2 + psa(i_a,2)**2     &
               +       psa(i_a,3)**2 + PYMASS(kf_a)**2 ) &
               + SQRT( psa(i_b,1)**2 + psa(i_b,2)**2     &
               +       psa(i_b,3)**2 + PYMASS(kf_b)**2 )
        ps0(5) = SQRT( ps0(4)**2 - ps0(1)**2 - ps0(2)**2- ps0(3)**2 )
        ps0(6) = ( PYCHGE( kf_a ) + PYCHGE( kf_b ) ) / 3D0
!       Transverse momenta.
        pT_a = SQRT( psa(i_a,1)**2 + psa(i_a,2)**2 )
        pT_b = SQRT( psa(i_b,1)**2 + psa(i_b,2)**2 )
!----------------------   Colllision System Identifying   ----------------------
!-------------------------------------------------------------------------------


!       To accelerate the calculation of PYTHIA 8, the forward instantiation
!        has been done by calling SUBROUTINE PAINST only once in main.f90.
!#TODO(Lei20241008): PYTHIA 8 mode cannot handle free 3-mom collisions here.

100     N = 2
        do i=1,5,1
            P(1,i) = psa(i_a,i)
            P(2,i) = psa(i_b,i)
        end do


!-------------------------------------------------------------------------------
!----------------------------   Generation Failed   ----------------------------
!       Mainly for PYTHIA 6. No more than 100 attempts.
        i_xevent = i_xevent + 1
        if( i_xevent > 100 )then
            i_error = 1
            N = 0
            ! Recovers.
            MSTJ(21)  = MSTJ21
            MSTP(111) = MSTP111
            return
        end if
!----------------------------   Generation Failed   ----------------------------
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
!--------------------------   PYTHIA 6 Initializing   --------------------------
        iikk = 0
        kkii = 0
        MSTP(111) = mstptj
!       For the debug.
        if( i_had_model == 1 .AND. i_deex > 99 )then
            MSTJ(21)  = 0
            MSTP(111) = 1
        end if
!       Initilizes the colllision. (In p_40.f)
        if( .NOT.IS_PYTHIA8(i_mode) )then
            ! Uses "3MOM" directly.
            i_frame_3MOM = 3
            ! Just a dummy value.
            ss_dummy = ps0(4)
            MSTP(127) = 1
            call PAINIT_KF( i_frame_3MOM, kf_a, kf_b, ss_dummy )
            ! Too low energy for the PYTHIA 6 running.
            if( MSTI(53) == 1 )then
                i_error = 1
                N = 0
                ! Recovers.
                MSTJ(21)  = MSTJ21
                MSTP(111) = MSTP111
                return
            end if
        end if
!       The initialization of PYTHIA 8 has been done by calling
!        PAINIT_KF in main.f90.
!--------------------------   PYTHIA 6 Initializing   --------------------------
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
!------------------------   Core Processes Generating   ------------------------
!       Executes the collision in PYTHIA 6 or 8.
        if( IS_PYTHIA8(i_mode) )then
            MINT(5)   = iii
            MINT(355) = i_NN_pair
            MINT(11)  = kf_a
            MINT(12)  = kf_b
        end if
        ! In Pythia8_fort_interface.f90.
        call PAEVNT
!------------------------   Core Processes Generating   ------------------------
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
!----------------------------   PYTHIA 6 Checking   ----------------------------
!       Sums of px, py, pz, E, inv. m, and charge after the excution.
        ps1 = 0D0
        do i=1,6,1
            ps1(i) = PAPYP(0,i)
        end do
!       Mainly for collisions in PYTHIA 6 mode, especially when eCM < 10 GeV.
        if( .NOT.IS_PYTHIA8(i_mode) )then
!       Charge is not conserved. Re-generates the event.
            if( ABS(ps0(6)-ps1(6)) > 1D-10 ) goto 100
!       4-momentum is not conserved. Re-generates the event.
            do i=1,4,1
                if( ABS(ps0(i)-ps1(i)) > 1D-5 ) goto 100
            end do
!       Re-generates the event if any errors exist during the excution.
            ! if( MSTU(23) > 0 .OR. MSTU(30) > 0 ) goto 100
            if( MSTU(24) > 0 ) goto 100
!       Error due to the "iikk, kkii" in PYTHIA 6. Re-generates the event.
            if( iikk == 2 .OR. kkii == 2 ) goto 100
!       Re-generates the event if there are junctions when consider inelatic
!        processes in parton rescattering via SFM. Reduces complexity.
            if( ( i_had_model == 0 .OR. i_had_model == 3) &
                .AND. ( iparres == 1 .OR. iparres == 3 ) )then
                    if( MSTU(29) == 1 ) goto 100
            end if
!       Uses D-framework ?
            IF( .FALSE. )THEN
!       Re-generates the event if we consider the leading-particle
!        reconstructions (only in pA/Ap now). Reduces complexity.
!#TODO(Lei20240218): need new model for the initial multiple nucleon interaction
                if( ( ipden == 0 .AND. itden == 1 ) &
                    .OR. ( ipden == 1 .AND. itden == 0 ) )then
                    if( MSTU(29) == 1 ) goto 100
                end if
            END IF
        end if
!----------------------------   PYTHIA 6 Checking   ----------------------------
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
!-----------------------------   Data Recording   ------------------------------
!       Records the original contents of each collision after PYTHIA execution.
        ! In Pythia8_fort_interface.f90.
        call PARECO( i_NN_pair )

!       Removes unnecessary entries (mainly for speeding up "parcas").
!       Uses V(i,5) as index pointers for the further recovering after parcas.
!       0.1D0 is just a number used to ensure accurate conversion from INT().
        if( .NOT.IS_PYTHIA8(i_mode) )then
            do i=1,N,1
                if( mstptj == 0 ) V(i,5) = ( N_TOT_AA - N + i )*1D0 + 0.1D0
                ! Removes junctions for Coal, if any.
                if( i_had_model == 1 .AND. K(i,2) == 88 ) K(i,1) = 51
            end do
            ! In Pythia8_fort_interface.f90.
            call PAEDIT(1)
!       Moves non-historical from /PYJETS_PY8/ to /PYJETS/.
        else if( IS_PYTHIA8(i_mode) )then
            N = 0
            do i=1,N_PY8,1
                KS = K_PY8(i,1)
                if( .NOT.IS_EXIST(KS,i_mode) ) cycle
                N = N + 1
                do j=1,5,1
                    K(N,j) = K_PY8(i,j)
                    P(N,j) = P_PY8(i,j)
                    V(N,j) = V_PY8(i,j)
                end do
                if( mstptj == 0 ) V(N,5) = ( N_TOT_AA - N_PY8 + i )*1D0 + 0.1D0
            end do
        end if

!       Records the numbers of the participant nucleons (kpar),
!        binary pp (kpp), pn/np (knp), and nn collisions (knn),
!        or lepton-nucleon collisions (kep).
        if( (ipden < 2) .AND. (itden < 2) )then
            if( (kf_a == kf_b) .AND. (kf_a == 2212) ) kpp = kpp + 1
            if( (kf_a == kf_b) .AND. (kf_a == 2112) ) knn = knn + 1
            if(  kf_a /= kf_b  )                      knp = knp + 1
            if(  pT_a <= 1D-4 ) kpar = kpar + 1
            if(  pT_b <= 1D-4 ) kpar = kpar + 1
        else if( (ipden >= 11) .AND. (itden < 2) )then
            kep = kep + 1
        end if
        ! Deduction due to the regenerated event.
        if( i_NN_pair_0 == i_NN_pair )then
            if( (ipden < 2) .AND. (itden < 2) )then
                if( (kf_a == kf_b) .AND. (kf_a == 2212) ) kpp = kpp - 1
                if( (kf_a == kf_b) .AND. (kf_a == 2112) ) knn = knn - 1
                if(  kf_a /= kf_b  )                      knp = knp - 1
                if(  pT_a <= 1D-4 ) kpar = kpar - 1
                if(  pT_b <= 1D-4 ) kpar = kpar - 1
            else if( (ipden >= 11) .AND. (itden < 2) )then
                kep = kep - 1
            end if
        end if
!-----------------------------   Data Recording   ------------------------------
!-------------------------------------------------------------------------------


!       Collects lost 4-momentum.
        do i=1,4,1
            throe_p(i) = throe_p(i) + ( ps0(i) - ps1(i) )
            ! Deduction due to the regenerated event.
            if( i_NN_pair_0 == i_NN_pair ) throe_p(i) = throe_p(i) - dp0(i)
            dp0(i) = ps0(i) - ps1(i)
        end do

        i_NN_pair_0 = i_NN_pair

        ! Recovers.
        MSTJ(21)  = MSTJ21
        MSTP(111) = MSTP111


        return
        end



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine D_framework( iiii, naf00, i_a, i_b, i_error )
!!      Simulates D-framwork, i.e. does "parcas", "sfm/coal"
!!       NN pair by pari inside "parini".
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        INTEGER PYCHGE
        PARAMETER (KSZJ=80000)
        COMMON/PYDAT3/MDCY(500,3),MDME(8000,2),BRAT(8000),KFDP(8000,5)
        COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
        COMMON/PYSUBS/MSEL,MSUB(500),KFIN(2,-40:40),NON,CKIN(200)
        COMMON/PYPARS/MSTP(200),PARP(200),MSTI(200),PARI(200)
        COMMON/PYJETS/N,NPAD,K(KSZJ,5),P(KSZJ,5),V(KSZJ,5)
        COMMON/PYINT1/MINT(400),VINT(400)
        common/papr/t0,cspipiKK,dep,ddt,edipi,epin,ecsnn,ekn,ecspsn,ecspsm, &
         rnt,rnp,rao,rou0,vneu,vneum,ecsspn,ecsspm,ecsen
        common/syspar/ipden,itden,suppm,suptm,suppc,suptc,r0p,r0t, &
         nap,nat,nzp,nzt,pio
        common/wz/c17(500,3),ishp(kszj),tp(500),coor(3),p17(500,4)
        common/sa1/kjp21,non1,bp,iii,neve,nout,nosc
        common/sa2/nsa,non2,ksa(kszj,5),psa(kszj,5),vsa(kszj,5)
        common/sa4/tau(kszj),tlco(kszj,4)
        common/sa5/kfmax,kfaco(100),numb(100),numbs(100),non5, &
         disbe(100,100)
        common/sa6/kfmaxi,nwhole
        common/sa6_c/ithroq,ithrob,ich,non6_c,throe(4)
        common/sa6_p/ithroq_p,ithrob_p,ich_p,non6_p,throe_p(4)
        common/sa7/ispmax,isdmax,iflmax,ispkf(100),n_bin_hist,asd(10), &
         afl(100,10,2),i_y_or_eta
        common/sa10/csnn,cspin,cskn,cspipi,cspsn,cspsm,rcsit,ifram, &
         iabsb,iabsm,i_sigma_AQM,ajpsi,csspn,csspm,csen
        common/sa12/ppsa(5),nchan,nsjp,sjp,taup,taujp
        common/sa13/kjp20,non13,vjp20,vjp21,vjp22,vjp23
        common/sa15/nps,npsi,pps(5000,5),ppsi(5000,5)
        common/sa16/x_ratio,decpro,dni(10),dpi(10),edi(10),bmin,bmax, &
         bar(10),abar(10),barf(10),abarf(10),emin(10),eminf(10), &
         eplu(10),epluf(10),nmax
        common/sa21/pincl(5),pscal(5),pinch(5),vnu,fq2,w2l,yyl,zl,xb, &
          pph,vnlep
        common/sa23/kpar,knn,kpp,knp,kep,NON_p,skpar,sknn,skpp,sknp,skep
        common/sa24/adj1(40),nnstop,non24,zstop
        common/sa25/i_inel_proc,i_time_shower,para1_1,para1_2
        common/sa26/ndiq(kszj),npt(kszj),ifcom(kszj),idi,idio
        common/sa27/itime,kjp22,gtime,astr,akapa(6),parj1,parj2,parj3, &
         parj21,parj4,adiv,gpmax,nnc
        common/sa28/nstr,nstra(kszj),nstrv(kszj),nstr0, &
         nstr1,nstr1a(kszj),nstr1v(kszj)
        common/sa30/vneump,vneumt,mstptj
        common/sa33/smadel,ecce,secce,parecc,iparres
        common/sa34/itorw,iikk,cp0,cr0,kkii
        common/sbe/nbe,nonbe,kbe(kszj,5),pbe(kszj,5),vbe(kszj,5)
        common/saf/naf,nonaf,kaf(kszj,5),paf(kszj,5),vaf(kszj,5)
        common/aaff/naff,nonff,kaff(kszj,5),paff(kszj,5),vaff(kszj,5)
        common/sbh/nbh,nonbh,kbh(kszj,5),pbh(kszj,5),vbh(kszj,5)
        common/sa1_h/nn,non1_h,kn(kszj,5),pn(kszj,5),rn(kszj,5)
        common/ctllist/nctl,noinel(600),nctl0,nctlm
        common/ctllist_p/nreac(9),nrel
        common/sgam/ngam,nongam,kgam(kszj,5),pgam(kszj,5),vgam(kszj,5)
        common/coal1/bmrat,i_mm
!       For the evolution time.
        common/anly_time/ time_ini,  time_par,  time_had, &
                         stime_ini, stime_par, stime_had
!       For the calculation of tje single, multiple string tension.
        common/string_kapa1/ parj10, ww0, seco
!       For the independent Optical Clauber calculation.
        common/Glauber5/ pathn, evbin, pirr, tirr, spathn, sevbin
!       For the gamma related statistics.
        common/anly_gamma1/ egam,  egam1,  egam2,  egam3, &
                           segam, segam1, segam2, segam3
!       For the simulation control.
        COMMON/SA1_PY8/ i_mode, i_tune, KF_woDecay(100), &
               KF_proj, KF_targ, win, energy_B, psno, b_min, b_max
        dimension i_leading(2), Q_leading(2), Q_leading_0(2), quantity(10)
        data time_par_sum / 0D0 /
        data time_par_max / 0D0 /


        i_error = 0

!       Diffractive hadrons/leptons (without partons).
        if( N == 0 ) goto 998
        i_pure_hadron = 0


!-------------------------------------------------------------------------------
!--------------------------   Partonic Rescattering   --------------------------
        call do_partonic_rescattering( i_error, i_pure_hadron )
        if( i_error == 1 )then
            ! Counteracts the operation inside it.
            iii = iii + 1
            return
        end if
        ! Inaccurate in this mode.
        egam2 = 0D0
!--------------------------   Partonic Rescattering   --------------------------
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
!------------------------------   Hadronization  -------------------------------
        call do_hadronization( i_error, i_pure_hadron )
        if( i_error == 1 )then
            ! Counteracts the operation inside it.
            iii = iii + 1
            return
        end if

!------------------------   Failed Parton Collecting   -------------------------
!       Special treatment for D-framework using /sa6_c/.
        do i=1,nbe,1
            if( kbe(i,2) > 0 )then
                ithroq = ithroq + 1
            elseif( K(i,2) < 0 )then
                ithrob = ithrob + 1
            end if
            ich   = ich   + PYCHGE( kbe(i,2) )
            do j=1,4,1
                throe(j) = throe(j) + pbe(i,j)
            end do
        end do
        nbe = 0
!------------------------   Failed Parton Collecting   -------------------------

!------------------------------   Hadronization  -------------------------------
!-------------------------------------------------------------------------------


998     continue

!       "sbh" to "PYJETS". Hadrons from the diffractive event and remnants etc.
        if( nbh > 0 )then
            do jj=1,nbh,1
                ii = N + jj
                do m=1,5,1
                    K(ii,m) = kbh(jj,m)
                    P(ii,m) = pbh(jj,m)
                    V(ii,m) = vbh(jj,m)
                end do
            end do
            N = N + nbh
            nbh = 0
        end if


!-------------------------------------------------------------------------------
!------------------------   Leading Nucleon Filtering   ------------------------
        if( N > 0 )then
            ! Filters and orders particles.
            call filt
!           For pA, Ap and AA, returns nucleons only.
                ! Upto p.
                ! nup0 = numbs(1)
            ! Uptp p and n.
            nup0 = numbs(2)
!           "nup" is the number of leading particles.
            nup = nup0

            i_leading_mode = 0
            ! i_leading_mode:
            !       =0, all newly generated N with the arbitrary 4-momentum
            !       =1, one or two newly generated N with the min pT for
            !           pA/Ap (requires pz <0 / >0) or AB collisions
            !       =2, ... with the max +pz or/and -pz
            !       =3, ... with the max |p| ...
            !       =4, ... with the max E ...
            if( i_leading_mode < 0 .OR. i_leading_mode > 4 ) i_leading_mode = 1
            if( i_leading_mode /= 0 )then
                ! Index of the leading proton (nucleon).
                i_leading = 0
                ! Quantity used for filtering.
                Q_leading = 0D0
                Q_leading_0 = 0D0
                quantity = 0D0
                ! Finds the leading nucleon(s) "i_leading".
                do i=1,nup0,1
                    ! pT = PAPYP(i,10)
                    quantity(1) = PAPYP(i,10)
                    pz = P(i,3)
                    quantity(2) = ABS( pz )
                    ! abs_P = PAPYP(i,8)
                    quantity(3) = PAPYP(i,8)
                    ! ee = P(i,4)
                    quantity(4) = P(i,4)
                    ! ! pA, requires pz > 0.
                    ! if( ipden == 0 .AND. itden == 1 .AND. pz < 0D0 ) cycle
                    ! ! Ap, requires pz < 0.
                    ! if( ipden == 0 .AND. itden == 0 .AND. pz > 0D0 ) cycle
                    if( pz > 0D0 )then
                        Q_leading(1) = quantity( i_leading_mode )
                        ! Min pT.
                        if( i_leading_mode == 1 )then
                            if( Q_leading(1) < Q_leading_0(1) )then
                                i_leading(1) = i
                                Q_leading_0(1) = Q_leading(1)
                            end if
                        ! Max pz, |p| or E.
                        else
                            if( Q_leading(1) > Q_leading_0(1) )then
                                i_leading(1) = i
                                Q_leading_0(1) = Q_leading(1)
                            end if
                        end if
                    else if( pz < 0D0 )then
                        Q_leading(2) = quantity( i_leading_mode )
                        ! Min pT.
                        if( i_leading_mode == 1 )then
                            if( Q_leading(2) < Q_leading_0(2) )then
                                i_leading(2) = i
                                Q_leading_0(2) = Q_leading(2)
                            end if
                        ! Max pz, |p| or E.
                        else
                            if( Q_leading(2) > Q_leading_0(2) )then
                                i_leading(2) = i
                                Q_leading_0(2) = Q_leading(2)
                            end if
                        end if
                    end if
                end do
                ! Fills the leading nucleon(s) into the 1-st (2-nd) entry.
                nup = 0
                ! The leading nucleon with +pz.
                if( i_leading(1) > 0 )then
                    nup = nup + 1
                    i_exchange = i_leading(1)
                    do j=1,5,1
                        ! nup -> N+1
                        K( N+1, j ) = K( nup, j )
                        P( N+1, j ) = P( nup, j )
                        V( N+1, j ) = V( nup, j )
                        ! i_exchange (leading nucleon) -> nup
                        K( nup, j ) = K( i_exchange, j )
                        P( nup, j ) = P( i_exchange, j )
                        V( nup, j ) = V( i_exchange, j )
                        ! N+1 -> i_exchange
                        K( i_exchange, j ) = K( N+1, j )
                        P( i_exchange, j ) = P( N+1, j )
                        V( i_exchange, j ) = V( N+1, j )
                    end do
                    ! Let's set its position as the old collided nucleon.
                    i_position = i_a
                    if( psa( i_b, 3 ) > 0D0 ) i_position = i_b
                    do j=1,3,1
                        V( nup, j ) = vsa( i_position, j )
                    end do
                    ! Avoids two new leading nucleons colliding with each other.
                    if( ipden == 1 .AND. itden == 1 )then
                        do j=1,3,1
                            V( nup, j ) = V( nup, j ) + P(nup,j)/P(nup,4) * ddt
                        end do
                    end if
                end if
                ! The leading nucleon with -pz.
                if( i_leading(2) > 0 )then
                    nup = nup + 1
                    i_exchange = i_leading(2)
                    do j=1,5,1
                        ! nup -> N+1
                        K( N+1, j ) = K( nup, j )
                        P( N+1, j ) = P( nup, j )
                        V( N+1, j ) = V( nup, j )
                        ! i_exchange (leading nucleon) -> nup
                        K( nup, j ) = K( i_exchange, j )
                        P( nup, j ) = P( i_exchange, j )
                        V( nup, j ) = V( i_exchange, j )
                        ! N+1 -> i_exchange
                        K( i_exchange, j ) = K( N+1, j )
                        P( i_exchange, j ) = P( N+1, j )
                        V( i_exchange, j ) = V( N+1, j )
                    end do
                    ! Let's set its position as the one of collided nucleon.
                    i_position = i_a
                    if( psa( i_b, 3 ) < 0D0 ) i_position = i_b
                    do j=1,3,1
                        V( nup, j ) = vsa( i_position, j )
                    end do
                    ! Avoids two new leading nucleons colliding with each other.
                    if( ipden == 1 .AND. itden == 1 )then
                        do j=1,3,1
                            V( nup, j ) = V( nup, j ) + P(nup,j)/P(nup,4) * ddt
                        end do
                    end if
                end if
            end if
!           "nup" is the number of leading particles.
            naf = naf00
!           Energy cut.
            call PASORT( 1, nup, "pyjets", 'e', "min_to_max" )
            nup0 = nup
            ! nup  = 0
            E_cut = 0D0
            ! E_cut = 10D0
            ! do i=1,nup0,1
            !     if( P(i,4) > E_cut ) nup = nup + 1
            ! end do
!           "PYJETS" above "nup" to "saf".
            do i = nup+1, N, 1
                naf = naf + 1
                do j=1,5,1
                    kaf(naf,j) = K(i,j)
                    paf(naf,j) = P(i,j)
                    vaf(naf,j) = V(i,j)
                end do
            end do
!           "PYJETS" below "nup" to "sbh".
            if( nup > 0 )then
                do j=1,5,1
                    do i=1,nup,1
                        kbh(i,j) = K(i,j)
                        pbh(i,j) = P(i,j)
                        vbh(i,j) = V(i,j)
                    end do
                end do
            end if
            nbh = nup
            nup = 0
            N = 0
        end if
!------------------------   Leading Nucleon Filtering   ------------------------
!-------------------------------------------------------------------------------


!       Defineds the time of parcas as the averaged one in D-framework.
        time_par_sum = time_par_sum + time_par
        time_par = time_par_sum / iiii * 1D0
        ! Or the max time?
        ! if( time_par < time_par_max ) time_par = time_par_max
        ! time_par_max = time_par


        return
        end



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine chargecme
!!      The CME-induced charge initial charge separation by switching the
!!       py values of a fraction of the downward(upward) moving(u,d,s,c)quarks
!!       for symmetrical collision systems,i.e., Ru&Ru Zr&Zr at RHIC and LHC.
!       Here in symmetrical systems, nap=nat,nzp=nzt, and the fraction and
!        magnetic field function is A*bp-B*bp^3 type.  by shezl 2021
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        PARAMETER (KSZJ=80000)
        LOGICAL IS_EXIST
        COMMON/PYJETS/N,NONJ,K(KSZJ,5),P(KSZJ,5),V(KSZJ,5)
        common/sa1/kjp21,non1,bp,iii,neve,nout,nosc
        common/syspar/ipden,itden,suppm,suptm,suppc,suptc,r0p,r0t, &
         nap,nat,nzp,nzt,pio
        common/schuds/schun,schudn,schudsn,sfra
        common/sa24/adj1(40),nnstop,non24,zstop
        COMMON/SA1_PY8/ i_mode, i_tune, KF_woDecay(100), &
               KF_proj, KF_targ, win, energy_B, psno, b_min, b_max
        common/sa6_p/ithroq_p,ithrob_p,ich_p,non6_p,throe_p(4)
        common/cme_local/ numk(kszj)
        real(kind=8) :: p2u,erhic,erela,ruzcp


        ! RHIC energy 200.
        erhic = 200D0
        ! RHIC energy as a base.
        erela = 0.45D0 + 0.55D0 * ( win / erhic )
        rnzp  = DBLE(nzp)
        rnap  = DBLE(nap)
        rerz  = rnzp/rnap
        ! isobar Zr Ru(96,42) as a base.
        ruzcp = ( ( 96D0 / 42D0 ) * rerz )**(0.667D0)

        sfra  = 3.1D0 * ( 2448.135D0 * nap**(-1.667D0) * bp &
              - 160.810D0 * nap**(-2.333D0) * bp**3 ) * erela * ruzcp * 0.01D0

        do i=1,N,1
!           Only for the particles that still exist.
            KSi = K(i,1)
            if( .NOT.IS_EXIST(KSi,i_mode) ) cycle
            if( ABS( K(i,2) ) == 1 .OR. ABS( K(i,2) ) == 2 .OR. &
                ABS( K(i,2) ) == 3 .OR. ABS( K(i,2) ) == 4        )then
                schun = schun + 1D0
                if( PYR(1) <= sfra )then
                    numk(i) = 0
                    schudn = schudn + 1D0
                    do ii=1,N,1
!                   Only for the particles that still exist.
                        KSii = K(ii,1)
                        if( .NOT.IS_EXIST(KSii,i_mode) ) cycle
                        if( numk(ii) == 1 ) cycle
                        do jj=ii+1,N,1
!                       Only for the particles that still exist.
                            KSjj = K(jj,1)
                            if( .NOT.IS_EXIST(KSjj,i_mode) ) cycle
                            if( numk(jj) == 1 ) cycle
                            if( ( K(ii,2) + K(jj,2) ) == 0 .AND. &
                                ( K(ii,2) * P(ii,2) < 0 )  .AND. &
                                ( K(jj,2) * P(jj,2) < 0 )          )then
                                p2u = P(ii,2)
                                P(ii,2)  = P(jj,2)
                                P(jj,2)  = p2u
                                schudsn  = schudsn + 1D0
                                numk(ii) = 1
                                numk(jj) = 1
!       Recalculates energies to ensure on-shell then collects the lost 4-mom.
                                eii = p(ii,4)
                                ejj = p(jj,4)
                                p(ii,4) = SQRT( p(ii,5)**2 + p(ii,1)**2 &
                                        +       p(ii,2)**2 + p(ii,3)**2 )
                                p(jj,4) = SQRT( p(jj,5)**2 + p(jj,1)**2 &
                                        +       p(jj,2)**2 + p(jj,3)**2 )
                                throe_p(4) = throe_p(4) &
                                           + ( eii + ejj - p(ii,4) - p(jj,4) )
                            end if
                        end do
                    end do
                end if
            end if
        end do


        return
        end



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine lorntz( ilo, b, pi, pj )
!!      Performs Lorentz (or inverse Lorentz) transformation.
!       ilo: =0, Lorentz transformation.
!            =1, inverse Lorentz transformation.
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        dimension pi(4), pj(4), b(3), dpi(4), dpj(4)


        bb = b(1)*b(1) + b(2)*b(2) + b(3)*b(3)
        DB = SQRT(bb)
        eps1 = 1D0 - 1D-12
        if( DB > eps1 )then
            do i=1,3
!       Rescales boost vector if too close to unity.
                b(i) = b(i)*( eps1 / DB )
            end do
            DB = eps1
            bb = DB**2
        end if
        bbb = 1D0 - bb
        gam = 1D0 / SQRT(bbb)
        ga  = gam*gam / ( gam + 1D0 )
        do i=1,4,1
            dpi(i) = pi(i)
            dpj(i) = pj(i)
        end do

        if( ilo == 0 )then
!       Lorentz transformation.
            pib = dpi(1)*b(1) + dpi(2)*b(2) + dpi(3)*b(3)
            pjb = dpj(1)*b(1) + dpj(2)*b(2) + dpj(3)*b(3)
            do i=1,3
                pi(i) = dpi(i) + b(i)*( ga*pib - gam*dpi(4) )
                pj(i) = dpj(i) + b(i)*( ga*pjb - gam*dpj(4) )
            end do
            pi(4) = gam*( dpi(4) - pib )
            pj(4) = gam*( dpj(4) - pjb )
        else
!       inverse Lorentz transformation
            pib = dpi(1)*b(1) + dpi(2)*b(2) + dpi(3)*b(3)
            pjb = dpj(1)*b(1) + dpj(2)*b(2) + dpj(3)*b(3)
            do i=1,3
                pi(i) = dpi(i) + b(i)*( ga*pib + gam*dpi(4) )
                pj(i) = dpj(i) + b(i)*( ga*pjb + gam*dpj(4) )
            end do
            pi(4) = gam*( dpi(4) + pib )
            pj(4) = gam*( dpj(4) + pjb )
        end if


        return
        end



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine his(t1,istop)
!!      Classical Newton motion in cms system of colliding pair.
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        PARAMETER (KSZJ=80000)
        PARAMETER (NSIZE=750000)
        common/sa2/nsa,non2,ksa(kszj,5),psa(kszj,5),vsa(kszj,5)
        common/sa4/tau(kszj),tlco(kszj,4)
        common/papr/t0,cspipiKK,dep,ddt,edipi,epin,ecsnn,ekn,ecspsn,ecspsm, &
         rnt,rnp,rao,rou0,vneu,vneum,ecsspn,ecsspm,ecsen
        common/ctllist/nctl,noinel(600),nctl0,nctlm
        common/ctllist_t/ lc(nsize,5), tc(nsize), tw(nsize)
        common/wz/c17(500,3),ishp(kszj),tp(500),coor(3),p17(500,4)


        istop = 1
        i_n   = 0
        r0 = rao * dmax1( rnt, rnp )
        do 200 i=1,nsa,1
!           if(t1 <= tau(i)) goto 100
!       do move particles which have not produced
            if( ishp(i) == 1 ) goto 10
!           pp4 = psa(i,4)
!           do j=1,3
!               vp = psa(i,j) / pp4
!               vsa(i,j) = vsa(i,j) + vp * (t1 - vsa(i,4) )
!           end do
            i_n = i_n + 1
            goto 200

10          aa  = 0D0
            pp4 = psa(i,4)
!       Due to the fast speed of baryons, we could not use a limited interaction
!        region.
!           if( ABS( K(i,2)) > 1000 ) r0 = 1D0 + 10D0 * r0
            do j=1,3,1
                vp = psa(i,j) / pp4
                vsa(i,j) = vsa(i,j) + vp * ( t1 - vsa(i,4) )
                aa = aa + ( vsa(i,j) - coor(j) )**2
            end do
!           vsa(i,4)=t1
            aa = SQRT(aa)
            if( aa < r0 ) goto 100
!       If freeze-out occurs, deducts the distance between the last collision
!        and now
            do j=1,3,1
                vp = psa(i,j) / pp4
                vsa(i,j) = vsa(i,j) - vp * ( t1 - vsa(i,4) )
            end do
            ishp(i) = 0
            do il=1,nctl,1
                if( lc(il,1) == i .OR. lc(il,2) == i ) tc(il) = 0D0
            end do
            goto 200
100         continue
            vsa(i,4) = t1
200     continue
        if( i_n == nsa ) return
        istop = 0


        return
        end



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine ptcre(l,l1,time)
!!      Gives four position to the particles after calling PYTHIA.
!       l and l1 are colliding particles.
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        LOGICAL IS_EXIST,IS_NUCLEUS
        PARAMETER (KSZJ=80000)
        COMMON/PYJETS/N,NPAD,K(KSZJ,5),P(KSZJ,5),V(KSZJ,5)
        common/sa2/nsa,non2,ksa(kszj,5),psa(kszj,5),vsa(kszj,5)
        common/sa4/tau(kszj),tlco(kszj,4)
        common/syspar/ipden,itden,suppm,suptm,suppc,suptc,r0p,r0t, &
         nap,nat,nzp,nzt,pio
        COMMON/SA1_PY8/ i_mode, i_tune, KF_woDecay(100), &
               KF_proj, KF_targ, win, energy_B, psno, b_min, b_max


        do i=1,N,1
!       Only for the particles that still exist.
            KS = K(i,1)
            if( .NOT.IS_EXIST(KS,i_mode) ) cycle
!       From hh/lh collisions, generated particles  are distributed on the
!        surface with unit radius.
            if( .NOT.IS_NUCLEUS(KF_proj) .AND. .NOT.IS_NUCLEUS(KF_targ) )then
                COS_theta = 2D0*PYR(1) - 1D0
                SIN_theta = SQRT( 1D0 - COS_theta**2 )
                phi = 2D0 * pio * PYR(1)
                v(i,1) = SIN_theta * COS(phi)
                v(i,2) = SIN_theta * SIN(phi)
                v(i,3) = COS_theta
!       For hA/lA/AB, between two collided particles.
            else
                rl = PYR(1)
                do m=1,3,1
                    v(i,m) = v(i,m) + vsa(l,m) * rl + vsa(l1,m) * ( 1D0 - rl )
                end do
            end if
            V(i,4) = time
        end do


        return
        end



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine ptcre_n(time,gamt)
!!      Arranges particles (quark,diquark, and gluon mainly) after
!        calling PYTHIA into the overlap region randomly.
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        PARAMETER (KSZJ=80000)
        COMMON/PYJETS/N,NPAD,K(KSZJ,5),P(KSZJ,5),V(KSZJ,5)
        common/sa1/kjp21,non1,bp,iii,neve,nout,nosc
        common/syspar/ipden,itden,suppm,suptm,suppc,suptc,r0p,r0t, &
         nap,nat,nzp,nzt,pio


        if(ipden /= 0 .or. itden /= 0)then
            b=bp/r0t
            do i=1,N,1
                iiii=0
54              iiii=iiii+1
                if(iiii == 10000)then
                    write(22,*) 'difficult to arrange produced ' // &
                                 'particles in subroutine ptcre,' // &
                                 'infinitive loop may occur'
                endif
!       sample a point in the unit sphere of target
                x=1.-2.*pyr(1)
                y=1.-2.*pyr(1)
                z=1.-2.*pyr(1)
                rr=dsqrt(x*x+y*y+z*z)
                if(rr > 1) goto 54
!       x and y components of that point in the system of unit sphere of
!        projectile are x and y-b, respectively. Adjudge that does (x,y-b) is
!        in the sphere of projectile
                r1=r0p*dsqrt(x*x+(b-y)*(b-y))
                if(r1 > r0p) goto 54
                xx=x*r0t
                yy=y*r0t
                zz=z*r0t/gamt
                v(i,1)=xx
                v(i,2)=yy
                v(i,3)=zz
            enddo
        endif
        if(ipden == 0 .and. itden == 0)then
            do i=1,N,1
                cita=2*pyr(1)-1.
                fi=2.*pio*pyr(1)
                sita=dsqrt(1.-cita**2)
                v(i,1)=sita*dcos(fi)
                v(i,2)=sita*dsin(fi)
                v(i,3)=cita
                v(i,3)=v(i,3)/gamt
            enddo
        endif
        do i=1,N,1
            v(i,4) = time
        enddo


        return
        end



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine updpip(l,l1,time)
!!      update hadron list 'sa2' after calling PYTHIA ('sbh' to 'sa2')
!!       remove collision pair composed of l and/or l1, remove l (l1)
!!       from 'sa2'.
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        PARAMETER (KSZJ=80000)
        PARAMETER (NSIZE=750000)
        COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
        common/wz/c17(500,3),ishp(kszj),tp(500),coor(3),p17(500,4)
        common/sbh/nbh,nonbh,kbh(kszj,5),pbh(kszj,5),vbh(kszj,5)
        common/ctllist/nctl,noinel(600),nctl0,nctlm
        common/ctllist_t/ lc(nsize,5), tc(nsize), tw(nsize)
        common/sa1/kjp21,non1,bp,iii,neve,nout,nosc
        common/sa2/nsa,non2,ksa(kszj,5),psa(kszj,5),vsa(kszj,5)
        common/sa4/tau(kszj),tlco(kszj,4)
        common/sa5/kfmax,kfaco(100),numb(100),numbs(100),non5, &
         disbe(100,100)
        common/sa6/kfmaxi,nwhole
        common/sa12/ppsa(5),nchan,nsjp,sjp,taup,taujp
        common/sa14/ipyth(2000),idec(2000),iwide(2000)
        common/sa30/vneump,vneumt,mstptj
        common/papr/t0,cspipiKK,dep,ddt,edipi,epin,ecsnn,ekn,ecspsn,ecspsm, &
         rnt,rnp,rao,rou0,vneu,vneum,ecsspn,ecsspm,ecsen
        common/syspar/ipden,itden,suppm,suptm,suppc,suptc,r0p,r0t, &
         nap,nat,nzp,nzt,pio
!       ipyth: stord line number of produced hadron in hadron list (sa2)


        do m=1,2000
            ipyth(m)=0
        enddo
!       'sbh' to 'sa2' (i.e. produced hadrons-> hadron list 'sa2')
        ll  = l
        ll1 = l1
        if(nbh == 0) goto 200
        do 500 i=1,nbh
            kf = kbh(i,2)
            do 600 j=1,kfmax,1
                if(kf /= kfaco(j)) goto 600
                jj = numb(j) + 1
!       update particle list etc.
                do m=nsa,jj,-1
                    mm = m + 1
!080104             ksa(mm,2) = ksa(m,2)
!080104             ksa(mm,1) = 1
!080104             ksa(mm,3) = ksa(m,3)
                    do m1=1,5
                        ksa(mm,m1) = ksa(m,m1)   ! 080104
                        psa(mm,m1) = psa(m,m1)
                        vsa(mm,m1) = vsa(m,m1)
                    enddo
                    ishp(mm) = ishp(m)
                    tau(mm)  = tau(m)
                enddo
                do m=1,2000
                    ipym = ipyth(m)
                    if(ipym >= jj) ipyth(m) = ipym + 1
                enddo
                if(ll  >= jj) ll  = ll  + 1
                if(ll1 >= jj) ll1 = ll1 + 1
!       update the values of lc(m,1-2) with value >= jj
                do m=1,nctl
                    lc1 = lc(m,1)
                    if(lc1 >= jj) lc(m,1) = lc1 + 1
                    lc2 = lc(m,2)
                    if(lc2 >= jj) lc(m,2) = lc2 + 1
                enddo
!       give proper values to particle jj.
!221203         ksa(jj,2) = kf
!221203         ksa(jj,1) = 1
!221203         ksa(jj,3) = 0
                do m=1,5
                    ksa(jj,m) = kbh(i,m)   ! 221203
                    psa(jj,m) = pbh(i,m)
                    vsa(jj,m) = vbh(i,m)
                enddo
                ishp(jj) = 1
                tau(jj) = time + t0 * pbh(i,4)/pbh(i,5)
!       the values of 'ishp' and 'tau' for hadrons from 'PYTHIA'
!        are given here, the proper formation time of 'PYTHIA' particle
!        is assume to be equal to t0 fm/c, except nucleon and j/psi
                if(kf == 2212 .or. kf == 2112)then
                    tau(jj) = time + t0 * pbh(i,4)/pbh(i,5) * taup
                elseif(kf == 443.or.kf == 100443)then
                    tau(jj) = time + t0 * pbh(i,4)/pbh(i,5) * taujp
                endif
                ipyth(i) = jj
                do m=j,kfmax
                    numb(m) = numb(m) + 1
                enddo
                nsa = nsa + 1
                goto 500
600         enddo   ! 040223
!040223 if produced hadron is not in given hadron classification
            nsa = nsa + 1
            do m=1,5
                ksa(nsa,m) = kbh(i,m)
                psa(nsa,m) = pbh(i,m)
                vsa(nsa,m) = vbh(i,m)
            enddo
            ishp(nsa) = 0
            tau(nsa)  = 0D0
            ipyth(i)  = nsa
500     enddo   ! 040223
200     continue   ! 241110
        l  = ll
        l1 = ll1
!       remove colli. pair composed of l or l1
        jj = 0
        do 300 ii=1,nctl
            i1 = lc(ii,1)
            j1 = lc(ii,2)
            if(i1 == l .or. i1 == l1) goto 300
            if(j1 == l .or. j1 == l1) goto 300
            jj = jj + 1
            tc(jj) = tc(ii)
            tw(jj) = tw(ii)
            do m=1,5
                lc(jj,m) = lc(ii,m)
            enddo
300     continue
        do ii = jj+1, nctl+1
            tc(ii) = 0.0
            tw(ii) = 0.0
            do m=1,5
                lc(ii,m) = 0
            enddo
        enddo
        nctl = jj
!       remove hadrons l and l1 from 'sa2'
        kf1 = ksa(l,2)
        kf2 = ksa(l1,2)
        kf = kf1
        ll = l
        do 700 i=1,2
            if(ll == nsa)then   !
                do i1=1,kfmax
                    if(kf /= kfaco(i1)) goto 400
!241110             numbm = numb(i1)
!241110             do i2=1,i1
!241110                 if(numb(i2) == numbm) numb(i2) = numb(i2) - 1
!241110             enddo
                    do m=i1,kfmax
                        numb(m) = numb(m) - 1
                    enddo
                    if(i1 > 1)then
                        numba = numb(i1)
                        do m=1,i1-1
                            if(numb(m) == numba) numb(m) = numb(m) - 1
                        enddo
                    endif
!241110
                    goto 100
400             enddo
            endif   !
            do j=ll+1,nsa
                jj = j - 1
!080504         ksa(jj,2) = ksa(j,2)
!080504         ksa(jj,1) = 1
!080504         ksa(jj,3) = ksa(j,3)
                do m=1,5
                    ksa(jj,m) = ksa(j,m)   ! 080504
                    psa(jj,m) = psa(j,m)
                    vsa(jj,m) = vsa(j,m)
                enddo
                ishp(jj) = ishp(j)
                tau(jj)  = tau(j)
            enddo
            if(nctl == 0) goto 900
            do m=1,nctl
                lc1 = lc(m,1)
                lc2 = lc(m,2)
                if(lc1 > ll) lc(m,1) = lc1 - 1
                if(lc2 > ll) lc(m,2) = lc2 - 1
            enddo
900         do 800 j=1,kfmax
                if(kf /= kfaco(j)) goto 800
                do m=j,kfmax
                    numb(m) = numb(m) - 1
                enddo
                if(j > 1)then
                    numba = numb(j)
                    do m=1,j-1
                        if(numb(m) == numba) numb(m) = numb(m) - 1
                    enddo
                endif
                goto 100
800         continue
100         continue
            nsa = nsa - 1
            if(l1 > ll) l1 = l1 - 1
            do m=1,2000
                ipym = ipyth(m)
                if(ipym > ll) ipyth(m) = ipym - 1
            enddo
            if(i == 2) goto 700
            ll = l1
            kf = kf2
700     continue


        return
        end



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine coelas(ic,jc,eij,pi,pj)
!!      Performs elastic scattering.
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        PARAMETER (KSZJ=80000)
        common/sa2/nsa,non2,ksa(kszj,5),psa(kszj,5),vsa(kszj,5)
        dimension pi(4),pj(4)


        iic = ksa(ic,2)
        jjc = ksa(jc,2)
        d   = 3.65D0 * ( eij - PYMASS(iic) - PYMASS(jjc) )
        if( d < 1D-10 ) return
        ! pt = 0.2D0
        pt = 0.5D0
        A  = MIN( 10.3D0, 1.D0 / ( 1.12D0 * pt ) / ( 1.12D0 * pt ) )
        d6 = d**6
        B  = d6 * A / ( 1D0 + d6 )
        if( B < 1D-20 ) B = 1.D-20
        pm2 = pi(1)**2 + pi(2)**2 + pi(3)**2
        pm  = SQRT(pm2)
        t0  = - 4D0 * pm2
        if( ABS(t0) < 1.D-20 )then
            cctas = 1.D0
            goto 100
        end if
        cc = PYR(1)
        if( ABS( B*t0 ) < 0.0001D0 )then
            abt = 1.D0
!       els eif( B*t0 < -50.D0 )then
!           abt = 0.D0
        else
            abt = EXP( MAX( -7.0D2, B*t0 ) )
        end if
        tt1 = LOG( cc + (1D0 - cc) * abt )
        if( ABS(tt1) < 1.D-30 .AND. B <= 1.D-20 )then
            cctas = 1.D0
            goto 100
        end if
        tt = tt1 / b
        if( ABS(tt) < 1.D-20 )then
            cctas = 1.D0
            goto 100
        end if
        cctas = 1.D0 - tt * 2.D0 / t0
        if( ABS(cctas) > 1.D0 )then
            cctas = SIGN( 1.D0, cctas )
        end if
100     continue
        sctas = SQRT( 1.D0 - cctas**2 )
        fis   = 2.D0 * 3.141592653589793D0 * PYR(1)
        cfis  = COS(fis)
        sfis  = SIN(fis)
        call rotate( cctas, sctas, cfis, sfis, pm, pi, pj )


        return
        end



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine rotate( cctas, sctas, cfis, sfis, pp3, pi, pj )
!!      Performs rotation.
!       pi,pj: input, four momentum of colliding pair before scattering
!              output,four momentum of scattered particles after rotation
!       pp3: momentum modulus of pi or pj, both are equal in their cms,
!        after scattering
!       cctas, sctas, cfis, sfis: direction cosines/sines of the momentum of
!        one of scattered particles relative to the momentum of the
!        corresponding particle before scattering.
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        dimension pi(4), pj(4)


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
        subroutine updple(ic,jc,b,pi,pj)
!!      Updates particle list for elastic scattering.
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        PARAMETER (KSZJ=80000)
        common/sa2/nsa,non2,ksa(kszj,5),psa(kszj,5),vsa(kszj,5)
!       note the name of the arrays in 'sa2'
        dimension pi(4),pj(4),b(3)


        ilo=1
!       ilo=1 for inverse Lorentz transformation
        call lorntz(ilo,b,pi,pj)
        do i=1,4
            psa(ic,i)=pi(i)
            psa(jc,i)=pj(i)
        enddo


        return
        end



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine tran_saf
!!      'saf' to 'PYJETS'
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        PARAMETER (KSZJ=80000)
        COMMON/PYJETS/N,NPAD,K(KSZJ,5),P(KSZJ,5),V(KSZJ,5)
        common/saf/nsa,nonsa,ksa(kszj,5),psa(kszj,5),vsa(kszj,5)


        do m=1,5
            do l=1,nsa
                k(l,m)=ksa(l,m)
                p(l,m)=psa(l,m)
                v(l,m)=vsa(l,m)
            enddo
        enddo
        n=nsa


        return
        end



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine tran_sbe
!!      'sbe' to 'PYJETS'
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        PARAMETER (KSZJ=80000)
        COMMON/PYJETS/N,NPAD,K(KSZJ,5),P(KSZJ,5),V(KSZJ,5)
        common/sbe/nbe,nonbe,kbe(kszj,5),pbe(kszj,5),vbe(kszj,5)


        do m=1,5
            do l=1,nbe
                k(l,m)=kbe(l,m)
                p(l,m)=pbe(l,m)
                v(l,m)=vbe(l,m)
            enddo
        enddo
        n=nbe


        return
        end



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine tran_sa2   ! 280524
!!      'sa2' to 'PYJETS'
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        PARAMETER (KSZJ=80000)
        COMMON/PYJETS/N,NPAD,K(KSZJ,5),P(KSZJ,5),V(KSZJ,5)
        common/sa2/nbe,nonbe,kbe(kszj,5),pbe(kszj,5),vbe(kszj,5)


        do m=1,5
            do l=1,nbe
                k(l,m)=kbe(l,m)
                p(l,m)=pbe(l,m)
                v(l,m)=vbe(l,m)
            enddo
        enddo
        n=nbe


        return
        end



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine conse(np,pp,ps,jj)
!!      keep four momentum conservation
!       np : the # of particles
!       ps : four momentum to which the four momenta of particles should
!             conserve
!       pp : four momenta of particles
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
        common/papr/t0,cspipiKK,dep,ddt,edipi,epin,ecsnn,ekn,ecspsn,ecspsm, &
         rnt,rnp,rao,rou0,vneu,vneum,ecsspn,ecsspm,ecsen
        dimension pp(250,5),ps(4),ff(250),pxyz(3),arp(3)


        ps4=ps(4)
        do i=1,3
            pxyz(i)=0.
        enddo
        jj=0
100     es=0.
        do i=1,np
            es=es+pp(i,4)
        enddo
        fr=es/ps4
        if(dabs(1.-fr)  <=  dep) goto 200
        do i=1,np
            ppm=pp(i,4)/0.938
            ppf=ppm/fr
            ff(i)=dsqrt(dabs(ppf*ppf-1.)/(ppm*ppm-1.))
            do j=1,3
                ppp=ff(i)*pp(i,j)
                pp(i,j)=ppp
                pxyz(j)=pxyz(j)+ppp
            enddo
        enddo
        do i=1,3
            arp(i)=dabs(1.-pxyz(i)/ps(i))
            pxyz(i)=pxyz(i)-ps(i)
        enddo
        if(dabs(1.-fr) <= dep .and.arp(1) <= dep .and. arp(2) <= dep &
         .and. arp(3) <= dep) goto 200
        do i=1,3
            pxyz(i)=pxyz(i)/np
        enddo
        do i=1,np
            do j=1,3
                pp(i,j)=pp(i,j)-pxyz(j)
            enddo
            pp(i,4)=dsqrt(0.880+pp(i,1)**2+pp(i,2)**2+pp(i,3)**2)
!       0.880 = 0.938*0.938
        enddo
        jj=jj+1
        if(jj == 4000)then
            write(9,*) 'infinitive loop may occur in subroutine ' // &
                       'conse(), which means four-momentum ' // &
                       'conservation needed is hard to be achieved,' // &
                       'check value of PARAM(9)'
            return
        endif
        goto 100


200     return
        end



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine codi(pis,cfi1,sfi1,ccta1,scta1)
!!      Calculates the 'orientation' of the vector pis.
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        dimension pis(4),pi(4)


        pi    = pis
        fi1s  = PYANGL( pis(1), pis(2) )
        cta1s = PYANGL( pis(3), SQRT( pis(1)**2 + pis(2)**2) )
        fi1   = fi1s
        cta1  = cta1s
        cfi1  = COS(fi1)
        sfi1  = SIN(fi1)
        ccta1 = COS(cta1)
        scta1 = SIN(cta1)


        return
        end



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine rosa( cfi1,sfi1, ccta1,scta1, cfis,sfis, cctas,sctas, pis)
!!      Performs rotatation for produced particles from 'PYTHIA'.
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        dimension pis(4),pi(4)


        pi = pis
        pp = SQRT( pi(1)*pi(1) + pi(2)*pi(2) + pi(3)*pi(3) )
        call codi( pis, cfis, sfis, cctas, sctas )
        pi(1) = pp * ( cfi1*( ccta1*sctas*cfis +scta1*cctas) - sfi1*sctas*sfis )
        pi(2) = pp * ( sfi1*( ccta1*sctas*cfis +scta1*cctas) + cfi1*sctas*sfis )
        pi(3) = pp * ( ccta1*cctas - scta1*sctas*cfis )
        pis = pi


        return
        end



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine remo
!!      Moves hadron (lepton) and junction from "PYJETS" to "sbh".
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        LOGICAL IS_EXIST,IS_PARTON
        PARAMETER (KSZJ=80000)
        COMMON/PYJETS/N,NPAD,K(KSZJ,5),P(KSZJ,5),V(KSZJ,5)
        common/sbh/nbh,nonbh,kbh(kszj,5),pbh(kszj,5),vbh(kszj,5)
        common/sa1/kjp21,non1,bp,iii,neve,nout,nosc
        common/sa12/ppsa(5),nchan,nsjp,sjp,taup,taujp
        common/sa24/adj1(40),nnstop,non24,zstop
        COMMON/SA1_PY8/ i_mode, i_tune, KF_woDecay(100), &
               KF_proj, KF_targ, win, energy_B, psno, b_min, b_max


        nbh = 0
!       Hadronization model.
        i_had_model = INT( adj1(12) )
        jb = 0
201     continue
        do i1 = jb+1, N ,1
            KS = K(i1,1)
            KF = K(i1,2)
            if( i_had_model == 0 .OR. i_had_model == 3 )then
                ! Keeps junctions in "PYJETS" for sfm.
                if( IS_PARTON(KF) .OR. KF == 88 &
                    .OR. .NOT.IS_EXIST(KS,i_mode)  ) then
                    jb = jb + 1
                    cycle
                end if
            else if( i_had_model /= 0 .AND. i_had_model /= 3 )then
                ! Removes junctions from "PYJETS" for coal.
                if( IS_PARTON(KF) .OR. .NOT.IS_EXIST(KS,i_mode) ) then
                    jb = jb + 1
                    cycle
                end if
            end if
            ! Throws away junctions for Coal.
            if( KF /= 88 )then
                nbh = nbh + 1
                do i2=1,5,1
                    kbh(nbh,i2) = K(i1,i2)
                    pbh(nbh,i2) = P(i1,i2)
                    vbh(nbh,i2) = V(i1,i2)
                end do
            end if
            ! Moves the particle list one step downward from i1+1 to N.
            do jj=1,5,1
                do j = i1+1,N, 1
                    k(j-1,jj) = K(j,jj)
                    p(j-1,jj) = P(j,jj)
                    v(j-1,jj) = V(j,jj)
                end do
            end do
            N = N - 1
            goto 201
        end do


        return
        end



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine break
!!      Breaks up diquark (anti-diquark), gives four momenta
!!       and four positions to the broken objects
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        LOGICAL IS_EXIST,IS_DIQUARK
        PARAMETER (KSZJ=80000)
        COMMON/PYJETS/N,NPAD,K(KSZJ,5),P(KSZJ,5),V(KSZJ,5)
        common/sa26/ndiq(kszj),npt(kszj),ifcom(kszj),idi,idio
        common/saf/naf,nonaf,kaf(kszj,5),paf(kszj,5),vaf(kszj,5)
        COMMON/SA1_PY8/ i_mode, i_tune, KF_woDecay(100), &
               KF_proj, KF_targ, win, energy_B, psno, b_min, b_max


        ii = idio
        N0 = N

        do i = 1, N0, 1
            KS = K(i,1)
            KF = K(i,2)
            if( .NOT.IS_DIQUARK(KF) .OR. .NOT.IS_EXIST(KS,i_mode) ) cycle
!       Converts KF code of diquark to correspoding ones of quarks.
            KF1 =  MOD( KF/1000, 10 )
            KF2 =  MOD( KF/100,  10 )
!       Gives properties to the broken quarks.
            K( i, 2 )   = KF1
            K( N+1, 1 ) = K( i, 1 )
            K( N+1, 2 ) = KF2
            K( N+1, 3 ) = i + naf
            K( N+1, 4 ) = 0
            K( N+1, 5 ) = 0
            ii = ii + 1
            npt(ii) = N + 1 + naf
            ifcom(ii) = i + naf
!       Gives four momentum to the broken quarks.
            call bream( i, KF1, KF2 )
!       Gives four coordinate to the broken quarks.
            call coord(i)
            N = N + 1
        end do


        return
        end



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine bream( ii, kf1, kf2 )
!!      Gives four momenta to the broken quarks.
!       ii: line number of diquark in "PYJETS"
!       kf1, kf2: flavor codes of broken quarks
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        PARAMETER (KSZJ=80000)
        COMMON/PYJETS/N,NPAD,K(KSZJ,5),P(KSZJ,5),V(KSZJ,5)
        common/sa6_p/ithroq_p,ithrob_p,ich_p,non6_p,throe_p(4)
        dimension pi(4), ps(4), pp(20,5)
        integer decsuc


        am1 = amass(kf1)
        am2 = amass(kf2)
!       If zero mass diquark from PYTHIA.
        if( P(ii,5) <= 1D-15 )then
            am1 = 0D0
            am2 = 0D0
        end if
!       pp : four momenta & mass of broken quarks, local variable
        pp(1,5) = am1
        pp(2,5) = am2
!       ps : four momentum of diquark, local variable
        do i1=1,4,1
            ps(i1) = P(ii,i1)
        end do

!       Activates it for "decay method".
        goto 400

!       Broken quarks share out diquark three momentum randomly,
!        denoted as 'random three momentum method'
401     do i1=1,3,1
            pi(i1)   = PYR(1) * P(ii,i1)
            pp(2,i1) = ps(i1) - pi(i1)
            pp(1,i1) = pi(i1)
        end do
        pp11 = pp(1,1)
        pp12 = pp(1,2)
        pp13 = pp(1,3)
        pp14 = am1**2 + pp11**2 + pp12**2 + pp13**2
        if( pp14 <= 0D0 ) pp14 = 1D-20
        pp(1,4) = SQRT(pp14)
        pp21 = pp(2,1)
        pp22 = pp(2,2)
        pp23 = pp(2,3)
        pp24 = am2**2 + pp21**2 + pp22**2 + pp23**2
        if( pp24 <= 0D0 ) pp24 = 1D-20
        pp(2,4) = SQRT(pp24)

!       Activates it for 'random three momentum method'.
        goto 300

!       For 'decay method'.
400     continue
!       Decay method.
        decsuc = 1
        call decmom( ps, pp, am1, am2, decsuc )
        if( decsuc == 0 ) goto 401

300     continue
!       Adjusts four momentum conservation by iteration, no more than
!        4000 iterations
!       call conser(2,pp,ps)
        do i1=1,4,1
            p(ii,i1) = pp(1,i1)
        end do
        P(ii,5) = am1
        do i1=1,4,1
            P( N+1, i1 ) = pp(2,i1)
        enddo
        P( N+1, 5 ) = am2

!       Recalculates energies to ensure on-shell then collects the lost
!        4-momentum.
        P( ii,  4 ) = SQRT(P(ii,5)**2  + P(ii,1)**2  +P(ii,2)**2  + P(ii,3)**2)
        P( N+1, 4 ) = SQRT(P(N+1,5)**2 + P(N+1,1)**2 +P(N+1,2)**2 + P(N+1,3)**2)
        do i2=1,4,1
            throe_p(i2) = throe_p(i2) + ( ps(i2) - P(ii,i2) - P( N+1, i2 ) )
        end do


        return
        end



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine decmom( ps, pp, am1, am2, decsuc )
!!      Calculates four momentum of decayed particles.
!       ps: four momentum of decaying particle
!       pp: four momentum of decayed particles
!       am1 (am2): mass of decayed particle
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        PARAMETER (KSZJ=80000)
        COMMON/PYJETS/N,NPAD,K(KSZJ,5),P(KSZJ,5),V(KSZJ,5)
        dimension pi(4), pj(4), ps(4), pp(20,5), bb(3)
        integer decsuc


        decsuc = 1
!       Calculates the E and |p| of decayed particle in rest frame of
!        decaying particle.
        sm2 = ps(4)**2 - ps(1)**2 - ps(2)**2 - ps(3)**2
        if( sm2 < 1D-5 )then
            decsuc = 0
            return
        end if
        sm = dsqrt( sm2 )
!       pp(1,4) = ( sm2 - am2*am2 + am1*am1 ) / 2D0 / sm
!       pp(2,4) = ( sm2 - am1*am1 + am2*am2 ) / 2D0 / sm
        ppp = ( sm2 - (am1 + am2)**2 ) * ( sm2 - (am1 - am2)**2 )
        if( ppp <= 1D-5 )then
            decsuc = 0
            return
        end if
        ppp = SQRT( ppp ) / 2D0 / sm

!       Activate it for the exponential cos(seta) distribution.
!       goto 500

!       The direction of broken particle is sampled isotropically in "4pi".
        coset = 1D0 - 2D0 * PYR(1)
        if( ABS(coset) > 1D0 )then
            coset = coset / ABS( coset )
        end if
        siset = 1D0 - coset**2
        if( siset < 1D-28 ) siset = 1D-28
        siset = SQRT( siset )
100     cosi1 = PYR(1)
        cosi12 = cosi1**2
        eta2   = 2D0*PYR(1) - 1D0
        eta22  = eta2**2
        coseta = cosi12 + eta22
        if( coseta > 1D0 ) goto 100
        if( coseta < 1D-28 ) coseta = 1D-28
        cofi = ( cosi12 - eta22 ) / coseta
        sifi = 2D0 * cosi1 * eta2 / coseta
        goto 600

! 500   continue
!       cos(seta) is sampled from the exponential distribution when
!        0 < seta < pi/2 and its absolute value is assumed to be symmetry
!        about seta = pi/2. "fi" is assumed to be isotropic in "2pi".
        coset = LOG( 1D0 + 1.7183*PYR(1) )
        if( PYR(1) < 0.5D0 ) coset = -coset
        siset=1.-coset*coset
        if( siset < 1D-28 ) siset = 1D-28
        siset = SQRT( siset )
        fi = 2D0 * 3.1416 * PYR(1)
        cofi = COS( fi )
        sifi = SIN( fi )

600     continue
        pi(1) = ppp * siset * cofi
        pi(2) = ppp * siset * sifi
        pi(3) = ppp * coset
        pi4 = ppp**2 + am1**2
        if( pi4 < 1D-28 ) pi4 = 1D-28
        pi(4) = SQRT(pi4)
        pj(1) = -pi(1)
        pj(2) = -pi(2)
        pj(3) = -pi(3)
        pj4 = ppp**2 + am2**2
        if( pj4 < 1D-28 ) pj4 = 1D-28
        pj(4) = SQRT(pj4)
!       Rotates to the frame where the decaying particle, ps, is described.
!       Calculates the direction cosines of ps.
        fi1   = PYANGL( ps(1), ps(2) )
        ps12  = ps(1)**2 + ps(2)**2
        if( ps12 < 1D-28 ) ps12 = 1D-28
        ps12  = SQRT(ps12)
        cta1  = PYANGL( ps(3), ps12 )
        cfi1  = COS( fi1 )
        sfi1  = SIN( fi1 )
        ccta1 = COS( cta1 )
        scta1 = SIN( cta1 )
        sctas = siset
        cctas = coset
        sfis  = sifi
        cfis  = cofi
        pi(1) = cfi1 * ( ccta1*sctas*cfis + scta1*cctas ) - sfi1*sctas*sfis
        pi(2) = sfi1 * ( ccta1*sctas*cfis + scta1*cctas ) + cfi1*sctas*sfis
        pi(3) = ccta1 * cctas - scta1*sctas*cfis
        pi(1) = ppp * pi(1)
        pi(2) = ppp * pi(2)
        pi(3) = ppp * pi(3)
        do i=1,3,1
            pj(i) = - pi(i)
        end do
        pi4 = pi(1)**2 + pi(2)**2 + pi(3)**2 + am1**2
        if( pi4 < 1D-28 ) pi4 = 1D-28
        pj4 = pj(1)**2 + pj(2)**2 + pj(3)**2 + am2**2
        if( pj4 < 1D-28 ) pj4 = 1D-28
        pi(4) = SQRT(pi4)
        pj(4) = SQRT(pj4)
!       Boosts to the moving frame of decaying particle.
        ee = ps(4)
        if( ee < 1D-14 ) ee = 1D-14
        do i1=1,3,1
            bb(i1) = ps(i1) / ee
        end do
        call lorntz( 1, bb, pi, pj )
        pp(1,1) = pi(1)
        pp(1,2) = pi(2)
        pp(1,3) = pi(3)
        pp(1,4) = pi(4)
        pp(2,1) = pj(1)
        pp(2,2) = pj(2)
        pp(2,3) = pj(3)
        pp(2,4) = pj(4)


        return
        end



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine conser(np,pp,ps)
!!      Adjusts four momentum conservation.
!       np: the # of particles
!       pp: four momenta of particles have to be conserved
!       ps: the four momentum should be conserved to
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        common/papr/t0,cspipiKK,dep,ddt,edipi,epin,ecsnn,ekn,ecspsn,ecspsm, &
         rnt,rnp,rao,rou0,vneu,vneum,ecsspn,ecsspm,ecsen
        dimension pp(20,5),ps(4),ff(20),pxyz(3),arp(3)


        ps4 = ps(4)
        do i=1,3,1
            pxyz(i) = 0D0
        end do
        jj = 0
100     es = 0D0
        do i=1,np,1
            es = es + pp(i,4)
        end do
        fr = es / ps4

        if( ABS(1D0 - fr) <= dep) return

        do i=1,np,1
            amas  = pp(i,5)
            ppm   = pp(i,4) / amas
            ppf   = ppm / fr
            ff(i) = SQRT( ABS( ppf*ppf - 1D0 ) / ( ppm*ppm - 1D0 ) )
            do j=1,3,1
                ppp = ff(i) * pp(i,j)
                pp(i,j) = ppp
                pxyz(j) = pxyz(j) + ppp
            end do
        end do
        do i=1,3,1
            arp(i)  = ABS( 1D0 - pxyz(i) / ps(i) )
            pxyz(i) = pxyz(i) - ps(i)
        end do
        if( ABS(1D0 - fr) <= dep .AND. arp(1) <= dep .AND. arp(2) <= dep &
            .AND. arp(3) <= dep ) return
        do i=1,3,1
            pxyz(i) = pxyz(i) / np
        end do
        do i=1,np,1
            do j=1,3,1
                pp(i,j) = pp(i,j) - pxyz(j)
            end do
            pp5  = pp(i,5)
            pp52 = pp5*pp5
            pp(i,4) = SQRT( pp52 + pp(i,1)**2 + pp(i,2)**2 + pp(i,3)**2 )
        end do
        jj = jj + 1

        if( jj == 4000 )then
            write(9,*) 'infinitive loop may occur in subroutine ' // &
                       'conser(), which means four-momentum ' // &
                       'conservation needed is hard to be achieved, ' // &
                       'check value of PARAM(9)'
            return
        endif

        goto 100


        return
        end



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine coord(ii)
!!      Gives four positions to broken quarks
!       First broken quark takes the four position of diquark.
!       Second broken quark is arranged around first ones within
!        0.5 fm randomly in each of three position coordinates and has
!        same fourth position coordinate as diquark.
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        PARAMETER (KSZJ=80000)
        COMMON/PYJETS/N,NPAD,K(KSZJ,5),P(KSZJ,5),V(KSZJ,5)
        dimension rr(3)


        do i1=1,3,1
            rr(i1) = PYR(1) * 0.5D0
            V( N+1, i1 ) = V(ii,i1) + rr(i1)
            if( PYR(1) > 0.5D0 ) V( N+1, i1 ) = V(ii,i1) - rr(i1)
        end do
        V( N+1, 4 ) = V(ii,4)


        return
        end



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine copl(time)
!!      Calculates coordinate of center of mass of non-freeze-out system
!!      position of a particle, checking is it freezes out or not, is
!!       calculated with respect to this origin.
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        PARAMETER (KSZJ=80000)
        common/sa2/nsa,non2,ksa(kszj,5),psa(kszj,5),vsa(kszj,5)
        common/sa4/tau(kszj),tlco(kszj,4)
        common/wz/c17(500,3),ishp(kszj),tp(500),coor(3),p17(500,4)


        tt = time
        coor = 0D0
        samass = 0D0
        do ii=1,nsa,1
!           if( tau(ii) > tt ) cycle
            if( ishp(ii) == 0 ) cycle
            kf = ksa(ii,2)
            damass = PYMASS(KF)
            if( ABS(kf) == 213 .OR. kf == 113 ) amass1 = psa(ii,5)
            if( ABS(kf) == 213 .OR. kf == 113 &
                .AND. ABS(damass - amass1) > 0.001D0 ) damass = amass1
            samass = samass + damass
            do jj=1,3,1
                coor(jj) = coor(jj) + damass * vsa(ii,jj)
            end do
        end do
        coor = coor / MAX( 0.14D0, samass )


        return
        end



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine ctlcre
!!      Creates initial collision list.
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        PARAMETER (KSZJ=80000,NSIZE=750000)
        COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
        common/sa5/kfmax,kfaco(100),numb(100),numbs(100),non5, &
         disbe(100,100)
        common/sa35/ncpart,ncpar(kszj)
        common/ctllist/nctl,noinel(600),nctl0,nctlm
        common/ctllist_t/ lc(nsize,5), tc(nsize), tw(nsize)
        common/syspar/ipden,itden,suppm,suptm,suppc,suptc,r0p,r0t, &
         nap,nat,nzp,nzt,pio


!       time = 0D0
        nctl = 1
        dminf = 100D0
        ! In order to consider ppbar or pbarp.
        nzpab = ABS(nzp)
        nztab = ABS(nzt)
        nzpt  = nzpab + nztab
        napt  = nap   + nat
        ncpart = 0
        do i1=1,napt,1
            ncpar(i1) = 0
        end do
        ! projectile proton or lepton
        do 10 l=1,nzpab,1
            ! target proton
            do l1 = nzpab+1, nzpt, 1
                tc(nctl) = 0D0
                mtc = 0
                call coij( l, l1, nctl, mtc, dminf, iif, jf )
                if(mtc > 0)then
                    nctl = nctl + 1
                    mtc = 0
                end if
            end do
            ! target neutron
            do l1 = nap+nztab+1 , napt, 1
                tc(nctl) = 0D0
                mtc = 0
                call coij( l, l1, nctl, mtc, dminf, iif, jf )
                if(mtc > 0)then
                    nctl = nctl + 1
                    mtc = 0
                end if
            end do
10      continue
        ! projectile neutron
        do 20 l = nzpt+1, nap+nztab, 1
            ! target proton
            do l1 = nzpab+1, nzpt, 1
                tc(nctl) = 0D0
                mtc = 0
                call coij( l, l1, nctl, mtc, dminf, iif, jf )
                if(mtc > 0)then
                    nctl = nctl + 1
                    mtc = 0
                end if
            end do
            ! target neutron
            do l1= nap+nztab+1, napt, 1
                tc(nctl) = 0D0
                mtc = 0
                call coij( l, l1, nctl, mtc, dminf, iif, jf )
                if(mtc > 0)then
                    nctl = nctl + 1
                    mtc = 0
                end if
            end do
20      continue
        if(mtc == 0) nctl = nctl - 1
        if( nctl == 0 )then
!       at least one collision should occur. this collision has the smallest
!        'least approaching distance', that is guaranteed by the variable
!        'dminf'
!           lc(1,1) = iif
!           lc(1,2) = jf
!           tc(1) = 0.02D0
!           tw(1) = 0D0
!           nctl = 1
!       Regenerates a new event now.
            return
        end if
        do i1=1,napt,1
            if( ncpar(i1) >= 1 ) ncpart = ncpart + 1
        end do


        return
        end



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine coij(i,j,icp,mtc,dminf,iif,jf)
!!      Calculates collision time & fill up lc(J,1-2) as well as tc(J)
!!       for creating the particle initial collsion time list
!!       with particle linear trajectory assumption.
!       i (j): line # of colliding particle in 'sa2'
!       i: prijectile particle moving toward z direction
!       j: target particle moving toward -z direction
!       mtc=1: calculation is success; =0: calculation is fail
!       J: order # in collision time list
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        PARAMETER (KSZJ=80000)
        PARAMETER (NSIZE=750000)
        common/wz/c17(500,3),ishp(kszj),tp(500),coor(3),p17(500,4)
        common/papr/t0,cspipiKK,dep,ddt,edipi,epin,ecsnn,ekn,ecspsn,ecspsm, &
         rnt,rnp,rao,rou0,vneu,vneum,ecsspn,ecsspm,ecsen   ! 060813
        common/ctllist_t/ lc(nsize,5), tc(nsize), tw(nsize)
        common/syspar/ipden,itden,suppm,suptm,suppc,suptc,r0p,r0t, &
         nap,nat,nzp,nzt,pio
        common/sa2/nsa,non2,ksa(kszj,5),psa(kszj,5),vsa(kszj,5)
        common/sa35/ncpart,ncpar(kszj)   ! 280722
        dimension dr(3),db(3),dv(3),px(4),py(4),vi(4),vj(4), &
         pi(4),pj(4),b(3)


        mtc = 0
        ki=ksa(i,2)
        kj=ksa(j,2)
        pi(4)=psa(i,4)   ! energy of particle i
        if(pi(4) < 1.e-20) pi(4)=1.e-20
        pj(4)=psa(j,4)   ! energy of particle j
        if(pj(4) < 1.e-20) pj(4)=1.e-20
        deno6=pi(4)+pj(4)   ! total energy of ij-th colliding pair
        do k=1,3
            pi(k)=psa(i,k)   ! three momentum of particle i
            pj(k)=psa(j,k)   ! three momentum of particle j
            b(k)=(pi(k)+pj(k))/deno6   ! velocity of i and j cms
        enddo
        ilo=0
        call lorntz(ilo,b,pi,pj)   ! boost to cms of i and j for momentum
        do l=1,3
            px(l)=vsa(i,l)    ! three coordinate of particle i
            py(l)=vsa(j,l)    ! three coordinate of particle j
        enddo
        px(4)=0.   ! time of particle i
        py(4)=0.   ! time of particle j
        call lorntz(ilo,b,px,py)   ! boost to cms of i and j for coordinate
!       calculate instant relative distance between i and j
        rb=0.
        bb=0.
        rr=0.
        rtai=0.
        do k=1,3
            vi(k)=pi(k)/pi(4)   ! three velocity of i
            vj(k)=pj(k)/pj(4)   ! three velocity of j
        enddo
        do k=1,3
            dr(k)=px(k)-py(k)-(vi(k)*px(4)-vj(k)*py(4))
!       instant relative distance between i and j in k dimension (in above
!        equation i is target particle and j is projetile particle indeed,
!        but for calculation of relative distance, it does matter whether i
!        is target j is projectile or opposite
            db(k)=vi(k)-vj(k)   ! relative velocity between i and j
            dv(k)=db(k)
            rb=rb+dr(k)*db(k)
!       rb=(relative distance)*(relative velocity)=(relative distance)^2/t
            bb=db(k)**2+bb
!       bb: square of relative velocity
            rr=rr+dr(k)*dr(k)   ! (?)
!       dr(k)*dr(k): (relative distance between i and j)^2 in k dimension (?)
        enddo
        if(bb <= 1.e-10) return   ! (bb <= 1.e-10) it is true
        tcol=0.-rb/bb
!       rb/bb=(relative distance)^2/t/(relative velocity)^2
!        =(relative distance)^2/t/((relative distance)^2/t^2)=t
!       why tcol=-rb/bb ? it is because of i must be projectile particle
!        and j must be target particle
!        if(tcol-px(4)  <=  ddt) return
!        if(tcol-py(4)  <=  ddt) return
!       for collision to occur,time must one step ahead
!sa     if(tcol < 1.0e-7)return
        do iik=1,3
            dr(iik)=px(iik)-py(iik)-(vi(iik)* &
             px(4)-vj(iik)*py(4))+tcol*db(iik)
!sa     relative distance between i and j at current time
            rtai=rtai+dr(iik)*dr(iik)
        enddo
        sg=rtai   ! square of instant relative distance between i and j
!       sg=rr+tcol*rb
!       if(sg < 0)then
!           write(*,*)'sg=',sg   !
!           return
!       endif
        dmin=dsqrt(sg)   ! instant relative distance between i and j
        if(dmin < dminf)then
            dminf=dmin
            iif=i
            jf=j
        endif
        if(ipden < 2 .and. dmin > ecsnn) return
!       ecsnn: largest interaction distance of NN
        if(ipden > 2 .and. dmin > ecsen) return
!       ecsen: largest interaction distance of eN
!       distance between the two particles should be smaller than ecsnn (ecsen)
        do ik=1,3
            px(ik)=px(ik)+vi(ik)*(tcol-px(4))
            py(ik)=py(ik)+vj(ik)*(tcol-py(4))
        enddo
!       move along Newton trajectory in CMS
        px(4)=tcol
        py(4)=tcol
        ilo=1

        call lorntz(ilo,b,px,py)

!       transform back to Lab.
        if(px(4) > py(4)) px(4)=py(4)
        tcol=px(4)
        drmax=rao*dmax1(rnt,rnp)
        if(tcol <= drmax) goto 180
        return
180     tc(icp)=tcol
        mtc=1
        lc(icp,1)=i
        lc(icp,2)=j
        tw(icp) = 0D0
        ncpar(i) = ncpar(i) + 1
        ncpar(j) = ncpar(j) + 1


        return
        end



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine find(icp,tcp,ico)
!!      Finds out the binary collision with minimum collision time.
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        PARAMETER (NSIZE=750000)
        common/ctllist/nctl,noinel(600),nctl0,nctlm
        common/ctllist_t/ lc(nsize,5), tc(nsize), tw(nsize)


        icp = 0
        tcp = 20000D0
        do i=1,nctl,1
            if( ico == 0 ) goto 100
            if( tc(i) <= 1D-7 ) goto 241
100         if( tcp < tc(i) )  goto 241
            icp = i
            tcp = tc(i)
241         continue
        end do
        if( nctl == 0 ) icp=0


        return
        end



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine updtlp(time,i_mode)
!!      Updates collision list after calling 'PYTHIA' successfully.
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        PARAMETER (KSZJ=80000)
        PARAMETER (NSIZE=750000)
        common/sa1/kjp21,non1,bp,iii,neve,nout,nosc
        common/sa2/nsa,nonsa,ksa(kszj,5),psa(kszj,5),vsa(kszj,5)
        common/sa5/kfmax,kfaco(100),numb(100),numbs(100),non5, &
         disbe(100,100)
        common/sa12/ppsa(5),nchan,nsjp,sjp,taup,taujp
        common/sa14/ipyth(2000),idec(2000),iwide(2000)
        common/sa30/vneump,vneumt,mstptj
        common/sbh/nbh,non,kbh(kszj,5),pbh(kszj,5),vbh(kszj,5)
        common/ctllist/nctl,noinel(600),nctl0,nctlm
        common/ctllist_t/ lc(nsize,5), tc(nsize), tw(nsize)
        common/papr/t0,cspipiKK,dep,ddt,edipi,epin,ecsnn,ekn,ecspsn,ecspsm, &
         rnt,rnp,rao,rou0,vneu,vneum,ecsspn,ecsspm,ecsen
        common/syspar/ipden,itden,suppm,suptm,suppc,suptc,r0p,r0t, &
         nap,nat,nzp,nzt,pio
        common/coal1/bmrat,i_mm
!       ipyth: store line number of produced hardon in hadron list 'sa2'


!       loop over old colliding pairs
        j=0
        do i=1,nctl,1
            if( ABS(tc(i)-time) <= ddt ) cycle
!       through away the pair whih tc<= time
            j=j+1
            tc(j)=tc(i)
            tw(j)=tw(i)
            do m=1,5
                lc(j,m)=lc(i,m)
            end do
        end do

        nctl = j + 1

!       for B-framework (PYTHIA-like framework)
        if( mstptj == 1)then
            nctl = j
            return
        end if

        m2=numb(2)
        m4=numb(4)
        m6=numb(6)
        m7=numb(7)
!       m9=numb(9)
!       m17=numb(17)
!       m19=numb(19)
!       m25=numb(25)
!       m29=numb(29)
!       m32=numb(32)
!       m34=numb(34)
!       note: # of produced hadrons equal to zero (nbh=0) after call 'PYTHIA'
!        in case of w/o reconstruction leading proton
!       proceed for case of with reconstruction leading nucleon
!       constract hadron collision pair composed of one from produced hadrons
!        and another one in 'sa2'
!       loop over produced hadrons in 'sbh'
        do j11=1,nbh,1
            j1 = ipyth(j11)
            kfjab = ABS( ksa(j1,2) )
            if( ( i_mm  ==  2 .and. &
             ( kfjab /= 2212 .and. kfjab /= 2112 .and. &
               kfjab /= 11   .and. kfjab /= 12   .and. kfjab /= 13 .and. &
               kfjab /= 14   .and. kfjab /= 15   .and. kfjab /= 16 ) ) &
             .OR. &
             ( i_mm == 6.and. &
             ( kfjab /= 2212 .and. kfjab /= 2112 .and. kfjab /= 211.and. &
               kfjab /= 11   .and. kfjab /= 12   .and. kfjab /= 13 .and. &
               kfjab /= 14   .and. kfjab /= 15   .and. kfjab /= 16 ) ) &
            ) cycle
!       consider only the reinteraction among nucleons & nucleon with lepton
!       nucleon with pion, and among pions.

!       loop over particle list ('sa2')
            mm = numb(i_mm)
!       consider only the reinteraction of j11 with nucleons
            do i=1,mm,1
                do j22=1,nbh,1
                    j2 = ipyth(j22)
                    ! avoid particle collide with itself
                    if( i == j2 ) goto 600
                end do
                i1 = i
                kfi = ksa(i1,2)
                iflag = 0
                call rsfilt(j1,i1,iflag,i_mode)
                if(iflag == 0) cycle
                tc(nctl) = 0D0
                call tcolij(j1,i1,time,nctl,i_mode)
                if( tc(nctl) > 1D-7 ) nctl = nctl + 1
600         end do
        end do
        nctl = nctl - 1


        return
        end



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine rsfilt(l,l1,iflag,i_mode)
!!      Subroutine rsfilt plays the role of first range filter.
!!      Subroutine intdis plays the role of second range filter.
!!      Collision pairs not interested can not filter through both of
!!       rsfilt and intdis.
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        PARAMETER (KSZJ=80000)
        COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
        common/sa2/nsa,non2,ksa(kszj,5),psa(kszj,5),vsa(kszj,5)
        common/sa5/kfmax,kfaco(100),numb(100),numbs(100),non5, &
         disbe(100,100)
        common/wz/c17(500,3),ishp(kszj),tp(500),coor(3),p17(500,4)


        iflag = 0
        kl    = ksa(l,2)
        kl1   = ksa(l1,2)
        klab  = iabs(kl)
        kl1ab = iabs(kl1)
        if( l == l1 ) return
        if( ishp(l) == 0 .or. ishp(l1) == 0 ) return

!       high energy B-, C- and D-frameworks
        if( i_mode /= 1 )then
            ! consider NN collision
            if( klab == 2212 .and. (kl1ab == 2112 .or.kl1ab == 2212) ) iflag = 1
            if( klab == 2112 .and. (kl1ab == 2112 .or.kl1ab == 2212) ) iflag = 1
            ! consider lN collision
            if( (klab >= 11.and.klab <= 16).and.(kl1 == 2112.or.kl1 == 2212) ) &
                iflag=1
            if( (klab == 2112.or.klab == 2212) .and. (kl1ab >= 11.and. &
                kl1ab <= 16) ) iflag=1
        end if

!       low energy A-framework
        if( i_mode == 1 )then
        ! consider NN, (Delta)N and (pi)N collisions
            if( ( klab==2212 .or. klab==2112 .or. klab==211 .or. kl==1114 &
                .or.kl == 2114.or.kl == 2214.or.kl == 2224 ) .and. &
                (kl1ab == 2212.or.kl1ab == 2112.or.kl1ab == 211.or.kl1 == 1114 &
                .or.kl1 == 2114.or.kl1 == 2214.or.kl1 == 2224) ) iflag=1
        end if


        return
        end



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine tcolij(l,l1,time,icp,i_mode)
!!      Calculates collision time & fill up lc(i,1-2) as well as tc(i).
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        PARAMETER (KSZJ=80000)
        PARAMETER (NSIZE=750000)
        COMMON/PYPARS/MSTP(200),PARP(200),MSTI(200),PARI(200)
        common/sa1/kjp21,non1,bp,iii,neve,nout,nosc
        common/sa10/csnn,cspin,cskn,cspipi,cspsn,cspsm,rcsit,ifram, &
         iabsb,iabsm,i_sigma_AQM,ajpsi,csspn,csspm,csen
        common/sa4/tau(kszj),tlco(kszj,4)
        common/sa2/nsa,non2,ksa(kszj,5),psa(kszj,5),vsa(kszj,5)
        common/papr/t0,cspipiKK,dep,ddt,edipi,epin,ecsnn,ekn,ecspsn,ecspsm, &
         rnt,rnp,rao,rou0,vneu,vneum,ecsspn,ecsspm,ecsen   ! 060813
        common/wz/c17(500,3),ishp(kszj),tp(500),coor(3),p17(500,4)
        common/ctllist_t/ lc(nsize,5), tc(nsize), tw(nsize)
        dimension dr(3),db(3),pi(4),pj(4),vi(3),vj(3)
        dimension ri(4),rj(4),rfi(4),rfj(4),b(3)


        pel  = psa(l,4)
        pel1 = psa(l1,4)
        if( pel < 1D-20 ) pel   = 1D-20
        if( pel1 < 1D-20 ) pel1 = 1D-20
        pi(4) = pel
        pj(4) = pel1
        do i=1,3
            pi(i) = psa(l,i)
            pj(i) = psa(l1,i)
            b(i)  = ( pi(i) + pj(i) ) / ( pi(4) + pj(4) )
        enddo
        ilo=0
        call lorntz(ilo,b,pi,pj)
!       perform Lorentz transf. to CMS frame for momentum.
        bta = dsqrt( b(1)**2 + b(2)**2 + b(3)**2 )
!       if boost is too violent,put particles on mass shell by hand.
        if( bta > 0.99999D0 )then
            kl    = ksa(l,2)
            kl1   = ksa(l1,2)
            klab  = iabs(kl)
            kl1ab = iabs(kl1)
            bmi   = PYMASS(kl)
            bmj   = PYMASS(kl1)
            pi(4) = dsqrt( bmi**2 + pi(1)**2 + pi(2)**2 + pi(3)**2 )
            pj(4) = dsqrt( bmj**2 + pj(1)**2 + pj(2)**2 + pj(3)**2 )
        end if
        ss=pi(4)+pj(4)
!       pair do not enter the collision list if the threshold is too small.
!       if( ( (klab == 2211 .or. klab == 2112) .and. &
!            (kl1ab == 2112 .or. kl1ab == 2212) ) .and. ss <= parp(2) ) goto 10

        do i=1,4
            ri(i) = vsa(l,i)
            rj(i) = vsa(l1,i)
        enddo
!       ri(4) = time
!       rj(4) = time
        call lorntz(ilo,b,ri,rj)
!       perform Lorentz transf. to CMS frame for coordinate.
        rb = 0.
        bb = 0.
        rr = 0.
        rtai = 0.
        kflag = 0
        do ik=1,3
            vi(ik) = pi(ik) / pi(4)
            vj(ik) = pj(ik) / pj(4)
        enddo

        do i=1,3
            rfi(i) = vsa(l,i)  + ( tau(l)  - time) * ( psa(l,i)  / psa(l,4) )
            rfj(i) = vsa(l1,i) + ( tau(l1) - time) * ( psa(l1,i) / psa(l1,4) )
        enddo
        rfi(4) = tau(l)
        rfj(4) = tau(l1)
        call lorntz(ilo,b,rfi,rfj)
!       gamli = psa(l,4)  / psa(l,5)
!       gamlj = psa(l1,4) / psa(l1,5)
        ctaui = rfi(4)
        ctauj = rfj(4)
        tcol  = ctaui
        if(ctaui < ctauj) tcol = ctauj
        do ik=1,3
            db(ik) = ( vi(ik) - vj(ik) ) * tcol
            dr(ik) = ri(ik) - rj(ik) - ( vi(ik)*ri(4) - vj(ik)*rj(4) ) + db(ik)
            rtai   = rtai + dr(ik)*dr(ik)
        enddo
        dott = 0.
        do ik=1,3
            dott = dr(ik)*pi(ik) + dott
        enddo
!       dott = -1
        if(dott >= 0.)then
            kflag = 1
            if(tcol <= ri(4)) goto 10
            if(tcol <= rj(4)) goto 10
        else
            rtai = 0.
            do ik=1,3
                dr(ik) = ri(ik) - rj(ik) - ( vi(ik)*ri(4) - vj(ik)*rj(4) )
                db(ik) = vi(ik) - vj(ik)
                rb = rb + dr(ik)*db(ik)
                bb = bb + db(ik)*db(ik)
                rr = rr + dr(ik)*dr(ik)
            enddo
            if(bb <= 1.e-10) goto 10
            tcol = 0. - rb/bb
            if(tcol <= ri(4)) goto 10
            if(tcol <= rj(4)) goto 10
            if(tcol-ctaui <= 0.) goto 10
            if(tcol-ctauj <= 0.) goto 10
!       for collision to occur,time must one step ahead
!Tai
            do ik=1,3
                dr(ik) = ri(ik)-rj(ik) - (vi(ik)*ri(4)-vj(ik)*rj(4))+tcol*db(ik)
                rtai = rtai + dr(ik)*dr(ik)
            enddo
!           gamai = pi(4) / PYMASS( ksa(l,2)  )
!           gamaj = pj(4) / PYMASS( ksa(l1,2) )

! TAIAN

!       when collision happens,particles should already be produced
!       we give a zero formation time for particles produced after
!       calling 'PYTHIA'
            sg1 = rr + tcol*rb
!tai
        endif
        sg = rtai
!       if(sg1 < 0.and.iii <= 50)then
!           write(*,*)'sar,tair=',sg1,sg
!           dmin = 0.
!           tcol = -rr/rb
!           goto 20
!       endif

        dmin = dsqrt(sg)
!       calculate the interaction distance between particles l & l1.
! 20    call intdis(l,l1,rsig,i_mode)
        call intdis(l,l1,rsig,i_mode)
!       distance between the two particles should be smaller than rsig
        if(dmin > rsig) goto 10
!       move along Newton trajectory in CMS
        do ik=1,3
            ri(ik) = ri(ik) + vi(ik)*( tcol - ri(4) )
            rj(ik) = rj(ik) + vj(ik)*( tcol - rj(4) )
        enddo
        ri(4) = tcol
        rj(4) = tcol
        ilo=1

        call lorntz(ilo,b,ri,rj)

!       transform back to Lab.
        tcol1 = ri(4)
        tcol2 = rj(4)
        if(kflag == 0)then
            if( tcol1 - tau(l)  < 0. ) goto 10
            if( tcol2 - tau(l1) < 0. ) goto 10
        else
            if( tcol1 - tau(l)  < -1.E-4 ) goto 10
            if( tcol2 - tau(l1) < -1.E-4 ) goto 10
        endif
        if( ri(4) > rj(4) ) ri(4) = rj(4)
        tcol = ri(4)
        if(tcol <= time) goto 10
!       collision happens in the future
!       if(ifram == 0) coor(3) = coor(3) + rnt
        do i=1,3
            ri(i) = vsa(l,i)  + psa(l,i)*(tcol-time)  / pel  - coor(i)
            rj(i) = vsa(l1,i) + psa(l1,i)*(tcol-time) / pel1 - coor(i)
        enddo
        rri = dsqrt( ri(1)*ri(1) + ri(2)*ri(2) + ri(3)*ri(3) )
        rrj = dsqrt( rj(1)*rj(1) + rj(2)*rj(2) + rj(3)*rj(3) )
!       the rnt in rao*max(rnt,rnp)+rnt is due to the fact that
!        we could not know the postion of the mass-center in the future
        rrr = rao*dmax1(rnt,rnp)
        if(ifram == 0) rrr = rao*dmax1(rnt,rnp) + rnt
!       if( dabs(ksa(l,2)) > 1000 .or. dabs(ksa(l1,2)) > 1000 ) rrr = 1.E+10*rrr
        if(rri > rrr) goto 10
        if(rrj > rrr) goto 10
!       particles under consideration must be still within considered region
!        when the collision happens
        if(tcol <= rrr) goto 18

        return

18      tc(icp)   = tcol
        lc(icp,1) = l
        lc(icp,2) = l1


10      return
        end



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine intdis(l,l1,rsig,i_mode)
!!      Calculates interaction distance,rsig, between colliding pair of l and l1
!!      It plays also the role of second range filter
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        PARAMETER (KSZJ=80000)
        common/sa2/nsa,non2,ksa(kszj,5),psa(kszj,5),vsa(kszj,5)
        common/papr/t0,cspipiKK,dep,ddt,edipi,epin,ecsnn,ekn,ecspsn,ecspsm, &
         rnt,rnp,rao,rou0,vneu,vneum,ecsspn,ecsspm,ecsen
        common/sa13/kjp20,non13,vjp20,vjp21,vjp22,vjp23


        rsig  = 0D0
        kl    = ksa(l,2)
        kl1   = ksa(l1,2)
        klab  = iabs(kl)
        kl1ab = iabs(kl1)
        idpl  = 0
        idpl1 = 0

        if( i_mode /= 1 )then
            if( klab  == 2212 .or. klab  == 2112 ) idpl  = 1
            if( kl1ab == 2212 .or. kl1ab == 2112 ) idpl1 = 1
            if( klab  == 211 ) idpl  = 3
            if( kl1ab == 211 ) idpl1 = 3
            if( klab  >= 11 .and. klab  <= 16 ) idpl  = 8
            if( kl1ab >= 11 .and. kl1ab <= 16 ) idpl1 = 8
            if( idpl == 1 .and. idpl1 == 1) rsig = ecsnn
            if( idpl == 8 .and. idpl1 == 1) rsig = ecsen
            if( idpl == 1 .and. idpl1 == 8) rsig = ecsen
            if(idpl == 3 .and. idpl1 == 1) rsig = epin
            if(idpl == 1 .and. idpl1 == 3) rsig = epin
            if(idpl == 3 .and. idpl1 == 3) rsig = edipi
        end if

        if(i_mode == 1)then
            if(klab  == 2212 .or. klab  == 2112)  idpl = 1   ! N
            if(kl1ab == 2212 .or. kl1ab == 2112) idpl1 = 1   ! N
            if(klab  == 1114 .or. klab  == 2114 .or. klab  == 2214 .or. klab &
             == 2224 .or. klab == 2224) idpl = 2   ! Delta
            if(kl1ab == 1114 .or. kl1ab == 2114 .or. kl1ab == 2214 .or. kl1ab &
             == 2224.or.kl1ab == 2224) idpl1 = 2   ! Delta
            if(klab  == 211)  idpl = 3   ! pion
            if(kl1ab == 211) idpl1 = 3   ! pion
            if(idpl == 1 .and. idpl1 == 1) rsig = ecsnn   ! NN
            if(idpl == 2 .and. idpl1 == 1) rsig = ecsnn   ! (Delta)N
            if(idpl == 1 .and. idpl1 == 2) rsig = ecsnn   ! N(Delta)
            if(idpl == 3 .and. idpl1 == 1) rsig = epin    ! (pi)N
            if(idpl == 1 .and. idpl1 == 3) rsig = epin    ! N(pi)
        endif


        return
        end



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine updatl(ic,jc,time,i_mode)
!!      Updates collision time list after elastic scattering.
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        PARAMETER (KSZJ=80000)
        PARAMETER (NSIZE=750000)
        COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
        common/sa2/nsa,non2,ksa(kszj,5),psa(kszj,5),vsa(kszj,5)
        common/sa5/kfmax,kfaco(100),numb(100),numbs(100),non5, &
         disbe(100,100)
        common/sa12/ppsa(5),nchan,nsjp,sjp,taup,taujp
        common/ctllist/nctl,noinel(600),nctl0,nctlm
        common/ctllist_t/ lc(nsize,5), tc(nsize), tw(nsize)
        common/papr/t0,cspipiKK,dep,ddt,edipi,epin,ecsnn,ekn,ecspsn,ecspsm, &
         rnt,rnp,rao,rou0,vneu,vneum,ecsspn,ecsspm,ecsen
        common/syspar/ipden,itden,suppm,suptm,suppc,suptc,r0p,r0t, &
         nap,nat,nzp,nzt,pio
        common/coal1/bmrat,i_mm


!       loop over old colliding pairs
        j=0
        do i=1,nctl,1
            i1 = lc(i,1)
            j1 = lc(i,2)
!           ia=(i1-ic)*(j1-jc)*(i1-jc)*(j1-ic)
!           if(ia == 0) cycle
            if(i1 == ic .or. i1 == jc) cycle
            if(j1 == ic .or. j1 == jc) cycle
            if(i1 == j1) cycle
            if( ABS(tc(i)-time) <= ddt ) cycle
!       through away the pair whih tc<= time
            j = j + 1
            tc(j) = tc(i)
            tw(j) = tw(i)
            do m=1,5,1
                lc(j,m) = lc(i,m)
            end do
        end do

        nctl = j + 1
!       loop over particle list

        j1 = ic
        do ik=1,2,1
            kfjab = ABS( ksa(j1,2) )

            if( ( i_mm  ==  2 .and. &
             ( kfjab /= 2212 .and. kfjab /= 2112 .and. &
               kfjab /= 11   .and. kfjab /= 12   .and. kfjab /= 13 .and. &
               kfjab /= 14   .and. kfjab /= 15   .and. kfjab /= 16) ) &
             .OR. &
             ( i_mm == 6.and. &
             ( kfjab /= 2212 .and. kfjab /= 2112 .and. kfjab /= 211.and. &
               kfjab /= 11   .and. kfjab /= 12   .and. kfjab /= 13 .and. &
               kfjab /= 14   .and. kfjab /= 15   .and. kfjab /= 16) ) &
            )then
                j1 = jc
                cycle
            end if
!       consider only the reinteraction among nucleons & nucleon with lepton
!       nucleon with pion, and among pions.
            mm = numb(i_mm)
            do i=1,mm,1
!       forbiden scattered particles colliding with each other
                if(j1 == ic .and. i == jc) cycle
                if(j1 == jc .and. i == ic) cycle
                if( i == j1 ) cycle
                i1 = i
                iflag = 0
                call rsfilt(j1,i1,iflag,i_mode)
                if(iflag == 0) cycle
                tc(nctl) = 0D0
                call tcolij(i1,j1,time,nctl,i_mode)
                if(tc(nctl) > 1D-7) nctl = nctl + 1
            end do
            j1 = jc
        end do
        nctl = nctl - 1


        return
        end



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine pauli(ii,ppaul)
!!      Calculates the unoccupation probability (ppaul) of particle ii
!!       in 'PYJETS'.
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        PARAMETER (KSZJ=80000)
        COMMON/PYJETS/N,NPAD,K(KSZJ,5),P(KSZJ,5),V(KSZJ,5)
        common/saf/nsa,nonsa,ksa(kszj,5),psa(kszj,5),vsa(kszj,5)
        dimension rkk(1000,4),pkk(1000,4),rr(4),pp(4),b(3)


        kf=k(ii,2)   ! new
        xxii=v(ii,1)
        yyii=v(ii,2)
        zzii=v(ii,3)
        ttii=v(ii,4)
        pxii=p(ii,1)
        pyii=p(ii,2)
        pzii=p(ii,3)
        eeii=p(ii,4)
        b(1)=pxii/eeii
        b(2)=pyii/eeii
        b(3)=pzii/eeii
!       pick up the partons with same flavor as ii from 'PYJETS' and
!        'saf'
        nkk=0
!       the new produced partons, except ii, in 'PYJETS' are also to be
!        considered
        do i=1,n   ! loop over new
            if(i == ii) goto 100
            kfi1=k(i,2)
            if(kfi1 == kf)then
                nkk=nkk+1
                do j=1,4
                    rkk(nkk,j)=v(i,j)
                    pkk(nkk,j)=p(i,j)
                enddo
            endif
100     enddo
        do i1=1,nsa   ! loop over old
            kfi1=ksa(i1,2)
            if(kfi1 == kf)then
                nkk=nkk+1
                do j=1,4
                    rkk(nkk,j)=vsa(i1,j)
                    pkk(nkk,j)=psa(i1,j)
                enddo
            endif
        enddo
!       boost to the rest frame of ii
        ilo=0
        do 200 j2=1,nkk
            do j1=1,4
                rr(j1)=rkk(j2,j1)
                pp(j1)=pkk(j2,j1)
            enddo
            call lorntz(ilo,b,rr,pp)
            do j1=1,4
                rkk(j2,j1)=rr(j1)
                pkk(j2,j1)=pp(j1)
            enddo
200     enddo
        rr(1)=xxii
        rr(2)=yyii
        rr(3)=zzii
        rr(4)=ttii
        call lorntz(ilo,b,rr,rr)
        xxii=rr(1)
        yyii=rr(2)
        zzii=rr(3)
        ttii=rr(4)
!       calculate the number of partons occupied in or on the surface of
!        six dimension cub (around ii): (dr*dp)**3=h**3, dr*dp=h, h=1.24
!        GeV*fm/c, if dr=1.0 fm then dp=1.24 GeV/c
        anq=0   ! statistics of ocuupation number in (dr*dp)**3=h**3
        do i1=1,nkk
            dxx=xxii-rkk(i1,1)
            dyy=yyii-rkk(i1,2)
            dzz=zzii-rkk(i1,3)
!       following three staments for without boost
!           dpx=pxii-pkk(i1,1)
!           dpy=pyii-pkk(i1,2)
!           dpz=pzii-pkk(i1,3)
!       following three staments for with boost
            dpx=pkk(i1,1)
            dpy=pkk(i1,2)
            dpz=pkk(i1,3)
            dxx=dabs(dxx)
            dyy=dabs(dyy)
            dzz=dabs(dzz)
            dpx=dabs(dpx)
            dpy=dabs(dpy)
            dpz=dabs(dpz)
            if(dxx <= 0.5.and.dyy <= 0.5.and.dzz <= 0.5.and. &
             dpx <= 0.62.and.dpy <= 0.62.and.dpz <= 0.62) anq=anq+1.
        enddo
        proba=anq/6.
!       6=2*3, spin and colour degeneracies of quark (antiquark)
        ppaul=1.-proba


        return
        end



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine crosep(rots,csen)
!!exponential interpolation for ep total cross section (in mbarn)
!! calculated with herafitter by Xing-Long Li on 10/Dec./2013
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        IMPLICIT INTEGER(I-N)
        !!the lower limit of rots (GeV)
        parameter (rotsMin=2.d0)
        !!the upper limit of rots (GeV)
        parameter (rotsMax=1003)
        parameter (Ndat=315)
!!sect: ep cross section data for sqrtS range: 2~1003GeV
!!the unit is mbarn
!!when i=1,sqrtS=2.0 GeV; sqrtS=2.0*(1.02)^(i-1) GeV otherwise
        dimension sect(Ndat)
        parameter (sect=[ &
        0.125377D-04,0.136063D-04,0.147199D-04,0.158792D-04,0.170824D-04, &
        0.183314D-04,0.196228D-04,0.209579D-04,0.223386D-04,0.237574D-04, &
        0.252260D-04,0.267283D-04,0.282729D-04,0.298683D-04,0.314933D-04, &
        0.331676D-04,0.348724D-04,0.366252D-04,0.384099D-04,0.402340D-04, &
        0.421004D-04,0.439997D-04,0.459366D-04,0.479080D-04,0.499176D-04, &
        0.519633D-04,0.540432D-04,0.561548D-04,0.583026D-04,0.604863D-04, &
        0.626984D-04,0.649432D-04,0.672178D-04,0.695343D-04,0.718663D-04, &
        0.742366D-04,0.766461D-04,0.790697D-04,0.815331D-04,0.840216D-04, &
        0.865379D-04,0.890822D-04,0.916557D-04,0.942586D-04,0.968830D-04, &
        0.995384D-04,0.102212D-03,0.104924D-03,0.107655D-03,0.110411D-03, &
        0.113187D-03,0.115993D-03,0.118824D-03,0.121673D-03,0.124552D-03, &
        0.127461D-03,0.130369D-03,0.133312D-03,0.136294D-03,0.139258D-03, &
        0.142289D-03,0.145318D-03,0.148358D-03,0.151438D-03,0.154537D-03, &
        0.157648D-03,0.160778D-03,0.163928D-03,0.167102D-03,0.170296D-03, &
        0.173515D-03,0.176740D-03,0.179988D-03,0.183258D-03,0.186551D-03, &
        0.189855D-03,0.193175D-03,0.196515D-03,0.199880D-03,0.203265D-03, &
        0.206648D-03,0.210060D-03,0.213497D-03,0.216952D-03,0.220408D-03, &
        0.223888D-03,0.227392D-03,0.230909D-03,0.234442D-03,0.237992D-03, &
        0.241560D-03,0.245143D-03,0.248743D-03,0.252364D-03,0.256010D-03, &
        0.259659D-03,0.263324D-03,0.267009D-03,0.270715D-03,0.274444D-03, &
        0.278182D-03,0.281933D-03,0.285693D-03,0.289480D-03,0.293325D-03, &
        0.297112D-03,0.300951D-03,0.304824D-03,0.308704D-03,0.312598D-03, &
        0.316505D-03,0.320433D-03,0.324381D-03,0.328354D-03,0.332347D-03, &
        0.336349D-03,0.340362D-03,0.344403D-03,0.348459D-03,0.352536D-03, &
        0.356641D-03,0.360763D-03,0.364893D-03,0.369052D-03,0.373232D-03, &
        0.377425D-03,0.381629D-03,0.385859D-03,0.390121D-03,0.394418D-03, &
        0.398725D-03,0.403041D-03,0.407380D-03,0.411739D-03,0.416122D-03, &
        0.420532D-03,0.424966D-03,0.429429D-03,0.433908D-03,0.438403D-03, &
        0.442923D-03,0.447466D-03,0.452044D-03,0.456646D-03,0.461274D-03, &
        0.465905D-03,0.470575D-03,0.475274D-03,0.480000D-03,0.484746D-03, &
        0.489511D-03,0.494317D-03,0.499153D-03,0.504004D-03,0.508874D-03, &
        0.513782D-03,0.518727D-03,0.523706D-03,0.528696D-03,0.533711D-03, &
        0.538757D-03,0.543840D-03,0.548953D-03,0.554107D-03,0.559291D-03, &
        0.564487D-03,0.569720D-03,0.574975D-03,0.580275D-03,0.585612D-03, &
        0.590979D-03,0.596380D-03,0.601802D-03,0.607256D-03,0.612757D-03, &
        0.618283D-03,0.623860D-03,0.629465D-03,0.635094D-03,0.640766D-03, &
        0.646474D-03,0.652220D-03,0.658001D-03,0.663819D-03,0.669669D-03, &
        0.675552D-03,0.681488D-03,0.687464D-03,0.693475D-03,0.699515D-03, &
        0.705602D-03,0.711726D-03,0.717890D-03,0.724097D-03,0.730344D-03, &
        0.736634D-03,0.742995D-03,0.749397D-03,0.755772D-03,0.762216D-03, &
        0.768727D-03,0.775280D-03,0.781876D-03,0.788516D-03,0.795202D-03, &
        0.801928D-03,0.808704D-03,0.815513D-03,0.822366D-03,0.829276D-03, &
        0.836237D-03,0.843243D-03,0.850306D-03,0.857418D-03,0.864568D-03, &
        0.871777D-03,0.879027D-03,0.886321D-03,0.893664D-03,0.901059D-03, &
        0.908518D-03,0.916017D-03,0.923583D-03,0.931201D-03,0.938865D-03, &
        0.946581D-03,0.954347D-03,0.962169D-03,0.970049D-03,0.977973D-03, &
        0.985961D-03,0.993995D-03,0.100210D-02,0.101027D-02,0.101850D-02, &
        0.102677D-02,0.103510D-02,0.104349D-02,0.105192D-02,0.106040D-02, &
        0.106897D-02,0.107761D-02,0.108635D-02,0.109510D-02,0.110390D-02, &
        0.111274D-02,0.112163D-02,0.113063D-02,0.113971D-02,0.114886D-02, &
        0.115806D-02,0.116732D-02,0.117663D-02,0.118600D-02,0.119543D-02, &
        0.120495D-02,0.121453D-02,0.122419D-02,0.123391D-02,0.124371D-02, &
        0.125356D-02,0.126348D-02,0.127347D-02,0.128353D-02,0.129367D-02, &
        0.130387D-02,0.131415D-02,0.132451D-02,0.133493D-02,0.134544D-02, &
        0.135602D-02,0.136666D-02,0.137737D-02,0.138816D-02,0.139902D-02, &
        0.140998D-02,0.142099D-02,0.143209D-02,0.144327D-02,0.145453D-02, &
        0.146587D-02,0.147729D-02,0.148878D-02,0.150035D-02,0.151200D-02, &
        0.152373D-02,0.153554D-02,0.154745D-02,0.155945D-02,0.157155D-02, &
        0.158371D-02,0.159591D-02,0.160822D-02,0.162062D-02,0.163313D-02, &
        0.164572D-02,0.165839D-02,0.167115D-02,0.168400D-02,0.169693D-02, &
        0.170996D-02,0.172306D-02,0.173628D-02,0.174958D-02,0.176297D-02, &
        0.177646D-02,0.179003D-02,0.180371D-02,0.181748D-02,0.183135D-02, &
        0.184530D-02,0.185935D-02,0.187350D-02,0.188775D-02,0.190209D-02 &
        ])


!!check if rots is in range of data set,if not in range 2~1003,return -1
        if(rots < rotsMin.or.rots > rotsMax)then
            csen = -1D0
            return
        end if
!!calculate csen by Interpolation
        x = LOG( rots / 2D0 ) / LOG( 1.02D0) + 1D0
        i = floor(x)
        csen = ( sect(i+1) - sect(i) ) * (x - i) + sect(i)


        return
        end



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine ltof(ii)
!!      Moves ii-th particle (lepton) in PYJETS to first position.
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        PARAMETER (KSZJ=80000)
        COMMON/PYJETS/N,NONJ,K(KSZJ,5),P(KSZJ,5),V(KSZJ,5)
        common/sa5/kfmax,kfaco(100),numb(100),numbs(100),non5, &
         disbe(100,100)
        dimension kk(5), pp(5), vv(5)


        do jj=1,5
            kk(jj) = K(ii,jj)
            pp(jj) = P(ii,jj)
            vv(jj) = V(ii,jj)
        end do
!       Moves particle list (PYJETS) one step forward from ii-1 to 1.
        do j1 = ii-1, 1, -1
            do jj=1,5,1
                K( j1+1, jj ) = K(j1,jj)
                P( j1+1, jj ) = P(j1,jj)
                V( j1+1, jj ) = V(j1,jj)
            end do
        end do
        do jj=1,5
            K(1,jj) = kk(jj)
            P(1,jj) = pp(jj)
            V(1,jj) = vv(jj)
        end do
        ! now first particle in PYJETS is e-
        do i1=1,kfmax,1
            numbs(i1) = numbs(i1) + 1
        end do


        return
        end



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine simulate_NN_in_framework_A
!!      Independent simulation of NN collisions in A-framework.
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        PARAMETER (KSZJ=80000)
        PARAMETER (NSIZE=750000)
        COMMON/PYJETS/N,NPAD,K(KSZJ,5),P(KSZJ,5),V(KSZJ,5)
        common/sa1/kjp21,non1,bp,iii,neve,nout,nosc
        common/sa2/nsa,non2,ksa(kszj,5),psa(kszj,5),vsa(kszj,5)
        common/sa10/csnn,cspin,cskn,cspipi,cspsn,cspsm,rcsit,ifram, &
         iabsb,iabsm,i_sigma_AQM,ajpsi,csspn,csspm,csen
        common/sa16/x_ratio,decpro,dni(10),dpi(10),edi(10),bmin,bmax, &
         bar(10),abar(10),barf(10),abarf(10),emin(10),eminf(10), &
         eplu(10),epluf(10),nmax
        common/ctllist_t/ lc(nsize,5), tc(nsize), tw(nsize)
        common/syspar/ipden,itden,suppm,suptm,suppc,suptc,r0p,r0t, &
         nap,nat,nzp,nzt,pio
!       For the charge and 4-momentum check.
        common/cp_check/ p_init(4), p_final(4)
!       For the simulation control.
        COMMON/SA1_PY8/ i_mode, i_tune, KF_woDecay(100), &
               KF_proj, KF_targ, win, energy_B, psno, b_min, b_max
        dimension pi(4),pj(4),b(3)


        time = 0D0


!-------------------------------------------------------------------------------
!------------------------   NN Low Energy A-framework   ------------------------
!       Initialization for a nucleon-nucleon collision in 'PYJETS'
        N = 2
!       Initiation array 'K'.
        K(1,1) = 2   ! A
        K(2,1) = 1   ! V

        if(nzp == 1)then
            K(1,2) = 2212
        else if(nzp == -1)then
            K(1,2) = -2212
        else if(nzp == 0)then
            K(1,2) = 2112
        end if

        if(nzt == 1)then
            K(2,2) = 2212
        else if(nzt == -1)then
            K(2,2) = -2212
        else if(nzt == 0)then
            K(2,2) = 2112
        end if

!       Initialization in spatial phase space.
!       v = 0D0

!       Initialization in momentum phase space.

        ep1 = 0D0
        ep2 = 0D0
        et1 = 0D0
        et2 = 0D0
        pp1 = 0D0
        pp2 = 0D0
        pt1 = 0D0
        pt2 = 0D0
        ! Sets the default frame as CMS.
        if(ifram > 2) ifram = 1

        if(ifram == 0)then
            pp1 = win     ! momentum of projetile (if proton)
            pt1 = 1D-20   ! momentum of target    (if proton)
            pp2 = win     ! momentum of projetile (if neutron)
            pt2 = 1D-20   ! momentum of target    (if neutron)
            pm2 = PYMASS(2212)**2       ! square mass of proton
            ep1 = SQRT(pp1*pp1 + pm2)   ! energy of projetile (if proton)
            et1 = SQRT(pt1*pt1 + pm2)   ! energy of target    (if proton)
            pm2 = PYMASS(2112)**2       ! square mass of neutron
            ep2 = SQRT(pp2*pp2 + pm2)   ! energy of projetile (if neutron)
            et2 = SQRT(pt2*pt2 + pm2)   ! energy of target    (if neutron)
        endif

        if(ifram == 1)then
            ep1 = 0.5D0*win   ! energy of projetile (if proton)
            et1 = ep1         ! energy of target    (if proton)
            ep2 = 0.5D0*win   ! energy of projetile (if neutron)
            et2 = ep2         ! energy of target    (if neutron)
            pm2 = PYMASS(2212)**2          ! square mass of proton
            pp1 =  SQRT( ep1*ep1 - pm2 )   ! momentum of projetile (if proton)
            pt1 = -SQRT( et1*et1 - pm2 )   ! momentum of target    (if proton)
            pm2 = PYMASS(2112)**2          ! square mass of nucleon
            pp2 =  SQRT( ep2*ep2 - pm2 )   ! momentum of projetile (if neutron)
            pt2 = -SQRT( et2*et2 - pm2 )   ! momentum of target    (if neutron)
        endif

        if(ifram == 2)then
            eA = win
            eB = energy_B
            ep1 = eA   ! energy of projetile particle (if proton)
            et1 = eB   ! energy of target particle    (if proton)
            ep2 = eA   ! energy of projetile particle (if neutron)
            et2 = eB   ! energy of target particle    (if neutron)
            if( ep1 < PYMASS(2212) )then
                if( iii == 1 )then
                    write(22,*)
                    write(22,*) "Warning! win(energy_A) < m_pronton " // &
                                "in the back-to-back collision " // &
                                "frame of A-framework. " // &
                                "Re-evaluated the energy " // &
                                "according to m_proton."
                    write(22,*)
                end if
                ep1 = PYMASS(2212)
            end if
            if( et1 < PYMASS(2212) )then
                if( iii == 1 )then
                    write(22,*)
                    write(22,*) "Warning! energy_B < m_pronton " // &
                                "in the back-to-back collision " // &
                                "frame of A-framework. " // &
                                "Re-evaluated the energy " // &
                                "according to m_proton."
                    write(22,*)
                end if
                et1 = PYMASS(2212)
            end if
            if( ep2 < PYMASS(2112) )then
                if( iii == 1 )then
                    write(22,*)
                    write(22,*) "Warning! win(energy_A) < m_neutron " // &
                                "in the back-to-back collision " // &
                                "frame of A-framework. " // &
                                "Re-evaluated the energy " // &
                                "according to m_neutron."
                    write(22,*)
                end if
                ep2 = PYMASS(2112)
            end if
            if( et2 < PYMASS(2112) )then
                if( iii == 1 )then
                    write(22,*)
                    write(22,*) "Warning! energy_B < m_neutron " // &
                                "in the back-to-back collision " // &
                                "frame of A-framework. " // &
                                "Re-evaluated the energy " // &
                                "according to m_neutron."
                    write(22,*)
                end if
                et2 = PYMASS(2112)
            end if
            pm2 = PYMASS(2212)**2          ! square mass of proton
            pp1 =  SQRT( ep1*ep1 - pm2 )   ! momentum of projetile       (if proton)
            pt1 = -SQRT( et1*et1 - pm2 )   ! momentum of target particle (if proton)
            pm2 = PYMASS(2112)**2          ! square mass of nucleon
            pp2 =  SQRT( ep2*ep2 - pm2 )   ! momentum of projetile       (if neutron)
            pt2 = -SQRT( et2*et2 - pm2 )   ! momentum of target particle (if neutron)
        endif

!       Four momenta of projectile particle.
        P(1,1) = 0D0
        P(1,2) = 0D0
        ! projectile particle is p (pbar)
        if( ABS(nzp) == 1 )then
            P(1,3) = pp1
            P(1,4) = ep1
            P(1,5) = PYMASS(2212)
        ! projectile particle is neutron
        else if(nzp == 0)then
            P(1,3) = pp2
            P(1,4) = ep2
            P(1,5) = PYMASS(2112)
        end if
        ! The initial total 4-momentum.
        p_init(3) = p_init(3) + P(1,3)
        p_init(4) = p_init(4) + P(1,4)

!       Four momenta of target particle.
        P(2,1) = 0D0
        P(2,2) = 0D0
        ! target particle is p (pbar)
        if( ABS(nzt) == 1 )then
            P(2,3) = pt1
            P(2,4) = et1
            P(2,5) = PYMASS(2212)
        ! target particle is neutron
        else if( nzt == 0 )then
            P(2,3) = pt2
            P(2,4) = et2
            P(2,5) = PYMASS(2112)
        end if
        ! The initial total 4-momentum.
        p_init(3) = p_init(3) + P(2,3)
        p_init(4) = p_init(4) + P(2,4)
        if( ABS(P(1,3)) < 1D-10 .AND. ABS(P(2,3)) < 1D-10 )then
            write(22,*)
            write(22,*) "Abort! Too low win(energy_A) and energy_B!"
            write(22,*)
            stop
        end if

!       Initiatlizeion finished.

!       full_events_history of OSC1999A
        call oscar(0)

!       'PYJETS' to 'sa2'
        nsa = N
        do i2=1,5
            do i1=1,N
                ksa(i1,i2) = K(i1,i2)
                psa(i1,i2) = P(i1,i2)
                vsa(i1,i2) = V(i1,i2)
            end do
        end do

        do i2=1,4
            ! four momentum of one colliding particle
            pi(i2) = psa(1,i2)
            ! four momentum of another one colliding particle
            pj(i2) = psa(2,i2)
        enddo

        if(ifram == 0 .OR. ifram == 2)then
!       boost to CMS frame of colliding pair
            do i2=1,3
                b(i2) = ( pi(i2) + pj(i2) ) / ( pi(4) + pj(4) )
            enddo
            ilo=0
            call lorntz(ilo,b,pi,pj)
        end if
        if( ifram == 1 ) b = 0D0    ! ifram=1, b(i)=0 & \gamma=1

!       following is for case of NN collision only
        ww = rcsit
!       the cross section ratio of (ela.)/tot =1- rcsit (rcsit=inela./tot)
        rrlu = PYR(1)
        inorex = 1

!***************************   Elastic Scattering   ****************************
        if(rrlu > ww)then
!       perform nucleon-nucleon (NN) elastic scattering (in parini.f90).
            call coelas_nn(1,2,pi,pj,time,b,2)
            inorex=1
!***************************   Elastic Scattering   ****************************

!**************************   Inelastic Scattering   ***************************
        else if(rrlu <= ww)then
!       perform NN inelastic scattering
!     p: 2212; n:2112; delta0: 2114; delta+: 2214; delta-: 1114; delta++: 2224
!       consider following 2->2 NN ielastic channels
!     1       p + p to delta+ + p
!     2       p + p to delta++ + n
!     3       p + n to delta+ + n
!     4       p + n to delta0 + p
!     5       n + n to delta0 + n
!     6       n + n to delta- + p
!       reverse scattering is not considered

!       inorex: endothermic or exothermic
!       endothermic reaction, inorex=1 & ela. (cf. 'coinel_nn')
!       exothermic reaction, inorex=2 & inela. (cf. 'coinel_nn')

!       give flavor code to inelastically scattered particles
!           inorex = 1
            rpy = PYR(1)

            if(nzp == 1 .and. nzt  == 1)then   ! pp   !!
                if(rpy > 0.5D0 )then
                    ksa(1,2) = 2214
                    ksa(2,2) = 2212
                else
                    ksa(1,2) = 2224
                    ksa(2,2) = 2112
                end if

            else if(nzp == 0 .and. nzt  == 0)then   ! nn   !!
                if(rpy > 0.5D0)then
                    ksa(1,2) = 2114
                    ksa(2,2) = 2112
                else
                    ksa(1,2) = 1114
                    ksa(2,2) = 2212
                end if

            else if( (nzp == 1 .and. nzt  == 0) .or. &   ! pn   !!
                     (nzt == 1 .and. nzp  == 0) )then
                if(rpy  >  0.5D0)then
                    ksa(1,2) = 2214
                    ksa(2,2) = 2112
                else
                    ksa(1,2) = 2114
                    ksa(2,2) = 2212
                end if
            end if   !!

            KF_1 = ksa(1,2)
            KF_2 = ksa(2,2)

            call adjudg_nn(KF_1,KF_2,pi,pj,inorex)

            if(inorex == 1)then
                call coelas_nn(1,2,pi,pj,time,b,2)
            else if(inorex == 2)then
!       give four momentum to inelastically scattered particles
                call coinel_nn(1,2,KF_1,KF_2,pi,pj)
!               call padecy_A(1,2)
            end if
        end if
!**************************   Inelastic Scattering   ***************************

        if( inorex == 1 )then
!       'sa2' to 'PYJETS'
            N = nsa
            do i2=1,5,1
                do i1=1,nsa,1
                    K(i1,i2) = ksa(i1,i2)
                    P(i1,i2) = psa(i1,i2)
                    V(i1,i2) = vsa(i1,i2)
                end do
            end do
        else if( inorex == 2 )then
            call padecy_A(1,2)
        end if

!       Removes photons from "PYJETS" to "sgam" if any.
        KF_gamma_decay = 22
        call remo_gam( KF_gamma_decay )

!       low energy NN collision A-framework finished
!------------------------   NN Low Energy A-framework   ------------------------
!-------------------------------------------------------------------------------


        return
        end



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine simulate_hh_of_AB_in_framework_A( l, l1, pi, pj, &
                                            time, b, ss, inorex, jorn )
!!      The hh collision simulation of an AB collision in A-framework.
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        PARAMETER (KSZJ=80000)
        PARAMETER (NSIZE=750000)
        COMMON/PYJETS/N,NPAD,K(KSZJ,5),P(KSZJ,5),V(KSZJ,5)
        common/sa1/kjp21,non1,bp,iii,neve,nout,nosc
        common/sa2/nsa,non2,ksa(kszj,5),psa(kszj,5),vsa(kszj,5)
        common/sa6_p/ithroq_p,ithrob_p,ich_p,non6_p,throe_p(4)
        common/sa10/csnn,cspin,cskn,cspipi,cspsn,cspsm,rcsit,ifram, &
         iabsb,iabsm,i_sigma_AQM,ajpsi,csspn,csspm,csen
        common/sa16/x_ratio,decpro,dni(10),dpi(10),edi(10),bmin,bmax, &
         bar(10),abar(10),barf(10),abarf(10),emin(10),eminf(10), &
         eplu(10),epluf(10),nmax
        common/ctllist_t/ lc(nsize,5), tc(nsize), tw(nsize)
        common/syspar/ipden,itden,suppm,suptm,suppc,suptc,r0p,r0t, &
         nap,nat,nzp,nzt,pio
        COMMON/SA1_PY8/ i_mode, i_tune, KF_woDecay(100), &
               KF_proj, KF_targ, win, energy_B, psno, b_min, b_max
        dimension pi(4),pj(4),b(3)


        kfa  = ksa(l,2)
        kfb  = ksa(l1,2)
        ikfa = ABS( kfa )
        ikfb = ABS( kfb )


!-------------------------------------------------------------------------------
!-----------------------------   N & Delta & pi   ------------------------------
        if( (ikfa == 2212.or.ikfa == 2112.or.ikfa == 211.or.kfa == 1114 &
          .or.kfa == 2114.or.kfa == 2214.or.kfa == 2224).and. &
            (ikfb == 2212.or.ikfb == 2112.or.ikfb == 211.or.kfb == 1114 &
          .or.kfb == 2114.or.kfb == 2214.or.kfb == 2224) )then
            nsa0 = nsa
            ww = rcsit
!       active  for ela. scattering only
!           ww = 0D0
!       active  for inela. scattering only
!           ww = 1D0
!       the cross section ratio of (ela.)/tot =1- rcsit
            rrlu = PYR(1)

!***************************   Elastic Scattering   ****************************
        if(rrlu > ww)then
!       perform two-body elastic scattering
            call coelas_nn(l,l1,pi,pj,time,b,nsa0)
            inorex = 1
!***************************   Elastic Scattering   ****************************

!**************************   Inelastic Scattering   ***************************
        elseif(rrlu <= ww)then
!       inelastic scattering
!     p: 2212; n:2112; delta0: 2114; delta+: 2214; delta-: 1114; delta++: 2224
!     pi+: 211; pi-: -211; pi0: 111
!       NN ielas. channels (N including delta):
!     1       p + p to delta+ + p
!     2       p + p to delta++ + n
!     3       p + n to delta+ + n
!     4       p + n to delta0 + p
!     5       n + n to delta0 + n
!     6       n + n to delta- + p
!     7       delta+ + p to p + p
!     8       delta+ + n to p + n
!     9       delta0 + p to p + n
!     10      delta0 + n to n + n
!     11      delta++ + n to p + p
!     12      delta- + p to n + n

!       piN inelastic scattering:
!       consider following 2->2 pi-nucleon ielas. channels
!     13      pion- + p to delta- + pion+
!     14      pion- + p to rho0 + n
!     15      pion- + p to rho- + p
!     16      pion- + p to delta+ + pion-
!     17      pion- + p to delta0 + pion0
!     18      pion- + n to delta- + pion0
!     19      pion- + n to rho- + n
!     20      pion- + n to delta0 + pion-
!     21      pion+ + p to delta++ + pion0
!     22      pion+ + p to delta+ + pion+
!     23      pion+ + p to rho+ + p
!     24      pion+ + n to delta++ + pion-
!     25      pion+ + n to delta0 + pion+
!     26      pion+ + n to delta+ + pion0
!     27      pion+ + n to rho0 + p
!     28      pion+ + n to rho+ + n
!       assume pion0 decay instantly, following processes is not considered
!     29      pion0 + p to delta0 + pion+
!     30      pion0 + p to delta++ + pion-
!     31      pion0 + p to rho+ + n
!     32      pion0 + p to rho0 + p
!     33      pion0 + p to delta+ + pion0
!     34      pion0 + n to delta+ + pion-
!     35      pion0 + n to delta- + pion+
!     36      pion0 + n to delta0 + pion0

!     piN 2->1 resonance production:
!     37      pion+ + p to delta++

!       inorex: endothermic or exothermic
!       thres=ss-am3-am4, threshold energy, effective in 'adjudg_nn'
!       thres:  <= 0, treat as endothermic reaction, inorex=1, & ela.
!       thres:  > 0, treat as exothermic reaction and inorex=2, inela.
!       jorn=0 (3): inelastically scattered particles, or resonance particle,
!        not join (join) reconstruction of hh collision pair
        inorex = 0
!       give flavor code to inelastically scattered particles
        rpy  = PYR(1)
        rpy1 = PYR(1)

!-----------------------------   N+N Scattering   ------------------------------
!       treat NN->(delta)N
        if(kfa == 2212.and.kfb == 2212)then   ! pp   !!
            if(rpy  >  0.5)then
                ksa(l,2)=2214   ! delta+
                ksa(l1,2)=2212
            else
                ksa(l,2)=2224   ! delta++
                ksa(l1,2)=2112
            endif
            jorn=3

        elseif( (kfa == 2212.and.kfb == 2112).or. &   ! pn   !!
                (kfb == 2212.and.kfa == 2112) )then
            if(rpy  >  0.5)then
                ksa(l,2)=2214   ! delta+
                ksa(l1,2)=2112
            else
                ksa(l,2)=2114   ! delta0
                ksa(l1,2)=2212
            endif
            jorn=3

        elseif(kfa == 2112.and.kfb == 2112)then   ! nn   !!
            if(rpy  >  0.5)then
                ksa(l,2)=2114   !  delta0
                ksa(l1,2)=2112
            else
                ksa(l,2)=1114   ! delta-
                ksa(l1,2)=2212
            endif
            jorn=3
!-----------------------------   N+N Scattering   ------------------------------

!---------------------------   Delta+N Scattering   ----------------------------
!       treat (delta)N->NN
        elseif( (kfa == 2214.and.kfb == 2212).or. &   ! (delta+)p   !!
                (kfb == 2214.and.kfa == 2212) )then
            ksa(l,2)=2212
            ksa(l1,2)=2212
            rpy1 = -1D0   ! call sa2pyj(l,l1)
            jorn=3

        elseif( (kfa == 2214.and.kfb == 2112).or. &   ! (delta+)n   !!
                (kfb == 2214.and.kfa == 2112) )then
            ksa(l,2)=2212
            ksa(l1,2)=2112
            rpy1 = -1D0   ! call sa2pyj(l,l1)
            jorn=3

        elseif( (kfa == 2114.and.kfb == 2212).or. &   ! (delta0)p   !!
                (kfb == 2114.and.kfa == 2212) )then
            ksa(l,2)=2212
            ksa(l1,2)=2112
            rpy1 = -1D0   ! call sa2pyj(l,l1)
            jorn=3

        elseif( (kfa == 2114.and.kfb == 2112).or. &   ! (delta0)n   !!
                (kfb == 2114.and.kfa == 2112) )then
            ksa(l,2)=2112
            ksa(l1,2)=2112
            rpy1 = -1D0   ! call sa2pyj(l,l1)
            jorn=3

        elseif( (kfa == 2224.and.kfb == 2112).or. &   ! (delta++)n   !!
                (kfb == 2224.and.kfa == 2112) )then
            ksa(l,2)=2212
            ksa(l1,2)=2212
            rpy1 = -1D0   ! call sa2pyj(l,l1)
            jorn=3

        elseif( (kfa == 1114.and.kfb == 2212).or. &   ! (delta-)p   !!
                (kfa == 1114.and.kfb == 2212) )then
            ksa(l,2)=2112
            ksa(l1,2)=2112
            rpy1 = -1D0   ! call sa2pyj(l,l1)
            jorn=3
!---------------------------   Delta+N Scattering   ----------------------------

!-----------------------------   pi+N Scattering   -----------------------------
!       treat (pi)N->... (which is one order higher than NN->... in
!        approximation)
        elseif( (kfa == -211.and.kfb == 2212).or. &   ! (pi-)p   !!
                (kfa == 2212.and.kfb == -211) )then
            if(rpy  <=  0.2)then
                ksa(l,2)=1114   ! delta-
                ksa(l1,2)=211
            elseif(rpy  >  0.2 .and. rpy  <=  0.4)then
                ksa(l,2)=113   ! rho0
                ksa(l1,2)=2112
            elseif(rpy  >  0.4 .and. rpy  <=  0.6)then
                ksa(l,2)=-213   ! rho-
                ksa(l1,2)=2212
            elseif(rpy  >  0.6 .and. rpy  <=  0.8)then
                ksa(l,2)=2214   ! delta+
                ksa(l1,2)=-211
            elseif(rpy  >  0.8)then
                ksa(l,2)=2114   ! delta0
                ksa(l1,2)=111
            endif
            rpy1 = 1D0 + decpro   ! call padecy_A(l,l1)
            jorn=0

        elseif( (kfa == -211.and.kfb == 2112).or. &   ! (pi-)n   !!
                (kfb == -211.and.kfa == 2112) )then
            if(rpy  <=  0.3333)then
                ksa(l,2)=1114   ! delta-
                ksa(l1,2)=111
            elseif(rpy  >  0.3333 .and. rpy  <=  0.6666)then
                ksa(l,2)=-213   ! rho-
                ksa(l1,2)=2112
            elseif(rpy  >  0.6666)then
                ksa(l,2)=2114   ! delta0
                ksa(l1,2)=-211
            endif
            rpy1 = 1D0 + decpro   ! call padecy_A(l,l1)
            jorn=0

        elseif( (kfa == 211.and.kfb == 2212).or. &   ! (pi+)p   !!
                (kfb == 211.and.kfa == 2212) )then
!           percen   = 0.33333333   ! for the case without resonance production
            percen   = 0.00000001   ! for the case with resonance production only
!           percen=0.25  ! for the case with both
            percen_1 = percen
            percen_2 = percen_1 + percen
            percen_3 = percen_2 + percen
!       seting percen= 0.33333333 for the case without resonance production
!       seting percen= 0.00000001 for the case with resonance production only
            jorn=0
            if(rpy  <=  percen_1)then
                ksa(l,2)=2224   ! delta++
                ksa(l1,2)=111
                rpy1 = 1D0 + decpro   ! call padecy_A(l,l1)
            elseif(rpy  >  percen_1 .and. rpy  <=  percen_2)then
                ksa(l,2)=2214   ! delta+
                ksa(l1,2)=211
                rpy1 = 1D0 + decpro   ! call padecy_A(l,l1)
            elseif(rpy  >  percen_2 .and. rpy  <=  percen_3)then
                ksa(l,2)=213   ! rho+
                ksa(l1,2)=2212
                rpy1 = 1D0 + decpro   ! call padecy_A(l,l1)
            elseif(rpy  >  percen_3)then
                ksa(l,2)=2224   ! delta++
                ksa(l1,1)=0
                amas=pymass(2224)
                thres=ss-amas
                if(thres <= 0.) inorex=1
                if(thres > 0.) inorex=2
                if(inorex == 1)then   ! elas.
                    call coelas_nn(l,l1,pi,pj,time,b,nsa0)
                endif
                ! resonance production
                if(inorex == 2)then
!       give four momentum to resonance particle
                    do i2=1,3
                        psa(l,i2) = psa(l,i2) + psa(l1,i2)
                    enddo
                    psal4=amas*amas+psa(l,1)*psa(l,1)+psa(l,2)*psa(l,2)+ &
                          psa(l,3)*psa(l,3)
                    if(psal4 <= 1.d-20) psal4=1.d-20
                        psal4=dsqrt(psal4)
                        psa(l,4)=psal4
                        ! energy, throw away
                        throe_p(4)=throe_p(4)+(ss-psal4)
                    do i2=1,4
                        pi(i2)=psa(l,i2)
                        pj(i2)=0.
                        psa(l1,i2)=0.
                    enddo
                    ! delta++ not decay immediately
                    if(rpy1 <= decpro)then
                        call sa2pyj(l,0)
                    ! delta++ decays immediately
                    else
                        call padecy_1(l)
                    endif
                endif
                ! Direct return for this special channel.
                return
            endif

        elseif( (kfa == 211.and.kfb == 2112).or. &   ! (pi+)n   !!
                (kfb == 211.and.kfa == 2112) )then
            if(rpy  <=  0.2)then
                ksa(l,2)=2224   ! delta++
                ksa(l1,2)=-211
            elseif(rpy  >  0.2 .and. rpy  <=  0.4)then
                ksa(l,2)=2114   ! delta0
                ksa(l1,2)=211
            elseif(rpy  >  0.4 .and. rpy  <=  0.6)then
                ksa(l,2)=2214   ! delta+
                ksa(l1,2)=111
            elseif(rpy  >  0.6 .and. rpy  <=  0.8)then
                ksa(l,2)=113   ! rho0
                ksa(l1,2)=2212
            elseif(rpy  >  0.8)then
                ksa(l,2)=213   ! rho+
                ksa(l1,2)=2112
            endif
            rpy1 = 1D0 + decpro   ! call padecy_A(l,l1)
            jorn=0
        endif
!-----------------------------   pi+N Scattering   -----------------------------

!       More general code.
            KF_1 = ksa(l,2)
            KF_2 = ksa(l1,2)
            if(inorex > 0) call adjudg_nn(KF_1,KF_2,pi,pj,inorex)
!       adjudge current two-body inelas. scattering to be really treated as
!        elas. (if inorex=1) or inelas. (if inorex=2) due to threshol energy
!       For channels that have not yet been implemented. Treats them as
!        elastic scatterigs.
            if(inorex == 0) inorex = 1
!       perform two-body elastic scattering
            if(inorex == 1) call coelas_nn(l,l1,pi,pj,time,b,nsa0)
            if( inorex == 2 )then
                call coinel_nn(l,l1,KF_1,KF_2,pi,pj)
!       treat inelastic scattering of current hh collision pair,
!       i.e. give four momentum to scattered particles
                if(rpy1 <= decpro)then
                    call sa2pyj(l,l1)
!       a part of update particle list ('sa2' to 'PYJETS') after inela.
!        scattering in case of outgoing channel with \Delta particle
                else
                    call padecy_A(l,l1)
!       a part of update particle list ('sa2' to 'PYJETS') after inela.
!        scattering in case of outgoing channel without \Delta particle
                end if
            end if
        end if
!**************************   Inelastic Scattering   ***************************

        end if
!-----------------------------   N & Delta & pi   ------------------------------
!-------------------------------------------------------------------------------


        return
        end



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine adjudg_nn(kl,kl1,pi,pj,inorex)
!!      Adjudges current hadron-hadron collision pair to be treated as
!!       elas. or inelas. two-body scattering in A-framework.
!       kl (kl1): flavor code of scattered particle
!       pi (pj): four momentum of scattering particle before scatering
!       pi (pj): four momentum of scattered particle after scatering
!       inorex: endothermic or exothermic
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        PARAMETER (KSZJ=80000)
        common/sa2/nsa,non2,ksa(kszj,5),psa(kszj,5),vsa(kszj,5)
        dimension pi(4),pj(4)


        am3=pymass(kl)   ! mass of one scattered particle
        am4=pymass(kl1)   ! mass of another one scattered particle
        ss=pi(4)+pj(4)
!       ss: cms energy (total energy) of current hh collision pair
        thres=ss-am3-am4
        if(thres <= 0.) inorex=1
!       thres:  <= 0, treat as endothermic reaction, inorex=1 & ela.
        if(thres > 0.) inorex=2
!       thres:  > 0, treat as exothermic reaction, inorex=2 & inela.


        return
        end



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine coelas_nn(l,l1,pi,pj,time,b,nsa0)
!!      Treats two-body elas. scattering of current hh collision pair
!!       in A-framework.
!       input and output messages are all in 'sa2'
!       l (l1): order number in 'sa2' of current hh collision pair
!       pi (pj): four momentum of scattering particle before scatering
!       pi (pj): four momentum of scattered particle after scatering
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        PARAMETER (KSZJ=80000,NSIZE=750000)
        common/sa2/nsa,non2,ksa(kszj,5),psa(kszj,5),vsa(kszj,5)
        common/ctllist_t/ lc(nsize,5), tc(nsize), tw(nsize)
        common/syspar/ipden,itden,suppm,suptm,suppc,suptc,r0p,r0t, &
         nap,nat,nzp,nzt,pio
        dimension pi(4),pj(4),b(3)
        ss=pi(4)+pj(4)
!       ss: cms energy (total energy) of current hh collision pair


!       perform elastic scattering
        call coelas(l,l1,ss,pi,pj)

        call updple(l,l1,b,pi,pj)
!       update particle list, pi and pj have been boosted back
!        to Lab or cms frame of nucleus-nucleus collision system
        if( ipden == 0.and.itden == 0 ) return   ! for NN collision system

        call updatl_nn(l,l1,time,3,nsa0,i_mode)   ! 250423
!       update collision list after ela. scattering
!       jorn=0 (3): scattered particles not join (join) reconstruction
!        of hh collision pair
!       note: CME is not considered in the scattering channels of low
!        energy A-framework.


        return
        end



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine coinel_nn(l,l1,kl,kl1,pi,pj)
!!      Treats inelastic scattering of current hh collision pair,
!!       i.e. simulate four momentum of scattered particles in inela.
!!       hh collision according formula of (4. 33) in book of Ben-Hao Sa,
!!       etc., 'Simulation Physics for High Energy Nucleus Collisions'
!!       in A-framework
!       l (l1): order number in 'sa2' of current hadron-hadron collision pair
!       kl (kl1): flavor code of scattered particle
!       pi (pj): four momentum of scattering particle before scatering
!       pi (pj): four momentum of scattered particle after scatering
!       inorex: endothermic or exothermic
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        PARAMETER (KSZJ=80000)
        common/sa2/nsa,non2,ksa(kszj,5),psa(kszj,5),vsa(kszj,5)
        dimension pi(4),pj(4)


        am3=pymass(kl)   ! mass of one scattered particle
        am4=pymass(kl1)   ! mass of another one scattered particle
        ss=pi(4)+pj(4)
!       ss: cms energy (total energy) of current h-h collision pair
        pp=(ss*ss-(am3+am4)**2)*(ss*ss-(am3-am4)**2)/(4.*ss*ss)
        if(pp < 0.) pp=1.e-10
        pp=sqrt(pp)
        pi(4)=(ss*ss+am3**2-am4**2)/(2.*ss)
!       energy of one particle (between two) after scattering
        fis=2.*3.1415926*pyr(1)
        cfis=cos(fis)
        sfis=sin(fis)
        csita=2*pyr(1)-1.
        ssita=dsqrt(1.-csita**2)
        pi(1)=pp*ssita*cfis
        pi(2)=pp*ssita*sfis
        pi(3)=pp*csita
        do i=1,3
            pj(i)=0.-pi(i)
        enddo
        pj(4)=ss-pi(4)
        do i2=1,4
            psa(l,i2)=pi(i2)   ! four momentum of one collided particle
            psa(l1,i2)=pj(i2)   ! four momentum of another one collided particle
        enddo


        return
        end



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine padecy_A(l,l1)
!!      Deacy of a hardron.
!!      part of 'sa2' to 'PYJETS', i.e. a part of updating particle list after
!!       inela. scattering in case of outgoing channel without \Delta particle
!!       in A-framework
!       l: order # of colliding hadron in 'sa2'
!       l1: order # of another colliding hadron in 'sa2'
!       input messages in 'sa2'
!       output messages compose 'PYJETS'
!       Double precision and integer declarations.
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        INTEGER PYCOMP
        PARAMETER (KSZJ=80000)
        COMMON/PYJETS/N,NPAD,K(KSZJ,5),P(KSZJ,5),V(KSZJ,5)
        COMMON/PYDAT3/MDCY(500,3),MDME(8000,2),BRAT(8000),KFDP(8000,5)
        common/sa2/nsa,non2,ksa(kszj,5),psa(kszj,5),vsa(kszj,5)


        N = 1
        K(1,1) = 1
        kf = ksa(l,2)
        K(1,2) = kf
        do i2=1,4
            P(1,i2) = psa(l,i2)
        end do
        P(1,5) = PYMASS(kf)
        MDCY( PYCOMP(kf), 1 ) = 1

        call PYDECY(1)
        call PYEDIT(1)

        if(l1 == 0) return

        N = 3
        K(3,1) = 1
        K(3,2) = ksa(l1,2)
        do i1=1,5
            p(3,i1) = psa(l1,i1)
        end do


        return
        end



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine padecy_delt(l)
!       Deacy of a hardron in A-framework.
!       l: order # of colliding hadron in 'delt'
!       input messages in 'delt'
!       output messages compose 'PYJETS'
!       Double precision and integer declarations.
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        INTEGER PYCOMP
        PARAMETER (KSZJ=80000)
        COMMON/PYJETS/N,NPAD,K(KSZJ,5),P(KSZJ,5),V(KSZJ,5)
        COMMON/PYDAT3/MDCY(500,3),MDME(8000,2),BRAT(8000),KFDP(8000,5)
        common/delt/nsa,non2,ksa(kszj,5),psa(kszj,5),vsa(kszj,5)
        dimension r(4)


        r(1) = vsa(l,1)
        r(2) = vsa(l,2)
        r(3) = vsa(l,3)
        r(4) = vsa(l,4)
        rlu  = PYR(1)

        N = 1
        K(1,1) = 1
        kf = ksa(l,2)
        K(1,2) = kf
        do i2=1,4
            P(1,i2) = psa(l,i2)
        end do
        P(1,5) = PYMASS(kf)
        MDCY( PYCOMP(kf), 1 ) = 1

        call PYDECY(1)
        call PYEDIT(1)

        V(1,1) = r(1)
        V(1,2) = r(2)
        V(1,3) = r(3)
        V(1,4) = r(4)
        V(1,1) = r(1) * ( 1D0 + rlu )
        V(1,2) = r(2) * ( 1D0 + rlu )
        V(1,3) = r(3) * ( 1D0 + rlu )
        V(1,4) = r(4)


        return
        end



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine padecy_1(l)
!!      Deacy of a hardron in A-framework.
!       l: order # of colliding hadron in 'sa2'
!       input messages in 'sa2'
!       output messages compose 'PYJETS'
!       Double precision and integer declarations.
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        INTEGER PYCOMP
        PARAMETER (KSZJ=80000)
        COMMON/PYJETS/N,NPAD,K(KSZJ,5),P(KSZJ,5),V(KSZJ,5)
        COMMON/PYDAT3/MDCY(500,3),MDME(8000,2),BRAT(8000),KFDP(8000,5)
        common/sa2/nsa,non2,ksa(kszj,5),psa(kszj,5),vsa(kszj,5)


        N = 1
        K(1,1) = 1
        kf = ksa(l,2)
        K(1,2) = kf
        do i2=1,4
            P(1,i2) = psa(l,i2)
        end do
        P(1,5) = PYMASS(kf)
        MDCY( PYCOMP(kf), 1 ) = 1

        call PYDECY(1)
        call PYEDIT(1)


        return
        end



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine sa2pyj(l,l1)
!!      Part of 'sa2' to 'PYJETS', i.e. a part of updating particle list after
!!       inela. scattering in case of outgoing channel with \Delta particle
!!       in A-framework.
!       l: order # of colliding hadron in 'sa2'
!       l1: order # of another colliding hadron in 'sa2'
!       output hadrons compose 'PYJETS'
!       Double precision and integer declarations.
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        PARAMETER (KSZJ=80000)
        COMMON/PYJETS/N,NPAD,K(KSZJ,5),P(KSZJ,5),V(KSZJ,5)
        common/sa2/nsa,non2,ksa(kszj,5),psa(kszj,5),vsa(kszj,5)


        N = 2
        K(1,1) = 1
        K(1,2) = ksa(l,2)
        do i1=1,5
            P(1,i1) = psa(l,i1)
        end do
        if(l1 == 0)then
            N = 1
            return
        end if
        K(2,1) = 1
        K(2,2) = ksa(l1,2)
        do i1=1,5
            P(2,i1) = psa(l1,i1)
        end do


        return
        end



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine updpip_nn(l,l1)
!!      Updates particle list 'sa2' ('sbh' to 'sa2') after inela. scattering
!!       in A-framework.
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        PARAMETER (KSZJ=80000,NSIZE=750000)   ! 080223
        common/sa2/nsa,non2,ksa(kszj,5),psa(kszj,5),vsa(kszj,5)
        common/sbh/nbh,nonbh,kbh(kszj,5),pbh(kszj,5),vbh(kszj,5)
        common/ctllist_t/ lc(nsize,5), tc(nsize), tw(nsize)


!       Updates particle list 'sa2' ('sbh' to 'sa2').
        ! from call 'padecay'
        if(nbh == 3)then
            do i2=1,5
                ksa(l,i2)=kbh(1,i2)
                psa(l,i2)=pbh(1,i2)
                vsa(l,i2)=vbh(1,i2)
            enddo
            do i2=1,5
                ksa(l1,i2)=kbh(3,i2)
                psa(l1,i2)=pbh(3,i2)
                vsa(l1,i2)=vbh(3,i2)
            enddo
            nsa=nsa+1
            do i2=1,5
                ksa(nsa,i2)=kbh(2,i2)
                psa(nsa,i2)=pbh(2,i2)
                vsa(nsa,i2)=vbh(2,i2)
            enddo
        endif

        ! from 'call sa2pyj'
        if(nbh == 2)then
            do i2=1,5
                ksa(l,i2)=kbh(1,i2)
                psa(l,i2)=pbh(1,i2)
                vsa(l,i2)=vbh(1,i2)
            enddo
            do i2=1,5
                ksa(l1,i2)=kbh(2,i2)
                psa(l1,i2)=pbh(2,i2)
                vsa(l1,i2)=vbh(2,i2)
            enddo
        endif

        ! from 'call padecy_1' or 'call sa2pyj'
        if(nbh == 1)then
            do i2=1,5
                ksa(l,i2)=kbh(1,i2)
                psa(l,i2)=pbh(1,i2)
                vsa(l,i2)=vbh(1,i2)
            enddo
            nsa=nsa-1
        endif


        return
        end



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine updatl_nn(ic,jc,time,ik,nsa0,i_mode)
!!      Updates collision time list after two-body inelastic scattering
!!       in A-framework.
!       ic (jc): order number of scattering particle in 'sa2'
!       ik=3 (0): scattered hadron joins (not joins) reconstrution of hh
!        collision pair
!       time: current time
!       nsa0: nsa before two-body scattering
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        PARAMETER (KSZJ=80000)
        PARAMETER (NSIZE=750000)
        COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
        common/sa2/nsa,non2,ksa(kszj,5),psa(kszj,5),vsa(kszj,5)
        common/sa5/kfmax,kfaco(100),numb(100),numbs(100),non5, &
         disbe(100,100)
        common/sa12/ppsa(5),nchan,nsjp,sjp,taup,taujp
        common/ctllist/nctl,noinel(600),nctl0,nctlm
        common/ctllist_t/ lc(nsize,5), tc(nsize), tw(nsize)
        common/papr/t0,cspipiKK,dep,ddt,edipi,epin,ecsnn,ekn,ecspsn,ecspsm, &
         rnt,rnp,rao,rou0,vneu,vneum,ecsspn,ecsspm,ecsen
        common/syspar/ipden,itden,suppm,suptm,suppc,suptc,r0p,r0t, &
         nap,nat,nzp,nzt,pio


!       loop over old colliding pairs
        j=0
        do i=1,nctl,1
            i1 = lc(i,1)
            j1 = lc(i,2)
            if(i1 == ic .or. i1 == jc) cycle
            if(j1 == ic .or. j1 == jc) cycle
            if(i1 == j1) cycle
            if( ABS(tc(i)-time) <= ddt ) cycle
!      throw away the pair whih tc<= time
            j = j + 1
            tc(j) = tc(i)
            tw(j) = tw(i)
            do m=1,5,1
                lc(j,m) = lc(i,m)
            end do
        end do
        if(ik == 0)then
            nctl = j
            return
        end if

        nctl = j + 1

        do ii=1,2
!       loop over new generated hadron which has filled in 'sa2' at ic (jc)
            if(ii == 1) j1 = ic
            if(ii == 2) j1 = jc
            kfa  = ksa(j1,2)
            kfab = ABS(kfa)
            ! not join
            if( kfab /= 2212 .and. kfab /= 2112 .and. kfab /= 211 &
             .and. kfa /= 1114 .and. kfa /= 2114 .and. kfa /= 2214 &
             .and. kfa /= 2224 ) cycle
!       reconstruction of hadronic collision pair
!       consider only the reinteraction between NN, N(Delta), Npi+, Npi-
!       loop over old particle list
            do i=1,nsa,1
!       forbiden scattered particles colliding with each other
                if(j1 == ic .and. i == jc) cycle
                if(j1 == jc .and. i == ic) cycle
                if( j1 == i ) cycle
                i1 = i
                iflag = 0
                call rsfilt(j1,i1,iflag,i_mode)
                if(iflag == 0) cycle
                tc(nctl) = 0D0
                call tcolij(i1,j1,time,nctl,i_mode)
                if( tc(nctl) > 1D-7 ) nctl = nctl + 1
            end do
        end do

!       loop over other generated hadrons
        do ii = nsa0+1, nsa, 1
            j1   = ii
            kfa  = ksa(j1,2)
            kfab = ABS(kfa)
            ! not join
            if(kfab /= 2212.and.kfab /= 2112.and.kfab /= 211 &
             .and.kfa /= 1114.and.kfa /= 2114.and.kfa /= 2214 &
             .and.kfa /= 2224) cycle
!       reconstruction of hadronic collision pair
!       consider only the reinteraction between NN, N(Delta), Npi+, Npi-
!       loop over old particle list
            do i=1,nsa0,1
!       forbiden scattered particles colliding with each other
                if(j1 == ic .and. i == jc) cycle
                if(j1 == jc .and. i == ic) cycle
                if( j1 == i ) cycle
                i1 = i
                iflag = 0
                call rsfilt(j1,i1,iflag,i_mode)
                if(iflag == 0) cycle
                tc(nctl) = 0D0
                call tcolij(i1,j1,time,nctl,i_mode)
                if(tc(nctl) > 1D-7) nctl = nctl + 1
            end do
        end do
        nctl = nctl - 1


        return
        end



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine remo_delt
!!      Moves delta from 'PYJETS' to 'delt' in A-framework.
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        PARAMETER (KSZJ=80000)
        COMMON/PYJETS/N,NPAD,K(KSZJ,5),P(KSZJ,5),V(KSZJ,5)
        common/delt/ndel,nodel,kdel(kszj,5),pdel(kszj,5),vdel(kszj,5)


        ndel = 0
        jb = 0
201     do 202 i1 = jb+1, N, 1
            kf = K(i1,2)
            kfab = ABS(kf)
!       remove delta from 'PYJETS' to 'delt'
            if( kfab /= 1114 .and. kfab /= 2114 .and. kfab /= 2214 &
            .and. kfab /= 2224 )then
                jb = jb+1
                goto 202
            end if
            ndel = ndel + 1
            do i2=1,5
                kdel(ndel,i2) = K(i1,i2)
                pdel(ndel,i2) = P(i1,i2)
                vdel(ndel,i2) = V(i1,i2)
            end do
            if(i1 == N)then
                N = N - 1
                goto 203
            end if
!       move particle list 'PYJETS' one step downward from i1+1 to n
            do jj=1,5
                do j=i1+1,N,1
                    K(j-1,jj) = K(j,jj)
                    P(j-1,jj) = P(j,jj)
                    V(j-1,jj) = V(j,jj)
                end do
            end do
            N = N - 1
            goto 201
202     end do
203     continue


        return
        end



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine simulate_lh_and_angantyr_in_framework_BC
!!      Simulation of hh, ll & lh/hl collisions in B- and C-framework,
!!       including hA/Ah & AB collisions via Angantyr.
!!      Generates the initial state.
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        LOGICAL IS_PYTHIA8,IS_EXIST,IS_DIQUARK
        PARAMETER (KSZJ=80000,KSZJ_PY8=300000)
        COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
        COMMON/PYDAT3/MDCY(500,3),MDME(8000,2),BRAT(8000),KFDP(8000,5)
        COMMON/PYPARS/MSTP(200),PARP(200),MSTI(200),PARI(200)
        COMMON/PYINT1/MINT(400),VINT(400)
        COMMON/PYJETS/N,NPAD,K(KSZJ,5),P(KSZJ,5),V(KSZJ,5)
        common/sa1/kjp21,non1,bp,iii,neve,nout,nosc
        common/sa10/csnn,cspin,cskn,cspipi,cspsn,cspsm,rcsit,ifram, &
         iabsb,iabsm,i_sigma_AQM,ajpsi,csspn,csspm,csen
        common/sa18/i_deex,n_deex_step,i_pT_coal,i_pT_endpoint,a_FF,aPS_c,aPS_b
        common/sa24/adj1(40),nnstop,non24,zstop
        common/sa26/ndiq(kszj),npt(kszj),ifcom(kszj),idi,idio
        common/sa27/itime,kjp22,gtime,astr,akapa(6),parj1,parj2,parj3, &
         parj21,parj4,adiv,gpmax,nnc
        common/sa28/nstr,nstra(kszj),nstrv(kszj),nstr0, &
         nstr1,nstr1a(kszj),nstr1v(kszj)
        common/sa30/vneump,vneumt,mstptj
        common/sa33/smadel,ecce,secce,parecc,iparres
        common/sa34/itorw,iikk,cp0,cr0,kkii
        common/sbe/nbe,non_be,kbe(kszj,5),pbe(kszj,5),vbe(kszj,5)
        common/saf/naf,nonaf,kaf(kszj,5),paf(kszj,5),vaf(kszj,5)
        common/sbh/nbh,nonbh,kbh(kszj,5),pbh(kszj,5),vbh(kszj,5)
        common/syspar/ipden,itden,suppm,suptm,suppc,suptc,r0p,r0t, &
         nap,nat,nzp,nzt,pio
        common/sgam/ngam,nongam,kgam(kszj,5),pgam(kszj,5),vgam(kszj,5)
!       For the parini related statistics.
        common/anly_parini1/ woun, swoun, snpctl0, snpctlm, snpar
!       For the Angantyr related statistics.
        common/anly_angantyr2/ sum_bParam_ANG(3), sum_sigma_ANG(8,2), &
                               sum_Ncoll_ANG(8), sum_Npart_ANG(4,2)
!       For the gamma related statistics.
        common/anly_gamma1/ egam,  egam1,  egam2,  egam3, &
                           segam, segam1, segam2, segam3
!       For the charge and 4-momentum check.
        common/cp_check/ p_init(4), p_final(4)
!       For the simulation control.
        COMMON/SA1_PY8/ i_mode, i_tune, KF_woDecay(100), &
               KF_proj, KF_targ, win, energy_B, psno, b_min, b_max
!       Arrays of particle information of one NN pair for PYTHIA 8.
        COMMON/PYJETS_PY8/ N_PY8, NPAD_PY8, &
               K_PY8(KSZJ_PY8,8), P_PY8(KSZJ_PY8,7), V_PY8(KSZJ_PY8,5)
        DATA i_call /0/
        dimension ps0(6), ps1(6)


        time = 0D0
!       Hadronization model.
        i_had_model = INT( adj1(12) )
!       Backup.
        MSTJ21  = MSTJ(21)
        MSTP111 = MSTP(111)

!       Prints initialization and differential cross section maximum only once.
        i_call = i_call + 1
        MSTP(122) = 0
        if( i_call == 1 )then
            MSTP(122)=1
        endif

!       To accelerate the calculation of PYTHIA 8, the forward initialization
!        has been done by calling the "SUBROUTINE PAINIT_KF" indisde
!        "initialize_system" in main (only once).

        i_PAEVNT = 0
99999   continue
        i_PAEVNT  = i_PAEVNT + 1
        MSTP(111) = mstptj
!       For the debug.
        if( i_had_model == 1 .AND. i_deex > 99 )then
            MSTJ(21)  = 0
            MSTP(111) = 1
        end if
        VINT(139) = bp
        MINT(5)   = iii


!-------------------------------------------------------------------------------
!--------------------------   PYTHIA 6 Initializing   --------------------------
!       Initialization for PYTHIA 6.
        iikk = 0
        kkii = 0
        energy_A = win
        if( .NOT.IS_PYTHIA8(i_mode) )then
            ifram_in = ifram
            N = 2
            do i=1,5,1
                P(1,i) = 0D0
                P(2,i) = 0D0
            end do
            ! Fixed target (FXT) with the incident momentum +pz.
            if( ifram == 0 )then
                ifram_in = 2
            ! CMS.
            else if( ifram == 1 )then
                ! ifram_in = 1
            ! Back-to-back collisions with eA +Z and eB -Z.
            else if( ifram == 2 )then
                ifram_in = 3
                if( energy_A < PYMASS( KF_proj ) )then
                    if( iii == 1 )then
                        write(22,*)
                        write(22,*) "Warning! win(energy_A) < m_proj " // &
                                    "in the back-to-back collision " // &
                                    "frame of B/C-framework. " // &
                                    "Re-evaluated the energy " // &
                                    "according to m_proj."
                        write(22,*)
                        win = PYMASS( KF_proj )
                        energy_A = win
                    end if
                end if
                if( energy_B < PYMASS( KF_proj ) )then
                    if( iii == 1 )then
                        write(22,*)
                        write(22,*) "Warning! energy_B < m_targ " // &
                                    "in the back-to-back collision " // &
                                    "frame of B/C-framework. " // &
                                    "Re-evaluated the energy " // &
                                    "according to m_targ."
                        write(22,*)
                        energy_B = PYMASS( KF_targ )
                    end if
                end if
                P(1,3) =  SQRT( energy_A**2 - PYMASS( KF_proj )**2 )
                P(2,3) = -SQRT( energy_B**2 - PYMASS( KF_targ )**2 )
                if( ABS(P(1,3)) < 1D-10 .AND. ABS(P(2,3)) < 1D-10 )then
                    write(22,*)
                    write(22,*) "Abort! Too low win(energy_A) and energy_B!"
                    write(22,*)
                    stop
                end if
            end if
            ! In Pythia8_fort_interface.f90.
            call PAINIT_KF( ifram_in, KF_proj, KF_targ, win )
!       Sums of px, py, pz, E, inv. m, and charge before the excution.
            do i=1,6,1
                ps0(i) = PAPYP(0,i)
            end do
            if( ifram == 0 )then
                ps0(3) = win
                ps0(4) = SQRT( win**2 + P(1,5)**2 ) + P(2,5)
            else if( ifram == 2 )then
                ps0(3) = SQRT( energy_A**2 - P(1,5)**2 ) &
                       - SQRT( energy_B**2 - P(2,5)**2 )
                ps0(4) = energy_A + energy_B
            end if
        end if
!--------------------------   PYTHIA 6 Initializing   --------------------------
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
!------------------------   Core Processes Generating   ------------------------
!       Executes the collision in PYTHIA 6 or 8.
        ! In Pythia8_fort_interface.f90.
        call PAEVNT
!------------------------   Core Processes Generating   ------------------------
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
!----------------------------   PYTHIA 6 Checking   ----------------------------
!       For PYTHIA 6 mode, especially when the energy < 10 GeV.
        if( .NOT.IS_PYTHIA8(i_mode) )then
!       Sums of px, py, pz, E, inv. m, and charge after the excution.
            ps1=0.
            do i=1,6,1
                ps1(i)=PAPYP(0,i)
            end do
!       Charge is not conserved. Re-generates the event.
            if( ABS(ps0(6)-ps1(6)) > 1D-10 ) goto 99999
!       4-momentum is not conserved. Re-generates the event.
            do i=1,4,1
                if( ABS(ps0(i)-ps1(i)) > 1D-5 ) goto 99999
            end do
!       Error counters in PYTHIA 6. Re-generates the event if any errors exist.
            ! if( MSTU(23) > 0 .OR. MSTU(30) > 0 ) goto 99999
!       Error type in PYTHIA. Re-generates the event if any errors exist.
            if( MSTU(24) > 0 ) goto 99999
!       Error due to the "iikk, kkii" in PYTHIA 6. Re-generates the event.
            if( iikk == 2 .OR. kkii == 2 ) goto 99999
!       Re-generates the event if there are junctions when consider inelatic
!        processes in parton rescattering via SFM. Reduces complexity.
            if( ( i_had_model == 0 .OR. i_had_model == 3) &
                .AND. ( iparres == 1 .OR. iparres == 3 ) )then
                    if( MSTU(29) == 1 ) goto 99999
            end if
        end if
        ! Recovers.
        MSTJ(21)  = MSTJ21
        MSTP(111) = MSTP111
!----------------------------   PYTHIA 6 Checking   ----------------------------
!-------------------------------------------------------------------------------


!       Records the initial 4-momentum.
        if( ifram == 0 )then
            p_init(3) = win
            p_init(4) = SQRT( win**2 + PYMASS(KF_proj)**2 ) + PYMASS(KF_targ)
        else if( ifram == 1 )then
            p_init(4) = win
        else if( ifram == 2 )then
            p_init(3) = SQRT( energy_A**2 - PYMASS( KF_proj )**2 ) &
                      - SQRT( energy_B**2 - PYMASS( KF_targ )**2 )
            p_init(4) = energy_A + energy_B
        end if

!       Records the original contents of each collision after PYTHIA execution.
        ! In Pythia8_fort_interface.f90.
        call PARECO( 1 )
!       full_events_history of OSC1999A
        N0 = N
        N = 2
        if( IS_PYTHIA8(i_mode) )then
            do i=1,5,1
                K(1,i) = K_PY8(2,i)
                K(2,i) = K_PY8(3,i)
                P(1,i) = P_PY8(2,i)
                P(2,i) = P_PY8(3,i)
                V(1,i) = V_PY8(2,i)
                V(2,i) = V_PY8(3,i)
            end do
        end if
        call oscar(0)
        N = N0


!-------------------------------------------------------------------------------
!----------------------------   Entries Filtering   ----------------------------
!       Uses V(i,5) as index pointers for the further recovering after parcas.
!       0.1 is just a number used to ensure the accurate conversion from INT().
        if( .NOT.IS_PYTHIA8(i_mode) )then
            do i=1,N,1
                if( mstptj == 0 ) V(i,5) = i*1D0 + 0.1D0
                ! Removes junctions for Coal, if any.
                if( i_had_model == 1 .AND. K(i,2) == 88 ) K(i,1) = 51
            end do
            call PAEDIT(1)
!       Moves non-historical from /PYJETS_PY8/ to /PYJETS/ for parcas.
        else if( IS_PYTHIA8(i_mode) )then
            N = 0
            do i=1,N_PY8,1
                KS = K_PY8(i,1)
                if( .NOT.IS_EXIST(KS,i_mode) ) cycle
                N = N + 1
                do j=1,5,1
                    K(N,j) = K_PY8(i,j)
                    P(N,j) = P_PY8(i,j)
                    V(N,j) = V_PY8(i,j)
                end do
                if( mstptj == 0 ) V(N,5) = i*1D0 + 0.1D0
            end do
        end if
!----------------------------   Entries Filtering   ----------------------------
!-------------------------------------------------------------------------------


!       Gives four positions to the particles generated in PYTHIA ("PYJETS"),
!        except for those from Angantyr which has been given by Angantyr.
        ! in parini.f90
        if( i_mode /= 8 .AND. i_mode /= 9 ) call ptcre(1,2,time)


!-------------------------------------------------------------------------------
!-------------------------   Gamma & Hadron Removing   -------------------------

!----------------------   B-framework Gamma 22 Removing  -----------------------
!       Removes photons from "PYJETS" to "sgam".
        if( mstptj == 1  )then
            KF_gamma_decay = 22
            call remo_gam( KF_gamma_decay )
!       egam3: total gamma energy after hadronization
            call prt_sgam( ngam, egam3, 3 )
            goto 333
        end if
!----------------------   B-framework Gamma 22 Removing  -----------------------

!----------------------   C-framework Gamma 44 Removing  -----------------------
!       Removes gamma from "PYJETS" to "sgam".
!       Do not remove any entries for CR if C-framework via sfm.
        if( i_had_model /= 0 .OR. ( i_had_model == 0 .AND. kjp22 < 5 ) )then
            n44 = 0
            do j=1,N,1
                kf = k(j,2)
                if(kf == 22)then
                    k(j,2) = 44
                    n44 = n44 + 1
                end if
            end do
!           Moves "44" from "PYJETS" to "sgam".
            if(n44 > 0) call remo_gam(44)
!           egam1: energy of gamma after partonic initiation
            call prt_sgam( n44, egam1, 1 )
!----------------------   C-framework Gamma 44 Removing  -----------------------

!-----------------------------   debug   ------------------------------
            if( i_had_model == 1 .AND. i_deex > 99 )then
                ! In coales.f90
                call do_debug( i_deex )
!-----------------------------   debug   ------------------------------

!-----------------------   C-framework Hadron Removing   -----------------------
!           Removes hadrons/leptons from "PYJETS" to "sbh".
            else
            ! In parini.f90
                call remo
            end if
        end if
!-----------------------   C-framework Hadron Removing   -----------------------

!-------------------------   Gamma & Hadron Removing   -------------------------
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
!----------------------------   Diquark Locating   -----------------------------
!       Prepares for the parton rescattering.
        if( i_had_model == 0 .OR. i_had_model == 3 )then
            i3 = 0
            do i1=1,N,1
                i3 = i1
                KS = K(i1,1)
                KF = k(i1,2)
                ! Identifies diquarks.
                ndiq(i1) = 0
                if( IS_DIQUARK(KF) .AND. IS_EXIST(KS,i_mode) )then
                    idi = idi + 1
                    ndiq(i1) = idi
                end if
                ! 'PYJETS' to 'sbe'
                do i2=1,5
                    kbe(i3,i2) = K(i1,i2)
                    pbe(i3,i2) = P(i1,i2)
                    vbe(i3,i2) = V(i1,i2)
                end do
            end do
            nbe = i3
        end if
!----------------------------   Diquark Locating   -----------------------------
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
!---------------------------   Diquark Breaking-up   ---------------------------
!       Breaks up diquark and give four momentum and four position
!        to broken quarks (working in 'PYJETS', in parini.f90)
        if( iparres == 0 .OR. iparres == 1 )  call break
!---------------------------   Diquark Breaking-up   ---------------------------
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
!-----------------------------   String Locating   -----------------------------
!       Finds number of strings and line number of first and last components
!        of each string.
!#TODO(Lei20240218): not work well for PYTHIA 8, need improvement.
        nstr1 = 0
        if( .NOT.IS_PYTHIA8(i_mode) .AND. i_had_model /= 1 )then
            jb = 0
10000       continue
            do i1 = jb+1, N, 1
                ! i1 is 'A'
                if( K(i1,1) == 2 )then
                    do i2 = i1+1, N, 1
                        ! i2 is 'V'
                        if( K(i2,1) == 1 )then
                            nstr1 = nstr1 + 1
                            ! line number of first component of nstr1-th string
                            nstr1a(nstr1) = i1
                            ! line number of first component of nstr1-th string
                            nstr1v(nstr1) = i2
                            jb = i2
                            goto 10000
                        end if
                    end do
                end if
            end do
        end if
!       nstr1: number of strings after call break
        nstr0 = nstr1
        idio  = idi
!-----------------------------   String Locating   -----------------------------
!-------------------------------------------------------------------------------


333     continue


!-------------------------------------------------------------------------------
!-------------------------   Information Statistics   --------------------------
!       Counts single-event cross sections.
        call PASTAT(-1,1)

!---------------------------   Angantyr Statistics   ---------------------------
        if( i_mode == 8 .OR. i_mode == 9 )then
            ! The impact parameter.
            bp   = VINT(139)
            ! The impact parameter angle phi.
            bPhi = VINT(370)
            ! Avoids potential divergence from PYTHIA8/Angantyr pA/AA + HardQCD.
            ! We cannot obtain accurate results from psno=3 but psno /= 3 does.
            if( INT(psno) /= 3 )then
                VINT(98)  = iii*1D0
                VINT(99)  = 1D0
                VINT(100) = 1D0
                PARI(2)   = PARI(1) / VINT(98)
                PARI(9)   = 1D0
                PARI(10)  = 1D0
            end if
            weightSum = VINT(98)
            ! one_over_weight = VINT(99)   ! 1/weight
            bWeight   = VINT(100)
            ! Cross sections.
            glauberTot    = VINT(371)
            glauberND     = VINT(372)
            glauberINEL   = VINT(373)
            glauberEL     = VINT(374)
            glauberDiffP  = VINT(375)
            glauberDiffT  = VINT(376)
            glauberDDiff  = VINT(377)
            glauberBSlope = VINT(378)
            ! Error.
            glauberTotErr    = VINT(381)
            glauberNDErr     = VINT(382)
            glauberINELErr   = VINT(383)
            glauberELErr     = VINT(384)
            glauberDiffPErr  = VINT(385)
            glauberDiffTErr  = VINT(386)
            glauberDDiffErr  = VINT(387)
            glauberBSlopeErr = VINT(388)

            nCollTot    = MINT(385)
            nCollND     = MINT(386)
            nCollNDTot  = MINT(387)
            nCollSDP    = MINT(388)
            nCollSDT    = MINT(389)
            nCollDD     = MINT(390)
            nCollCD     = MINT(391)
            nCollEL     = MINT(392)

            nPartProj   = MINT(393)
            nAbsProj    = MINT(394)
            nDiffProj   = MINT(395)
            nElProj     = MINT(396)
            nPartTarg   = MINT(397)
            nAbsTarg    = MINT(398)
            nDiffTarg   = MINT(399)
            nElTarg     = MINT(400)
!           Event accumulation.
            sum_bParam_ANG(1) = sum_bParam_ANG(1) + bp   * bWeight
            sum_bParam_ANG(2) = sum_bParam_ANG(2) + bPhi * bWeight
            sum_bParam_ANG(3) = sum_bParam_ANG(3) + bWeight**2
            sum_sigma_ANG(1,1)  = sum_sigma_ANG(1,1)  + glauberTot      *bWeight
            sum_sigma_ANG(1,2)  = sum_sigma_ANG(1,2)  + glauberTotErr   *bWeight
            sum_sigma_ANG(2,1)  = sum_sigma_ANG(2,1)  + glauberND       *bWeight
            sum_sigma_ANG(2,2)  = sum_sigma_ANG(2,2)  + glauberNDErr    *bWeight
            sum_sigma_ANG(3,1)  = sum_sigma_ANG(3,1)  + glauberINEL     *bWeight
            sum_sigma_ANG(3,2)  = sum_sigma_ANG(3,2)  + glauberINELErr  *bWeight
            sum_sigma_ANG(4,1)  = sum_sigma_ANG(4,1)  + glauberEL       *bWeight
            sum_sigma_ANG(4,2)  = sum_sigma_ANG(4,2)  + glauberELErr    *bWeight
            sum_sigma_ANG(5,1)  = sum_sigma_ANG(5,1)  + glauberDiffP    *bWeight
            sum_sigma_ANG(5,2)  = sum_sigma_ANG(5,2)  + glauberDiffPErr *bWeight
            sum_sigma_ANG(6,1)  = sum_sigma_ANG(6,1)  + glauberDiffT    *bWeight
            sum_sigma_ANG(6,2)  = sum_sigma_ANG(6,2)  + glauberDiffTErr *bWeight
            sum_sigma_ANG(7,1)  = sum_sigma_ANG(7,1)  + glauberDDiff    *bWeight
            sum_sigma_ANG(7,2)  = sum_sigma_ANG(7,2)  + glauberDDiffErr *bWeight
            sum_sigma_ANG(8,1)  = sum_sigma_ANG(8,1)  + glauberBSlope   *bWeight
            sum_sigma_ANG(8,2)  = sum_sigma_ANG(8,2)  + glauberBSlopeErr*bWeight
            do i=1,8,1
                j = i + 384
                sum_Ncoll_ANG(i) = sum_Ncoll_ANG(i) + MINT(j) * bWeight
            end do
            sum_Npart_ANG(1,1) = sum_Npart_ANG(1,1) + nPartProj *bWeight
            sum_Npart_ANG(2,1) = sum_Npart_ANG(2,1) + nAbsProj  *bWeight
            sum_Npart_ANG(3,1) = sum_Npart_ANG(3,1) + nDiffProj *bWeight
            sum_Npart_ANG(4,1) = sum_Npart_ANG(4,1) + nElProj   *bWeight
            sum_Npart_ANG(1,2) = sum_Npart_ANG(1,2) + nPartTarg *bWeight
            sum_Npart_ANG(2,2) = sum_Npart_ANG(2,2) + nAbsTarg  *bWeight
            sum_Npart_ANG(3,2) = sum_Npart_ANG(3,2) + nDiffTarg *bWeight
            sum_Npart_ANG(4,2) = sum_Npart_ANG(4,2) + nElTarg   *bWeight

            ! The number of wounded and unwounded nucleons.
            woun = nAbsProj + nDiffProj + nAbsTarg + nDiffTarg
            unwoun = nap + nat - woun
        end if
!---------------------------   Angantyr Statistics   ---------------------------

!-------------------------   Information Statistics   --------------------------
!-------------------------------------------------------------------------------


        return
        end



!ccccccccccccccccccccccccccccccccccccc end ccccccccccccccccccccccccccccccccccccc