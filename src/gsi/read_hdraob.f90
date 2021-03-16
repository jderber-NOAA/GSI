subroutine read_hdraob(nread,ndata,nodata,infile,obstype,lunout,twindin,sis,&
     prsl_full,hgtl_full,nobs,nrec_start)
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:  read_hdraob                 read obs from hdraob file
!   prgmmr: derber          org: np22                date: 2020-02-12
!
! abstract:  This routine reads high resolution raob data found in the hdraob
!            file.  Specific observation types read by this routine 
!            include surface pressure, temperature, winds (components
!            and speeds), and moisture.  
!
!            When running the gsi in regional mode, the code only
!            retains those observations that fall within the regional
!            domain
!
! program history log:
!   2020-02-24  derber

!   input argument list:
!     infile   - unit from which to read BUFR data
!     obstype  - observation type to process
!     lunout   - unit to which to write data for further processing
!     prsl_full- 3d pressure on full domain grid
!     hgtl_full- 3d height on full domain grid
!     nrec_start - number of subsets without useful information
!
!   output argument list:
!     nread    - number of type "obstype" observations read
!     nodata   - number of individual "obstype" observations read
!     ndata    - number of type "obstype" observations retained for further processing
!     twindin  - input group time window (hours)
!     sis      - satellite/instrument/sensor indicator
!     nobs     - array of observations on each subdomain for each processor
!
! attributes:
!   language: f90
!   machine:  ibm RS/6000 SP
!
!$$$
  use kinds, only: r_single,r_kind,r_double,i_kind
  use constants, only: zero,one_tenth,one,deg2rad,fv,t0c,half,&
      three,four,rad2deg,tiny_r_kind,huge_r_kind,huge_i_kind,&
      r60inv,r10,r100,r2000,eps,omeps
  use constants,only: rearth,stndrd_atmos_ps,rd,grav
  use gridmod, only: diagnostic_reg,regional,nlon,nlat,nsig,&
      tll2xy,txy2ll,rotate_wind_ll2xy,rotate_wind_xy2ll,&
      rlats,rlons,twodvar_regional,fv3_regional
  use convinfo, only: nconvtype,ctwind, &
      ncmiter,ncgroup,ncnumgrp,icuse,ictype,icsubtype,ioctype, &
      use_prepb_satwnd

  use obsmod, only: iadate,oberrflg,perturb_obs,perturb_fact,ran01dom
  use obsmod, only: blacklst,offtime_data,bmiss,ext_sonde
  use converr,only: etabl
  use converr_ps,only: etabl_ps,isuble_ps,maxsub_ps
  use converr_q,only: etabl_q,isuble_q,maxsub_q
  use converr_t,only: etabl_t,isuble_t,maxsub_t
  use converr_uv,only: etabl_uv,isuble_uv,maxsub_uv
  use convb_ps,only: btabl_ps
  use convb_q,only: btabl_q
  use convb_t,only: btabl_t
  use convb_uv,only: btabl_uv
  use gsi_4dvar, only: l4dvar,l4densvar,time_4dvar,winlen
  use qcmod, only: errormod,errormod_hdraob,noiqc,newvad,njqc
  use qcmod, only: nrand
  use mpimod, only: mype
  use nltransf, only: nltransf_forward
  use blacklist, only : blacklist_read,blacklist_destroy
  use blacklist, only : blkstns,blkkx,ibcnt
  use ndfdgrids,only: adjust_error
  use jfunc, only: tsensible
  use deter_sfc_mod, only: deter_sfc_type,deter_sfc2
  use mpimod, only: npe
  use gsi_io, only: verbose

  implicit none

! Declare passed variables
  character(len=*)                      ,intent(in   ) :: infile,obstype
  character(len=20)                     ,intent(in   ) :: sis
  integer(i_kind)                       ,intent(in   ) :: lunout,nrec_start
  integer(i_kind)                       ,intent(inout) :: nread,ndata,nodata
  integer(i_kind),dimension(npe)        ,intent(inout) :: nobs
  real(r_kind)                          ,intent(in   ) :: twindin
  real(r_kind),dimension(nlat,nlon,nsig),intent(in   ) :: prsl_full,hgtl_full

! Declare local parameters
  real(r_kind),parameter:: r0_01 = 0.01_r_kind
  real(r_kind),parameter:: r0_75 = 0.75_r_kind
  real(r_kind),parameter:: r0_7 = 0.7_r_kind
  real(r_kind),parameter:: r1_2 = 1.2_r_kind
  real(r_kind),parameter:: r1_02 = 1.02_r_kind
  real(r_kind),parameter:: r3_33= three + one/three
  real(r_kind),parameter:: r6   = 6.0_r_kind
  real(r_kind),parameter:: r20  = 20.0_r_kind
  real(r_kind),parameter:: r50  = 50.0_r_kind
  real(r_kind),parameter:: r90  = 90.0_r_kind
  real(r_kind),parameter:: r360 = 360.0_r_kind
  real(r_kind),parameter:: r500 = 500.0_r_kind
  real(r_kind),parameter:: r999 = 999.0_r_kind
  real(r_kind),parameter:: r1200= 1200.0_r_kind
  real(r_kind),parameter:: convert= 1.0e-6_r_kind
  real(r_kind),parameter:: emerr= 0.2_r_kind
  real(r_kind),parameter:: r0_1_bmiss=one_tenth*bmiss
  real(r_kind),parameter:: r0_01_bmiss=r0_01*bmiss

!  integer(i_kind),parameter:: mxtb=5000000
!  integer(i_kind),parameter:: nmsgmax=100000 ! max message count
  integer(i_kind), parameter:: maxlevs=10000

! Declare local variables
  logical tob,qob,uvob,psob
  logical outside,convob
  logical sfctype
  logical luse,windcorr
  logical patch_fog
  logical acft_profl_file
  logical,allocatable,dimension(:,:):: lmsg           ! set true when convinfo entry id found in a message

  character(80) hdstr,hdstr2,hdstr3,levstr
  character(40) maxtmintstr,owavestr
  character(80) obstr
  character(10) date
  character(8) subset
  character(8) prvstr,sprvstr     
  character(8) c_station_id
  character(8) cc_station_id
  character(1) sidchr(8)
  character(8) stnid
  character(1) cb
  character(1) cdummy

  integer(i_kind) ireadmg,ireadsb,icntpnt,icntpnt2,icount,iiout
  integer(i_kind) lunin,i,maxobs,j,idomsfc,it29,nmsgmax,mxtb
  integer(i_kind) kk,klon1,klat1,klonp1,klatp1
  integer(i_kind) nc,nx,isflg,ntread,itx,ii,ncsave
  integer(i_kind) ihh,idd,idate,iret,im,iy,k,levs
  integer(i_kind) kx,kx0,nreal,nchanl,ilat,ilon
  integer(i_kind) lim_qm
  integer(i_kind) ntmp,iout
  integer(i_kind) pflag,irec,zflag
  integer(i_kind) ntest,nvtest,iosub,ixsub,isubsub
  integer(i_kind) kl,k1,k2,k1_ps,k1_q,k1_t,k1_uv,k2_q,k2_t,k2_uv,k2_ps
  integer(i_kind) itypex,itypey,id
  integer(i_kind) minobs,minan
  integer(i_kind) ntb,ntmatch,ncx,cat
  integer(i_kind) nmsg                ! message index
  integer(i_kind) jj,start,next,ncount_ps,ncount_q,ncount_uv,ncount_t
  integer(i_kind),dimension(5):: idate5
  integer(i_kind),dimension(nconvtype)::ntxall
  integer(i_kind),dimension(nconvtype+1)::ntx
  integer(i_kind),allocatable,dimension(:):: isort,iloc,nrep
  integer(i_kind),allocatable,dimension(:,:):: tab
  integer(i_kind) ibfms,thisobtype_usage
  integer(i_kind) iwmo,ios
  integer(i_kind) ntime,itime,igroup,istation
  integer(i_kind) ierr_ps,ierr_q,ierr_t,ierr_uv   !  the position of error table collum
  integer(i_kind),dimension(maxlevs):: pqm,qqm,tqm,wqm
  integer(i_kind) idummy1,idummy2,glret,lindx !glret>0 means GLERL code exists.Others are dummy variables
  real(r_kind) time,timex,timeobs,toff,t4dv,zeps
  real(r_kind) qtflg,tdry,ediff,usage,ediff_ps,ediff_q,ediff_t,ediff_uv
  real(r_kind) u0,v0,uob,vob,dx,dy,dx1,dy1,w00,w10,w01,w11
  real(r_kind) qoe,qobcon,dlnpob,ppb,poe,obheight
  real(r_kind) toe,woe,errout,oelev,dlat,dlon,dlat_earth,dlon_earth
  real(r_kind) dlat_earth_deg,dlon_earth_deg
  real(r_kind) elev,stnelev
  real(r_kind) cdist,disterr,disterrmax,rlon00,rlat00
  real(r_kind) vdisterrmax,u00,v00
  real(r_kind) rminobs,rminan,es,dummy
  real(r_kind) del,terrmin,werrmin,perrmin,qerrmin,del_ps,del_q,del_t,del_uv
  real(r_kind) pjbmin,qjbmin,tjbmin,wjbmin
  real(r_kind) tsavg,ff10,sfcr,zz
  real(r_kind) time_correction,fact
  real(r_kind),dimension(nsig):: presl,hgtl,lnpresl
  real(r_kind),dimension(nsig-1):: dpres
  real(r_kind),dimension(maxlevs)::plevs
  real(r_kind),allocatable,dimension(:,:):: cdata_all,cdata_out
  real(r_kind) :: zob,tref,dtw,dtc,tz_tr

  real(r_double) rstation_id
  real(r_double) vtcd,glcd !virtual temp program code and GLERL program code
  real(r_double),dimension(2):: hdr
  real(r_double),dimension(8):: hdr2
  real(r_double),dimension(1):: hdr3
  real(r_double),dimension(2,maxlevs):: levdat
  real(r_double),dimension(8,maxlevs):: var_jb,obserr
  real(r_double),dimension(8,maxlevs):: obsdat



!  equivalence to handle character names
  equivalence(rstation_id,c_station_id)
  equivalence(rstation_id,sidchr)

!  data statements
  data hdstr  /'WMOB WMOS'/
  data hdstr2 /'YEAR MNTH DAYS HOUR MINU SECO CLATH CLONH' /
  data hdstr3 /'HEIT'/
  data obstr  /'LTDS LATDH LONDH GP10 WSPD WDIR TMDB TMDP'/
  data levstr  /'PRLC GP07'/
  data lunin / 13 /
  !* for match loction station and time
!       character(7*2000) cstn_idtime,cstn_idtime2
!       character(7) stn_idtime(2000),stn_idtime2(2000)
!       equivalence (stn_idtime(1),cstn_idtime)
!       equivalence (stn_idtime2(1),cstn_idtime2)
!       integer :: ii1,atmp,btmp,mytimeyy,ibyte
!       character(4) stid
!       real(8) :: rval
!       character(len=8) :: cval
!       equivalence (rval,cval)
!       character(7) flnm

  logical print_verbose
  
! print_verbose=.false.
  print_verbose=.true.
  if(verbose) print_verbose=.true.

! Initialize variables

  tob = obstype == 't'
  uvob = obstype == 'uv' 
  qob = obstype == 'q'
  psob = obstype == 'ps'

  zflag=0
  nreal=0
  if(tob)then
     nreal=25
  else if(uvob) then 
     nreal=26
  else if(psob) then
     nreal=20
  else if(qob) then
     nreal=26
  else 
     write(6,*) ' illegal obs type in READ_HDROAB ',obstype
     call stop2(94)
  end if

  if(perturb_obs .and. (tob .or. psob .or. qob))nreal=nreal+1
  if(perturb_obs .and. uvob )nreal=nreal+2


  if (blacklst) call blacklist_read(obstype)

  lim_qm=1
  pqm=0
  tqm=0
  qqm=0
  wqm=0
  terrmin=half
  werrmin=one
  perrmin=0.3_r_kind
  qerrmin=0.05_r_kind
  tjbmin=zero
  qjbmin=zero
  wjbmin=zero
  pjbmin=zero
!------------------------------------------------------------------------
  ntread=1
  ntmatch=0
  ntx(ntread)=0
  ntxall=0
  var_jb=zero


  do nc=1,nconvtype
     if(trim(ioctype(nc)) == trim(obstype))then
          ntmatch=ntmatch+1
          ntxall(ntmatch)=nc
     end if
  end do

  irec = 0
!! get message and subset counts

  call getcount_bufr(infile,nmsgmax,mxtb)
  allocate(lmsg(nmsgmax,ntread),tab(mxtb,3),nrep(nmsgmax))

  lmsg = .false.
  maxobs=0
  tab=0
  nmsg=0
  nrep=0
  ntb = 0
  ncount_ps=0;ncount_q=0;ncount_t=0;ncount_uv=0

! Open, then read date from bufr data
! call closbf(lunin)
  open(lunin,file=trim(infile),form='unformatted')
  call openbf(lunin,'IN',lunin)
  call datelen(10)

!  This is the message creating error message
  msg_report: do while (ireadmg(lunin,subset,idate) == 0)
     irec = irec + 1
     if(irec < nrec_start) cycle msg_report

!    Time offset
     if(nmsg == 0) call time_4dvar(idate,toff)
     nmsg=nmsg+1
     if (nmsg>nmsgmax) then
        write(6,*)'READ_HDRAOB: messages exceed maximum ',nmsgmax
        call stop2(50)
     endif
     loop_report: do while (ireadsb(lunin) == 0)
        ntb = ntb+1
        if (ntb>mxtb) then
           write(6,*)'READ_HDRAOB: reports exceed maximum ',mxtb
           call stop2(50)
        endif

!       Extract type information
        if(psob .or. tob .or. qob)kx=119
        if(uvob)kx=219


!       call ufbint(lunin,hdr,2,1,iret,hdstr)
!       igroup=hdr(1)
!       istation=hdr(2)
!       write(6,*) igroup,istation
!       call ufbint(lunin,obsdat,13,maxlevs,levs,obstr)
!       Check for blacklisting of station ID
!       if (blacklst .and. ibcnt > 0) then
!          stnid = transfer(hdr(4),stnid)
!          do i = 1,ibcnt
!             if( kx == blkkx(i) .and. stnid == blkstns(i) ) then
!                write(6,*)'READ_HDRAOB: blacklist station ',stnid, &
!                   'for obstype ',trim(obstype),' and kx=',kx
!                cycle loop_report
!             endif
!          enddo
!       endif

        call ufbint(lunin,levdat,2,maxlevs,levs,levstr)
        if(uvob .and. levdat(1,2) > 1.e8)kx=218
!  Match ob to proper convinfo type
        ncsave=0
        matchloop:do ncx=1,ntmatch
           nc=ntxall(ncx)
           if (kx == ictype(nc))then 

              ncsave=nc
              exit matchloop
           end if

        end do matchloop

!       write(6,*) ' levs =',levs
        maxobs=maxobs+max(1,levs)
        nx=1
        tab(ntb,1)=ncsave
        tab(ntb,2)=nx
        tab(ntb,3)=levs

     end do loop_report
  enddo msg_report
  if (nmsg==0) then
     call closbf(lunin)
     close(lunin)
     if(print_verbose)write(6,*)'READ_HDRAOB: no messages/reports '
     return
  end if
  if(print_verbose)write(6,*)'READ_HDRAOB: messages/reports = ',nmsg,'/',ntb,' ntread = ',ntread


  if(qob .and. print_verbose) write(6,*)'READ_HDRAOB: time offset is ',toff,' hours.'
!------------------------------------------------------------------------


! loop over convinfo file entries; operate on matches
  
  allocate(cdata_all(nreal,maxobs),isort(maxobs))
  isort = 0
  cdata_all=zero
  nread=0
  ntest=0
  nvtest=0
  nchanl=0
  ilon=2
  ilat=3
  pflag=0
  loop_convinfo: do nx=1,ntread

     call closbf(lunin)
     open(lunin,file=infile,form='unformatted')
     call openbf(lunin,'IN',lunin)
     call datelen(10)

!    Big loop over hdraob file	

     ntb = 0
     nmsg = 0
     icntpnt=0
     icntpnt2=0
     disterrmax=-9999.0_r_kind
     irec = 0
     loop_msg: do while (ireadmg(lunin,subset,idate)== 0)
        irec = irec + 1
        if(irec < nrec_start) cycle loop_msg

        loop_readsb: do while(ireadsb(lunin) == 0)
!          use msg lookup table to decide which messages to skip
!          use report id lookup table to only process matching reports
           ntb = ntb+1
           nc=tab(ntb,1)
           kx = ictype(nc) 
           if(nc <= 0 .or. tab(ntb,2) /= nx) then
              write(6,*) nc,tab(ntb,2),ntb
              cycle loop_readsb
           end if
                 
!          Extract type, date, and location information
           call ufbint(lunin,hdr,2,1,iret,hdstr)
           call ufbint(lunin,hdr2,8,1,iret,hdstr2)
           call ufbint(lunin,hdr3,1,1,iret,hdstr3)
           igroup=hdr(1)
           istation=hdr(2)
           id=1000*igroup+istation
           call ufbint(lunin,levdat,2,maxlevs,levs,levstr)
           if(psob)write(6,*) igroup,istation,id,levs,levdat(1,1)
           if(igroup < 0 .or. istation < 0 .or. id >= 100000)then
               write(6,*) ' hdr ',hdr
               cycle loop_readsb
           end if 
           write(c_station_id,'(i5,3x)') id
           if(levdat(1,2) > 1.e8 .and. (obstype == 'q' .or. obstype == 't'))cycle loop_readsb
           call ufbint(lunin,obsdat,8,maxlevs,levs,obstr)

!     Combine height obs into single array and divide by g
           do k=1,levs
              levdat(1,k)=.01*levdat(1,k)
              if(levdat(2,k) > 1.e8)then
                if(obsdat(4,k) < 1.e8)levdat(2,k)=obsdat(4,k)/grav
              else
                levdat(2,k)=levdat(2,k)/grav
              end if
           end do
!          write(6,*) '1',(obsdat1(j),j=1,8),levs
!          do i=1,levs
!             write(6,*) i,'pres',levdat(1,i),levdat(2,i),hdr(10)
!             write(6,*) '0',(obsdat(j,i),j=1,8)
!             write(6,*) '2',(hdr3(j),j=1,8)
!          end do

!------------------------------------------------------------------------

           dlat_earth_deg=hdr2(7)
           dlon_earth_deg=hdr2(8)
           if(abs(dlat_earth_deg)>r90 ) then
              write(6,*) ' invalid lat ',id,dlat_earth_deg
                cycle loop_readsb
           end if
           if(dlon_earth_deg >= r360)dlon_earth_deg=dlon_earth_deg-r360
           if(dlon_earth_deg < zero)dlon_earth_deg=dlon_earth_deg+r360
           dlon_earth=dlon_earth_deg*deg2rad
           dlat_earth=dlat_earth_deg*deg2rad
           if(regional)then
              call tll2xy(dlon_earth,dlat_earth,dlon,dlat,outside)    ! convert to rotated coordinate
              if(outside) cycle loop_readsb   ! check to see if outside regional domain
           else
              dlat = dlat_earth
              dlon = dlon_earth
              call grdcrd1(dlat,rlats,nlat,1)
              call grdcrd1(dlon,rlons,nlon,1)
           endif
!          Interpolate guess profile to observation location
           klon1= int(dlon);  klat1= int(dlat)
           dx   = dlon-klon1; dy   = dlat-klat1
           dx1  = one-dx;     dy1  = one-dy
           w00=dx1*dy1; w10=dx1*dy; w01=dx*dy1; w11=dx*dy

           klat1=min(max(1,klat1),nlat); klon1=min(max(0,klon1),nlon)
           if (klon1==0) klon1=nlon
           klatp1=min(nlat,klat1+1); klonp1=klon1+1
           if (klonp1==nlon+1) klonp1=1
           do kk=1,nsig
              presl(kk)=w00*prsl_full(klat1 ,klon1 ,kk) +  &
                        w10*prsl_full(klatp1,klon1 ,kk) + &
                        w01*prsl_full(klat1 ,klonp1,kk) + &
                        w11*prsl_full(klatp1,klonp1,kk)
           end do

           if(levdat(1,2) < 1.e8)then
              do k=1,levs
                 plevs(k)=0.1*levdat(1,k)   ! convert mb to cb
              end do
           else
!  If pressure is missing get an approximate pressure from guess field.  This
!  pressure is used for estimating ob error only.
!  Note that balloon drift is not included in this calculation.
              do kk=1,nsig
                 hgtl(kk) =w00*hgtl_full(klat1 ,klon1 ,kk) +  &
                           w10*hgtl_full(klatp1,klon1 ,kk) + &
                           w01*hgtl_full(klat1 ,klonp1,kk) + &
                           w11*hgtl_full(klatp1,klonp1,kk)
              end do
              do k=1,levs
                 obheight=levdat(2,k)
                 if(obheight < hgtl(1))then 
                    plevs(k)=presl(1)
                 else if (obheight > hgtl(nsig))then
!  Reject ob by making pressure zero
                    plevs(k) = 0.0
                 else
                    do kk=1,nsig-1
                       if(hgtl(kk) <= obheight .and. hgtl(kk+1) > obheight)then
                          fact=(hgtl(kk+1)-obheight)/(hgtl(kk+1)-hgtl(kk))
                          plevs(k)=exp(fact*log(presl(kk))+(1.-fact)*log(presl(kk+1)))
                       end if
                    end do
                 end if
              end do
              if(uvob .and. kx /= 218)write(6,*) ' inconsistent kx ',kx
           end if
!             Compute depth of guess pressure layersat observation location
           if (.not.twodvar_regional .and. levs > 1) then
             do kk=1,nsig-1
                dpres(kk)=presl(kk)-presl(kk+1)
              end do
           endif
           idate5(1)=iadate(1)
           idate5(2)=iadate(2)
           idate5(3)=iadate(3)
           idate5(4)=iadate(4)
           idate5(5)=0
           call w3fs21(idate5,minan)    !  analysis ref time in minutes relative to historic date
           rminan=minan

           if(offtime_data) then

!             in time correction for observations to account for analysis
!                      time being different from obs file time.
              write(date,'( i10)') idate
              read (date,'(i4,3i2)') iy,im,idd,ihh
              idate5(1)=iy
              idate5(2)=im
              idate5(3)=idd
              idate5(4)=ihh
              idate5(5)=0
              call w3fs21(idate5,minobs)    !  obs ref time in minutes relative to historic date
!             Add obs reference time, then subtract analysis time to get obs
!             time relative to analysis

              time_correction=float(minobs-minan)*r60inv

           else
              time_correction=zero
           end if

 
!          get observation launch time relative to analysis time
           idate5(1)=hdr2(1)
           idate5(2)=hdr2(2)
           idate5(3)=hdr2(3)
           idate5(4)=hdr2(4)
           idate5(5)=hdr2(5)
           call w3fs21(idate5,minobs)    !  obs launch time in minutes relative to historic date
           rminobs=minobs+hdr2(6)*r60inv
 
           timeobs=(rminobs-rminan)*r60inv

           t4dv=timeobs + toff
           zeps=1.0e-8_r_kind
           if (t4dv<zero  .and.t4dv>      -zeps) t4dv=zero
           if (t4dv>winlen.and.t4dv<winlen+zeps) t4dv=winlen
           t4dv=t4dv + time_correction
           time=timeobs + time_correction

!          write(6,*) 'time ',obsdat,oberrflg
!          write(6,*) 'time ',rminobs,rminan,timeobs,time_correction,toff,levs

     
!          Extract data information on levels

!          If available, get obs errors from error table
           

!          Set lower limits for observation errors
           terrmin=half
           werrmin=one
           perrmin=0.3_r_kind
           qerrmin=0.05_r_kind
           tjbmin=zero
           qjbmin=zero
           wjbmin=zero
           pjbmin=zero
           itypey=kx
!          write(6,*) 'njqc ',njqc,maxsub_t
           if( njqc) then
              if (psob)  then
                 itypex=itypey
                 ierr_ps=0
                 do i =1,maxsub_ps
                    if(icsubtype(nc)==isuble_ps(itypex,i)) then
                       ierr_ps=i+1
                       exit
                    else if(i== maxsub_ps .and. icsubtype(nc) /= isuble_ps(itypex,i)) then
                       ncount_ps=ncount_ps+1
                       do j=1,maxsub_ps
                          if(isuble_ps(itypex,j) ==0 ) then
                             ierr_ps=j+1
                             exit
                          endif
                       enddo
                       if (ncount_ps ==1) then
                          write(6,*) 'READ_HDRAOB: WARNING!!psob: cannot find subtyep in the &
                                     error table,itype,iosub=',itypex,icsubtype(nc)
                          write(6,*) 'read error table at colomn subtype as 0, error table column=',ierr_ps
                       endif
                    endif
                 enddo
                 do k=1,levs
                    ppb=plevs(k)*10.
                    ppb=max(zero,min(ppb,r2000))
                    if(ppb>=etabl_ps(itypex,1,1)) k1_ps=1
                    do kl=1,32
                       if(ppb>=etabl_ps(itypex,kl+1,1).and.ppb<=etabl_ps(itypex,kl,1)) k1_ps=kl
                    end do
                    if(ppb<=etabl_ps(itypex,33,1)) k1_ps=5
                    k2_ps=k1_ps+1
                    ediff_ps = etabl_ps(itypex,k2_ps,1)-etabl_ps(itypex,k1_ps,1)
                    if (abs(ediff_ps) > tiny_r_kind) then
                       del_ps = (ppb-etabl_ps(itypex,k1_ps,1))/ediff_ps
                    else
                      del_ps = huge_r_kind
                    endif
                    del_ps=max(zero,min(del_ps,one))
                    if(oberrflg)then
!                      write(6,*) 'READ_HDRAOB_PS:',itypex,k1_ps,ierr_ps,k2_ps,ierr_ps
                       obserr(1,k)=(one-del_ps)*etabl_ps(itypex,k1_ps,ierr_ps)+del_ps*etabl_ps(itypex,k2_ps,ierr_ps)
!rface pressure error
                       obserr(1,k)=max(obserr(1,k),perrmin)
                    endif
! Surface pressure b
                   var_jb(1,k)=(one-del_ps)*btabl_ps(itypex,k1_ps,ierr_ps)+del_ps*btabl_ps(itypex,k2_ps,ierr_ps)
                    var_jb(1,k)=max(var_jb(1,k),pjbmin)
                    if (var_jb(1,k) >=10.0_r_kind) var_jb(1,k)=zero
                 enddo
              endif
              if (tob) then
                 itypex=itypey
                 ierr_t=0
                 do i =1,maxsub_t
                    if( icsubtype(nc) == isuble_t(itypex,i) ) then
                       ierr_t=i+1
                       exit
                    else if( i == maxsub_t .and. icsubtype(nc) /= isuble_t(itypex,i)) then
                       ncount_t=ncount_t+1
                       do j=1,maxsub_t
                          if(isuble_t(itypex,j) ==0 ) then
                             ierr_t=j+1
                             exit
                          endif
                       enddo
                       if( ncount_t ==1) then
                          write(6,*) 'READ_HDRAOB,WARNING!! tob:cannot find subtyep in the error,& 
                                      table,itype,iosub=',itypex,icsubtype(nc)
                          write(6,*) 'read error table at colomn subtype as 0,error table column=',ierr_t
                       endif
                    endif
                 enddo
                 do k=1,levs
                    ppb=plevs(k)*10.
                    if(ppb <= zero .or. ppb > r2000)cycle loop_readsb
                    if(ppb>=etabl_t(itypex,1,1)) k1_t=1
                       do kl=1,32
                       if(ppb>=etabl_t(itypex,kl+1,1).and.ppb<=etabl_t(itypex,kl,1)) k1_t=kl
                    end do
                    if(ppb<=etabl_t(itypex,33,1)) k1_t=5
                    k2_t=k1_t+1
                    ediff_t = etabl_t(itypex,k2_t,1)-etabl_t(itypex,k1_t,1)
                    if (abs(ediff_t) > tiny_r_kind) then
                       del_t = (ppb-etabl_t(itypex,k1_t,1))/ediff_t
                    else
                      del_t = huge_r_kind
                    endif
                    del_t=max(zero,min(del_t,one))
! Temperature error
                    if(oberrflg)then
!                      write(6,*) 'READ_HDRAOB_T:',itypex,k1_t,itypey,k2_t,ierr_t,nc,kx,ppb
                       obserr(3,k)=(one-del_t)*etabl_t(itypex,k1_t,ierr_t)+del_t*etabl_t(itypex,k2_t,ierr_t)
                       obserr(3,k)=max(obserr(3,k),terrmin)
                    endif
!Temperature b
                    var_jb(3,k)=(one-del_t)*btabl_t(itypex,k1_t,ierr_t)+del_t*btabl_t(itypex,k2_t,ierr_t)
                    var_jb(3,k)=max(var_jb(3,k),tjbmin)
                    if (var_jb(3,k) >=10.0_r_kind) var_jb(3,k)=zero
                 enddo
              endif
              if (qob) then
                 itypex=itypey
                 ierr_q=0
                 do i =1,maxsub_q
                    if( icsubtype(nc) == isuble_q(itypex,i) ) then
                       ierr_q=i+1
                       exit
                    else if( i == maxsub_q .and. icsubtype(nc) /= isuble_q(itypex,i)) then
                       ncount_q=ncount_q+1
                       do j=1,maxsub_q
                          if(isuble_q(itypex,j) ==0 ) then
                             ierr_q=j+1
                             exit
                          endif
                       enddo
                       if(ncount_q ==1 ) then
                          write(6,*) 'READ_HDRAOB,WARNING!! qob:cannot find subtyep in the & 
                                     error table,itype,iosub=',itypex,icsubtype(nc)
                          write(6,*) 'read error table at colomn subtype as 0,error table column=',ierr_q
                       endif
                    endif
                 enddo
                 do k=1,levs
                    ppb=plevs(k)
                    ppb=max(zero,min(ppb,r2000))
                    if(ppb>=etabl_q(itypex,1,1)) k1_q=1
                    do kl=1,32
                       if(ppb>=etabl_q(itypex,kl+1,1).and.ppb<=etabl_q(itypex,kl,1)) k1_q=kl
                    end do
                    if(ppb<=etabl_q(itypex,33,1)) k1_q=5
                    k2_q=k1_q+1
                    ediff_q = etabl_q(itypex,k2_q,1)-etabl_q(itypex,k1_q,1)
                    if (abs(ediff_q) > tiny_r_kind) then
                       del_q = (ppb-etabl_q(itypex,k1_q,1))/ediff_q
                    else
                      del_q = huge_r_kind
                    endif
                    del_q=max(zero,min(del_q,one))
! Humidity error
                    if(oberrflg)then
!                       write(6,*) 'READ_HDRAOB_Q:',itypex,k1_q,itypey,k2_q,ierr_q,nc,kx,ppb
                       obserr(2,k)=(one-del_q)*etabl_q(itypex,k1_q,ierr_q)+del_q*etabl_q(itypex,k2_q,ierr_q)
                       obserr(2,k)=max(obserr(2,k),qerrmin)
                    endif
!Humidity b
                    var_jb(2,k)=(one-del_q)*btabl_q(itypex,k1_q,ierr_q)+del_q*btabl_q(itypex,k2_q,ierr_q)
                    var_jb(2,k)=max(var_jb(2,k),qjbmin)
                    if (var_jb(2,k) >=10.0_r_kind) var_jb(2,k)=zero
!                   if(itypey==120  ) then
!                     write(6,*) 'READ_HDRAOB:120_q,obserr,var_jb=',obserr(2,k),var_jb(2,k),ppb
!                   endif
                 enddo
             endif
             if (uvob) then
                itypex=itypey
                ierr_uv=0
                do i =1,maxsub_uv
                   if( icsubtype(nc) == isuble_uv(itypex,i) ) then
                       ierr_uv=i+1
                       exit
                    else if( i == maxsub_uv .and. icsubtype(nc) /= isuble_uv(itypex,i)) then
                       ncount_uv=ncount_uv+1
                       do j=1,maxsub_uv
                          if(isuble_uv(itypex,j) ==0 ) then
                             ierr_uv=j+1
                             exit
                          endif
                       enddo
                       if( ncount_uv == 1) then
                          write(6,*) 'READ_HDRAOB,WARNING!! uvob:cannot find subtyep in the error,&
                                      table,itype,iosub=',itypex,icsubtype(nc)
                          write(6,*) 'read error table at colomn subtype as 0,error table column=',ierr_uv
                       endif
                    endif
                 enddo
!  This has to be redone since no pressure for the ob. available.
                 do k=1,levs
                   ppb=plevs(k)*10.
                   ppb=max(zero,min(ppb,r2000))
                   if(ppb>=etabl_uv(itypex,1,1)) k1_uv=1
                   do kl=1,32
                      if(ppb>=etabl_uv(itypex,kl+1,1).and.ppb<=etabl_uv(itypex,kl,1)) k1_uv=kl
                   end do
                   if(ppb<=etabl_uv(itypex,33,1)) k1_uv=5
                   k2_uv=k1_uv+1
                   ediff_uv = etabl_uv(itypex,k2_uv,1)-etabl_uv(itypex,k1_uv,1)
                   if (abs(ediff_uv) > tiny_r_kind) then
                      del_uv = (ppb-etabl_uv(itypex,k1_uv,1))/ediff_uv
                   else
                     del_uv = huge_r_kind
                   endif
                   del_uv=max(zero,min(del_uv,one))
! Wind error
!                         write(6,*) 'READ_HDRAOB_UV:',itypex,k1_uv,itypey,k2_uv,ierr_uv,nc,kx,ppb
                   obserr(5,k)=(one-del_uv)*etabl_uv(itypex,k1_uv,ierr_uv)+del_uv*etabl_uv(itypex,k2_uv,ierr_uv)
                   obserr(5,k)=max(obserr(5,k),werrmin)
!Wind b
                   var_jb(5,k)=(one-del_uv)*btabl_uv(itypex,k1_uv,ierr_uv)+del_uv*btabl_uv(itypex,k2_uv,ierr_uv)
                   var_jb(5,k)=max(var_jb(5,k),wjbmin)
                   if (var_jb(5,k) >=10.0_r_kind) var_jb(5,k)=zero
!                      if(itypey==220) then
!                         write(6,*) 'READ_HDRAOB:220_uv,obserr,var_jb=',obserr(5,k),var_jb(5,k),ppb,k2_uv,del_uv
!                      endif
                enddo
             endif
           else
             do k=1,levs
                itypex=kx
                ppb=plevs(k)*10.
                if(ppb <= zero .or. ppb > r2000)cycle loop_readsb
                if(ppb>=etabl(itypex,1,1)) k1=1
                do kl=1,32
                   if(ppb>=etabl(itypex,kl+1,1).and.ppb<=etabl(itypex,kl,1)) k1=kl
                end do
                if(ppb<=etabl(itypex,33,1)) k1=5
                k2=k1+1
                if (abs(ediff) > tiny_r_kind) then
                   del = (ppb-etabl(itypex,k1,1))/ediff
                else
                   del = huge_r_kind
                endif
                del=max(zero,min(del,one))
                obserr(3,k)=(one-del)*etabl(itypex,k1,2)+del*etabl(itypex,k2,2)
                obserr(2,k)=(one-del)*etabl(itypex,k1,3)+del*etabl(itypex,k2,3)
                obserr(5,k)=(one-del)*etabl(itypex,k1,4)+del*etabl(itypex,k2,4)
                obserr(1,k)=(one-del)*etabl(itypex,k1,5)+del*etabl(itypex,k2,5)
                obserr(7,k)=(one-del)*etabl(itypex,k1,6)+del*etabl(itypex,k2,6)

                obserr(3,k)=max(obserr(3,k),terrmin)
                obserr(2,k)=max(obserr(2,k),qerrmin)
                obserr(5,k)=max(obserr(5,k),werrmin)
                obserr(1,k)=max(obserr(1,k),perrmin)
             enddo
           endif      ! endif for njqc



!          If temperature ob, extract information regarding virtual
!          versus sensible temperature


!          write(6,*) tob,qob,levs
           stnelev=hdr3(1)
           if(abs(hdr2(7)) > r90 .or. abs(hdr2(8)) > 720._r_kind)then 
              write(6,*) ' invalid lat,lon ',id,obstype,hdr2(7),hdr2(8)
              cycle loop_readsb
           end if
           if(psob)levs=1
           LOOP_K_LEVS: do k=1,levs
              if(tob .or. qob .or. psob)then
                 if(levdat(1,k) < 1. .or. levdat(1,k) > 1500.)then
                    write(6,*) ' invalid pressure ',id,k,levs,levdat(1,k)
                    cycle LOOP_K_LEVS
                 end if
              end if
              if(obsdat(1,k) < 0. .or. obsdat(1,k) > 900000.) then
                  write(6,*) ' invalid change in time ',id,obstype,k,obsdat(1,k) 
                  obsdat(1,k) = 0.
              end if
              if(abs(obsdat(2,k)) < 10. .and. abs(obsdat(3,k)) < 10.) then
                 dlat_earth_deg=hdr2(7)+obsdat(2,k)
                 dlon_earth_deg=hdr2(8)+obsdat(3,k)
              else
                 write(6,*) ' invalid change in lat/lon ',id,k,obstype,obsdat(2,k),obsdat(3,k)
                 dlat_earth_deg=hdr2(7)
                 dlon_earth_deg=hdr2(8)
              end if
              if(dlon_earth_deg >= r360)dlon_earth_deg=dlon_earth_deg-r360
              if(dlon_earth_deg < zero)dlon_earth_deg=dlon_earth_deg+r360
              dlon_earth=dlon_earth_deg*deg2rad
              dlat_earth=dlat_earth_deg*deg2rad

!             write(6,*) k,dlat_earth_deg,dlon_earth_deg,hdr2(8),obsdat(3,k)
              if(regional)then
                 call tll2xy(dlon_earth,dlat_earth,dlon,dlat,outside)    ! convert to rotated coordinate
                 if(outside) cycle loop_readsb   ! check to see if outside regional domain
              else
                 dlat = dlat_earth
                 dlon = dlon_earth
                 call grdcrd1(dlat,rlats,nlat,1)
                 call grdcrd1(dlon,rlons,nlon,1)
              endif
!             Interpolate guess pressure profile to observation location
              klon1= int(dlon);  klat1= int(dlat)
              dx   = dlon-klon1; dy   = dlat-klat1
              dx1  = one-dx;     dy1  = one-dy
              w00=dx1*dy1; w10=dx1*dy; w01=dx*dy1; w11=dx*dy

              klat1=min(max(1,klat1),nlat); klon1=min(max(0,klon1),nlon)
              if (klon1==0) klon1=nlon
              klatp1=min(nlat,klat1+1); klonp1=klon1+1
              if (klonp1==nlon+1) klonp1=1


!             if(qob)write(6,*) tob,qob,dlnpob,dlat,dlon

              icntpnt=icntpnt+1

!             Set usage variable              
              usage = zero
              if(icuse(nc) <= 0)usage=100._r_kind

              if(plevs(k) < 0.0001_r_kind) then
                 write(*,*) 'warning: obs pressure is too small:',kx,k,plevs(k),id
                 cycle
              endif
              dlnpob=log(plevs(k))  ! ln(pressure in cb)

              ndata=ndata+1
              nodata=nodata+1
              if(uvob)nodata=nodata+1
              iout=ndata
              isort(icntpnt)=iout


              if(ndata > maxobs) then
                 write(6,*)'READ_HDRAOB:  ***WARNING*** ndata > maxobs for ',obstype
                 ndata = maxobs
              end if

              call deter_sfc2(dlat_earth,dlon_earth,t4dv,idomsfc,tsavg,ff10,sfcr,zz)

!             Extract pressure level and quality marks

!             Temperature
              if(tob) then
                 if(obsdat(7,k) < 100. .or. obsdat(7,k) > 400.) then
                    write(6,*)id,'invalid temp',k,levs,obsdat(7,k),plevs(k)
                    cycle LOOP_K_LEVS
                 end if
                 ppb=plevs(k)*10.
!  setup later
                 if(levs > 100 .or. plevs(1)-plevs(levs) < .01)then

                   call errormod_hdraob(pqm,tqm,levs,plevs,errout,k,presl,dpres,nsig,lim_qm)
      
                 else

                   call errormod(pqm,tqm,levs,plevs,errout,k,presl,dpres,nsig,lim_qm)
                 end if
                 toe=obserr(3,k)*errout
                 if(ppb < r100)toe=toe*r1_2
                 cdata_all(1,iout)=toe                     ! temperature error
                 cdata_all(2,iout)=dlon                    ! grid relative longitude
                 cdata_all(3,iout)=dlat                    ! grid relative latitude
                 cdata_all(4,iout)=dlnpob                  ! ln(pressure in cb)

                 cdata_all(5,iout)=obsdat(7,k)             ! temperature ob.
                 cdata_all(6,iout)=rstation_id             ! station id
                 cdata_all(7,iout)=t4dv+obsdat(1,k)*r60inv ! time
                 cdata_all(8,iout)=nc                      ! type
                 qtflg=one
                 cdata_all(9,iout)=qtflg                   ! qtflg (virtual temperature flag)
                 cdata_all(10,iout)=0                      ! quality mark
                 cdata_all(11,iout)=obserr(3,k)            ! original obs error            
                 cdata_all(12,iout)=usage                  ! usage parameter
                 cdata_all(13,iout)=idomsfc                ! dominate surface type
                 cdata_all(14,iout)=tsavg                  ! skin temperature
                 cdata_all(15,iout)=ff10                   ! 10 meter wind factor
                 cdata_all(16,iout)=sfcr                   ! surface roughness
                 cdata_all(17,iout)=dlon_earth_deg         ! earth relative longitude (degrees)
                 cdata_all(18,iout)=dlat_earth_deg         ! earth relative latitude (degrees)
                 cdata_all(19,iout)=stnelev                ! station elevation (m)
                 cdata_all(20,iout)=levdat(2,k)            ! observation height (m)
                 cdata_all(21,iout)=zz                     ! terrain height at ob location
                 cdata_all(22,iout)='88888888'             ! provider name
                 cdata_all(23,iout)='88888888'             ! subprovider name
                 cdata_all(24,iout)=2                      ! cat
                 cdata_all(25,iout)=var_jb(3,k)            ! non linear qc for T
                 if(perturb_obs)cdata_all(nreal,iout)=ran01dom()*perturb_fact ! t perturbation
                 if (twodvar_regional) &
                    call adjust_error(cdata_all(17,iout),cdata_all(18,iout),cdata_all(11,iout),cdata_all(1,iout))

!             Winds 
              else if(uvob) then 
                 if(obsdat(6,k) < 0. .or. obsdat(6,k) > 360.) then
                    write(6,*)id,'invalid dir ', k,levs,obsdat(5,k),obsdat(6,k),plevs(k)
                    cycle LOOP_K_LEVS
                 end if
                 if(obsdat(5,k) < 0. .or. obsdat(5,k) > 200.) then
                    write(6,*)id,'invalid spd ', k,levs,obsdat(5,k),obsdat(6,k),plevs(k)
                    cycle LOOP_K_LEVS
                 end if
                 if(levs > 100 .or. plevs(1)-plevs(levs) < .01)then
                    call errormod_hdraob(pqm,qqm,levs,plevs,errout,k,presl,dpres,nsig,lim_qm)
                 else

                    call errormod(pqm,wqm,levs,plevs,errout,k,presl,dpres,nsig,lim_qm)
                 end if
                 woe=obserr(5,k)*errout
                 if(obsdat(1,k) < r50)woe=woe*r1_2
                 oelev=levdat(2,k)

!                Rotate winds to rotated coordinate
                 uob  = -obsdat(5,k)*sin(obsdat(6,k)*deg2rad) ! u-wind component
                 vob  = -obsdat(5,k)*cos(obsdat(6,k)*deg2rad) ! v-wind component


                 if(regional .and. .not. fv3_regional)then
                    u0=uob
                    v0=vob
                    call rotate_wind_ll2xy(u0,v0,uob,vob,dlon_earth,dlon,dlat)
                    if(diagnostic_reg) then
                       call rotate_wind_xy2ll(uob,vob,u00,v00,dlon_earth,dlon,dlat)
                       nvtest=nvtest+1
                       disterr=sqrt((u0-u00)**2+(v0-v00)**2)
                       vdisterrmax=max(vdisterrmax,disterr)
                    end if
                 endif

                 cdata_all(1,iout)=woe                     ! wind error
                 cdata_all(2,iout)=dlon                    ! grid relative longitude
                 cdata_all(3,iout)=dlat                    ! grid relative latitude
                 cdata_all(4,iout)=dlnpob                  ! ln(pressure in cb)
                 cdata_all(5,iout)=oelev                   ! height of observation
                 cdata_all(6,iout)=uob                     ! u obs
                 cdata_all(7,iout)=vob                     ! v obs
                 cdata_all(8,iout)=rstation_id             ! station id
                 cdata_all(9,iout)=t4dv+obsdat(1,k)*r60inv ! time
                 cdata_all(10,iout)=nc                     ! type
                 cdata_all(11,iout)=stnelev                ! station elevation
                 cdata_all(12,iout)=0                      ! quality mark
                 cdata_all(13,iout)=obserr(5,k)            ! original obs error
                 cdata_all(14,iout)=usage                  ! usage parameter
                 cdata_all(15,iout)=idomsfc                ! dominate surface type
                 cdata_all(16,iout)=tsavg                  ! skin temperature
                 cdata_all(17,iout)=ff10                   ! 10 meter wind factor
                 cdata_all(18,iout)=sfcr                   ! surface roughness
                 cdata_all(19,iout)=dlon_earth_deg         ! earth relative longitude (degrees)
                 cdata_all(20,iout)=dlat_earth_deg         ! earth relative latitude (degrees)
                 cdata_all(21,iout)=zz                     ! terrain height at ob location
                 cdata_all(22,iout)='88888888'             ! provider name
                 cdata_all(23,iout)='HDRAOB'               ! subprovider name
                 cdata_all(24,iout)=2                      ! cat
                 cdata_all(25,iout)=var_jb(5,k)            ! non linear qc parameter
                 cdata_all(26,iout)=one                    ! hilbert curve weight, modified later 
                 if(perturb_obs)then
                    cdata_all(27,iout)=ran01dom()*perturb_fact ! u perturbation
                    cdata_all(28,iout)=ran01dom()*perturb_fact ! v perturbation
                 endif
 
!             Specific humidity 
              else if(qob) then
!  setup later
                 if(levs > 100 .or. plevs(1)-plevs(levs) < .01)then
                    call errormod_hdraob(pqm,qqm,levs,plevs,errout,k,presl,dpres,nsig,lim_qm)
                 else

                    call errormod(pqm,qqm,levs,plevs,errout,k,presl,dpres,nsig,lim_qm)
                 end if
                 qoe=obserr(2,k)*errout*one_tenth
                 if(obsdat(7,k) < 215.) then
                    exit LOOP_K_LEVS
                 end if
                 if(obsdat(8,k) < 100. .or. obsdat(8,k) > 300.) then
                    write(6,*)id,'invalid td ', k,levs,obsdat(7,k),obsdat(8,k),plevs(k)
                    cycle LOOP_K_LEVS
                 end if
!   Need to convert from td to q
                 call fpvsx_ad(obsdat(8,k),es,dummy,dummy,.false.)
                 qobcon = eps*es/(plevs(k)-omeps*es)
                 cdata_all(1,iout)=qoe                     ! q error   
                 cdata_all(2,iout)=dlon                    ! grid relative longitude
                 cdata_all(3,iout)=dlat                    ! grid relative latitude
                 cdata_all(4,iout)=dlnpob                  ! ln(pressure in cb)
                 cdata_all(5,iout)=qobcon                  ! q ob
                 cdata_all(6,iout)=rstation_id             ! station id
                 cdata_all(7,iout)=t4dv+obsdat(1,k)*r60inv ! time
                 cdata_all(8,iout)=nc                      ! type
                 cdata_all(9,iout)=0.2_r_kind              ! q max error
                 cdata_all(10,iout)=obsdat(7,k)            ! dry temperature (obs is tv)
                 cdata_all(11,iout)=0                      ! quality mark
                 cdata_all(12,iout)=obserr(2,k)*one_tenth  ! original obs error
                 cdata_all(13,iout)=usage                  ! usage parameter
                 cdata_all(14,iout)=idomsfc                ! dominate surface type
                 cdata_all(15,iout)=dlon_earth_deg         ! earth relative longitude (degrees)
                 cdata_all(16,iout)=dlat_earth_deg         ! earth relative latitude (degrees)
                 cdata_all(17,iout)=stnelev                ! station elevation (m)
                 cdata_all(18,iout)=levdat(2,k)            ! observation height (m)
                 cdata_all(19,iout)=zz                     ! terrain height at ob location
                 cdata_all(20,iout)='88888888'             ! provider name
                 cdata_all(21,iout)='HDRAOB'               ! subprovider name
                 cdata_all(22,iout)=2                      ! cat
                 cdata_all(23,iout)=var_jb(2,k)            ! non linear qc b parameter
                 if(perturb_obs)cdata_all(24,iout)=ran01dom()*perturb_fact ! q perturbation
                 if (twodvar_regional) &
                    call adjust_error(cdata_all(15,iout),cdata_all(16,iout),cdata_all(12,iout),cdata_all(1,iout))
 
              else if(psob) then
                 if(obsdat(7,k) > 400. .or. stnelev > 7000.) exit LOOP_K_LEVS
                 poe=obserr(1,k)*one_tenth                !  convert from mb to cb
                 cdata_all(1,iout)=poe                     ! surface pressure error (cb)
                 cdata_all(2,iout)=dlon                    ! grid relative longitude
                 cdata_all(3,iout)=dlat                    ! grid relative latitude

                 cdata_all(4,iout)=plevs(k)                ! pressure (in cb)

                 cdata_all(5,iout)=obsdat(4,k)/grav        ! surface height
                 cdata_all(6,iout)=obsdat(7,k)             ! surface temperature
                 cdata_all(7,iout)=rstation_id             ! station id
                 cdata_all(8,iout)=t4dv                    ! time
                 cdata_all(9,iout)=nc                      ! type
                 cdata_all(10,iout)=0                      ! quality mark
                 cdata_all(11,iout)=obserr(1,k)*one_tenth  ! original obs error (cb)
                 cdata_all(12,iout)=usage                  ! usage parameter
                 cdata_all(13,iout)=idomsfc                ! dominate surface type
                 cdata_all(14,iout)=dlon_earth_deg         ! earth relative longitude (degrees)
                 cdata_all(15,iout)=dlat_earth_deg         ! earth relative latitude (degrees)
                 cdata_all(16,iout)=stnelev                ! station elevation (m)
                 cdata_all(17,iout)=zz                     ! terrain height at ob location
                 cdata_all(18,iout)='88888888'             ! provider name
                 cdata_all(19,iout)='HDRAOB'               ! subprovider name
                 cdata_all(20,iout)=var_jb(1,k)            ! non linear qc b parameter
                 if(perturb_obs)cdata_all(21,iout)=ran01dom()*perturb_fact ! ps perturbation
                 if (twodvar_regional) &
                    call adjust_error(cdata_all(14,iout),cdata_all(15,iout),cdata_all(11,iout),cdata_all(1,iout)) 

                 write(6,*) k,kx,c_station_id,(cdata_all(i,iout),i=1,20)
              end if

!
!    End k loop over levs
           end do  LOOP_K_LEVS
        end do loop_readsb

!
!   End of bufr read loop
     enddo loop_msg
!    Close unit to bufr file

! Normal exit

  enddo loop_convinfo! loops over convinfo entry matches
  deallocate(lmsg,tab,nrep)
  call closbf(lunin)

  if(print_verbose)write(6,*)'READ_HDRAOB:  closbf(',lunin,')'

  close(lunin)

! Write header record and data to output file for further processing
  allocate(iloc(ndata))
  icount=0
  do i=1,maxobs
     if(isort(i) > 0)then
       icount=icount+1
       iloc(icount)=isort(i)
     end if
  end do
  write(6,*) 'ndata,icount',ndata,icount
  if(ndata /= icount)then
     write(6,*) ' HDRAOB: mix up in read_prepbufr ,ndata,icount ',ndata,icount
     call stop2(50)
  end if
  allocate(cdata_out(nreal,ndata))
  do i=1,ndata
     itx=iloc(i)
     do k=1,nreal
        cdata_out(k,i)=cdata_all(k,itx)
     end do
  end do
  deallocate(iloc,isort,cdata_all)

  call count_obs(ndata,nreal,ilat,ilon,cdata_out,nobs)
  write(lunout) obstype,sis,nreal,nchanl,ilat,ilon,ndata
  write(lunout) cdata_out
 

  deallocate(cdata_out)

  if(diagnostic_reg .and. ntest>0) write(6,*)'READ_HDRAOB:  ',&
     'ntest,disterrmax=',ntest,disterrmax
  if(diagnostic_reg .and. nvtest>0) write(6,*)'READ_HDRAOB:  ',&
     'nvtest,vdisterrmax=',ntest,vdisterrmax

! End of routine
  return

end subroutine read_hdraob

