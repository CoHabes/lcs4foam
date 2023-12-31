!
!Copyright (C) 2015-2016, Justin R. Finn.  All rights reserved.
!libcfd2lcs is distributed is under the terms of the GNU General Public License
!
module io_m
      use data_m
      use structured_m
      implicit none
      !----
      !Read and write routine for r0,r1,r2 thedatatypes
      !structured_io is a many to many write
      !write_ftle is a write routine for ftle and flow map resutls to a single file
      !----

      integer,parameter:: &
            IO_READ = 1, &
            IO_WRITE = 2, &
            IO_APPEND = 3

      integer,parameter:: &
            R0_DATA = 0, &
            R1_DATA = 1, &
            R2_DATA = 2

      character(len=25),parameter:: ACTION_STRING(3)= (/&
            ' READING         : ', &
            ' WRITING         : ', &
            ' WRITING (APPEND): ' /)

      character(len=4),parameter:: FILE_EXT = '.dat'  !DONT TOUCH THIS!!!

      !Set the format of the tmp files:
      !Choose ASCII for human readable with a tecplot header
      !Choose BINARY for speed
      integer,parameter::&
            BINARY = 0,&
            ASCII = 1
      integer,parameter:: TMP_FILE_FORMAT = BINARY

      !Parameters for writing a single ascii file:
      integer,parameter:: MAX_CHUNK = 100000
      integer,parameter:: NDATA_PER_LINE = 10  !MAKE SURE MAX_CHUNK is divisible by this
      character(len=6),parameter:: ASCII_FMT = 'ES14.6'

      contains

      subroutine write_lcs(lcs,time)
            implicit none
            !-----
            type(lcs_t):: lcs
            real(LCSRP):: time
            !-----
            integer:: gn(3),offset(3)
            character(len=128):: FMT1,fname
            integer(8):: findex
            real:: t0,t1
            !-----
            !Output a datafile containing the LCS diagnostic results
            !-----

            t0 = cputimer(lcscomm,SYNC_TIMER)

            !-----
            !Generate the filename.
            !Note, the index corresponds to time T=T0
            !-----
            findex = nint(time/lcs%h,8) !file index at current time
            !if(lcs%diagnostic == FTLE_FWD .or. lcs%diagnostic == LP_TRACER) then
            if(lcs%diagnostic == FTLE_FWD .or. lcs%diagnostic == FTLE_BKWD) then
                  findex = findex - nint(lcs%T/lcs%h,8) !File index at t=t0 for fwd time diagnostic
            endif

            select case(findex)
                  case(-999999999:-100000000)
                        FMT1 = "(a,a,a,i10.9,a)"
                  case(-99999999:-10000000)
                        FMT1 = "(a,a,a,i9.8,a)"
                  case(-9999999:-1000000)
                        FMT1 = "(a,a,a,i8.7,a)"
                  case(-999999:-100000)
                        FMT1 = "(a,a,a,i7.6,a)"
                  case(-99999:-10000)
                        FMT1 = "(a,a,a,i6.5,a)"
                  case(-9999:-1000)
                        FMT1 = "(a,a,a,i5.4,a)"
                  case(-999:-100)
                        FMT1 = "(a,a,a,i4.3,a)"
                  case(-99:-10)
                        FMT1 = "(a,a,a,i3.2,a)"
                  case(-9:-1)
                        FMT1 = "(a,a,a,i2.1,a)"
                  case(0:9)
                        FMT1 = "(a,a,a,i1.1,a)"
                  case(10:99)
                        FMT1 = "(a,a,a,i2.2,a)"
                  case(100:999)
                        FMT1 = "(a,a,a,i3.3,a)"
                  case(1000:9999)
                        FMT1 = "(a,a,a,i10.10,a)"
                  case(10000:99999)
                        FMT1 = "(a,a,a,i5.5,a)"
                  case(100000:999999)
                        FMT1 = "(a,a,a,i6.6,a)"
                  case(1000000:9999999)
                        FMT1 = "(a,a,a,i7.7,a)"
                  case(10000000:99999999)
                        FMT1 = "(a,a,a,i8.8,a)"
                  case(100000000:999999999)
                        FMT1 = "(a,a,a,i9.9,a)"
                  case default
                        if(lcsrank==0)&
                              write(*,*) 'ERROR: unsupported range for file index,',findex
                        CFD2LCS_ERROR = 1
                        return
            end select
            write(fname,trim(FMT1))'./cfd2lcs_output/',trim(lcs%label),'_',findex,FILE_EXT
            if(lcsrank==0)&
                  write(*,*) 'In write_lcs...',trim(fname)

            !-----
            !Ouptut the LCS.  Data depends on diagnostic type.
            !-----
            select case (lcs%diagnostic)
                  case(FTLE_FWD,FTLE_BKWD)
                        gn = (/lcs%sgrid%gni,lcs%sgrid%gnj,lcs%sgrid%gnk/)
                        offset = (/lcs%sgrid%offset_i,lcs%sgrid%offset_j,lcs%sgrid%offset_k/)
                        !Special call to write grid,fm,ftle to ascii file
                        call write_ftle(trim(fname),gn,offset,lcs)
                  case(LP_TRACER)
                        call write_tracers(trim(fname),lcs%lp)
                  case default
            end select

            t1 = cputimer(lcscomm,SYNC_TIMER)
            cpu_io = cpu_io + max(t1-t0,0.0)

      end subroutine write_lcs
      
      subroutine write_tracers(fname,lp)
            implicit none
            !-----
            character(len=*):: fname
            type(lp_t):: lp
            !-----
            integer:: funit = 22
            integer:: ierr,proc,ip
            integer,parameter:: NVAR = 3
            character(len = 12):: myfmt
            !-----
            
            !rank 0 writes the header:
            if (lcsrank==0) then
                  open(funit,file=fname,status='replace',form='formatted')
                  write(funit,'(a,a,a)') 'TITLE = "',trim(lp%label),'"'
                  write(funit,'(a)') 'VARIABLES ='
                  write(funit,'(a)') '"X"'
                  write(funit,'(a)') '"Y"'
                  write(funit,'(a)') '"Z"'
                  write(funit,'(a)') 'ZONE'
                  write(funit,'(a,a,a)') 'T = "',trim(lp%label),'"'
                  write(funit,'(a,i10.10,a,ES18.4)') 'STRANDID = ',lp%id+1,', SOLUTIONTIME = ',scfd%t_np1
                  write(funit,'(a,i10.10,a,i10.10,a,i10.10)') 'I = ',lp%npall,' J = ',1, ' K = ',1
                  write(funit,'(a)') 'DATAPACKING = POINT'
                  close(funit)
            endif

            !Write the data, beware this will get sloooow for lots of procs
            write(myfmt,'(a,i1,a,a)') '(',NVAR,trim(ASCII_FMT),')'
            do proc = 0,nprocs-1
                  if(lcsrank==proc) then
                        if(lp%np>0) then
                              open(funit,file=fname,status='unknown',form='formatted',position='append')
                              do ip = 1,lp%np
                                    write(funit,trim(myfmt)) lp%xp%x(ip),lp%xp%y(ip),lp%xp%z(ip)
                              enddo
                              close(funit)
                        endif
                  endif
                  call mpi_barrier(lcscomm,ierr)
            enddo
      end subroutine write_tracers

      subroutine write_ftle(fname,gn,offset,lcs)
            implicit none
            !-----
            character(len=*):: fname
            integer,dimension(3):: gn
            integer,dimension(3):: offset
            type(lcs_t):: lcs
            !-----
            integer:: ni,nj,nk,ng
            integer:: funit = 22
            integer,parameter:: FTLE_STRAND = 1 
            !-----

            !For brevity:
            ni = lcs%sgrid%grid%ni
            nj = lcs%sgrid%grid%nj
            nk = lcs%sgrid%grid%nk
            ng = lcs%sgrid%grid%ng

            !rank 0 writes the header:
            if (lcsrank==0) then
                  open(funit,file=fname,status='replace',form='formatted')
                  write(funit,'(a,a,a)') 'TITLE = "',trim(lcs%ftle%label),'"'
                  write(funit,'(a)') 'VARIABLES ='
                  write(funit,'(a)') '"X"'
                  write(funit,'(a)') '"Y"'
                  write(funit,'(a)') '"Z"'
                  if(FLOWMAP_IO) then
                        write(funit,'(a)') '"FM-X"'
                        write(funit,'(a)') '"FM-Y"'
                        write(funit,'(a)') '"FM-Z"'
                  endif
                  if(VELOCITY_IO) then
                        write(funit,'(a)') '"U-X"'
                        write(funit,'(a)') '"U-Y"'
                        write(funit,'(a)') '"U-Z"'
                  endif
                  if(BCFLAG_IO) then
                        write(funit,'(a)') '"FLAG"'
                  endif
                  if(VELOCITY_INVARIANTS) then
                        write(funit,'(a)') '"Q"'
                        write(funit,'(a)') '"C"'
                        write(funit,'(a)') '"D"'
                        write(funit,'(a)') '"H"'
                        write(funit,'(a)') '"L2"'
                  endif
                  write(funit,'(a)') '"FTLE"'
                  write(funit,'(a)') 'ZONE'
                  write(funit,'(a,a,a)') 'T = "',trim(lcs%ftle%label),'"'
                  write(funit,'(a,i10.10,a,ES18.4)') 'STRANDID = ',FTLE_STRAND,', SOLUTIONTIME = ',scfd%t_np1
                  write(funit,'(a,i10.10,a,i10.10,a,i10.10)') 'I = ',gn(1),' J = ',gn(2), ' K = ',gn(3)
                  write(funit,'(a)') 'DATAPACKING = BLOCK'
                  close(funit)
            endif

            !Write the data
            call write_data_from_master(fname,0,ni,nj,nk,ng,gn,offset,lcs%sgrid%grid%x)
            call write_data_from_master(fname,0,ni,nj,nk,ng,gn,offset,lcs%sgrid%grid%y)
            call write_data_from_master(fname,0,ni,nj,nk,ng,gn,offset,lcs%sgrid%grid%z)
            if(FLOWMAP_IO) then
                  call write_data_from_master(fname,0,ni,nj,nk,ng,gn,offset,lcs%lp%fm%x)
                  call write_data_from_master(fname,0,ni,nj,nk,ng,gn,offset,lcs%lp%fm%y)
                  call write_data_from_master(fname,0,ni,nj,nk,ng,gn,offset,lcs%lp%fm%z)
            endif
            if(VELOCITY_IO) then
                  call write_data_from_master(fname,0,ni,nj,nk,ng,gn,offset,lcs%lp%ugrid%x)
                  call write_data_from_master(fname,0,ni,nj,nk,ng,gn,offset,lcs%lp%ugrid%y)
                  call write_data_from_master(fname,0,ni,nj,nk,ng,gn,offset,lcs%lp%ugrid%z)
            endif
            if(BCFLAG_IO) then
                  call write_data_from_master(fname,0,ni,nj,nk,ng,gn,offset,real(lcs%sgrid%bcflag%i,LCSRP))
            endif
            if(VELOCITY_INVARIANTS) then
               call write_data_from_master(fname,0,ni,nj,nk,ng,gn,offset,lcs%inv%Q%r)
               call write_data_from_master(fname,0,ni,nj,nk,ng,gn,offset,lcs%inv%C%r)
               call write_data_from_master(fname,0,ni,nj,nk,ng,gn,offset,lcs%inv%D%r)
               call write_data_from_master(fname,0,ni,nj,nk,ng,gn,offset,lcs%inv%H%r)
               call write_data_from_master(fname,0,ni,nj,nk,ng,gn,offset,lcs%inv%L2%r)
            endif
            call write_data_from_master(fname,0,ni,nj,nk,ng,gn,offset,lcs%ftle%r)

      end subroutine write_ftle

      subroutine write_data_from_master(fname,master,ni,nj,nk,ng,gn,offset,phi)
            implicit none
            !-----
            character(len=*):: fname
            integer:: ni,nj,nk,ng
            integer:: master,gn(3),offset(3)
            real(LCSRP):: phi(1-ng:ni+ng,1-ng:nj+ng,1-ng:nk+ng)
            !-----
            real(LCSRP)::mybuf(1:MAX_CHUNK),buf(1:MAX_CHUNK)
            integer::i,j,k,ii,jj,kk,ind
            integer:: ierr,write_start, write_end,nwrite,nline,line
            character(len=12):: myfmt
            logical:: last
            integer:: funit=22
            !-----
            !Funnel all the data through a master rank:
            !Work in chunks so that we can limit the memory overhead.
            !-----

            if(lcsrank==master .and. LCS_VERBOSE) &
                  write(*,*) 'in write_data_from_master...'

            ind = 0
            mybuf = 0.0_LCSRP
            last = .false.
            do k = 1, gn(3)
            do j = 1, gn(2)
            do i = 1, gn(1)
                  ind = ind+1
                  if(i==gn(1) .and. j==gn(2) .and. k==gn(3)) last = .true.
                  if(   i > offset(1) .and. i <= offset(1) + ni &
                  .and. j > offset(2) .and. j <= offset(2) + nj &
                  .and. k > offset(3) .and. k <= offset(3) + nk ) then
                        ii = i - offset(1)
                        jj = j - offset(2)
                        kk = k - offset(3)
                        mybuf(ind) = phi(ii,jj,kk)
                  else
                        !do nothing...
                  endif
                  if(ind == MAX_CHUNK .OR. last) then
                        !Reduce to master:
                        call MPI_REDUCE(mybuf(1),buf(1),ind,MPI_LCSRP,MPI_SUM,master,lcscomm,ierr)
                        !Master writes:
                        if(lcsrank==master) then
                              open(funit,file=fname,status='unknown',form='formatted',position='append')
                              nline = ceiling(real(ind,LCSRP)/real(NDATA_PER_LINE))
                              do line = 1,nline
                                    write_start = (line-1)*NDATA_PER_LINE+1
                                    write_end = min(write_start+NDATA_PER_LINE-1,ind)
                                    nwrite = write_end - write_start + 1
                                    if(nwrite<10) write(myfmt,'(a,i1,a,a)') '(',nwrite,trim(ASCII_FMT),')'
                                    if(nwrite>=10) write(myfmt,'(a,i2,a,a)') '(',nwrite,trim(ASCII_FMT),')'
                                    write(funit,trim(myfmt)) buf(write_start:write_end)
                              enddo
                              close(funit)
                        endif
                        !reset the buf:
                        ind = 0
                        mybuf = 0.0_LCSRP
                  endif
            enddo
            enddo
            enddo

      end subroutine write_data_from_master

      subroutine structured_io(fname,IO_ACTION,global_size,offset,r0,r1,r2)
      IMPLICIT NONE
            !-----
            character(len=*):: fname
            integer:: IO_ACTION
            integer,dimension(3):: global_size
            integer,dimension(3):: offset
            type(sr0_t),optional:: r0
            type(sr1_t),optional:: r1
            type(sr2_t),optional:: r2
            !-----
            integer:: NVAR, WORK_DATA
            character(len=LCS_NAMELEN),allocatable:: dataname(:)
            character(len=LCS_NAMELEN):: groupname
            real(LCSRP),allocatable :: thedata(:,:,:)  ! Write buffer
            integer,parameter:: NDIM = 3  !all data considered 3 dimensional
            integer:: ivar
            integer:: ni,nj,nk
            character(len=128):: junk
            integer:: funit
            integer,parameter:: funit_start = 1234
            character(len=20):: myfmt
            real:: t0,t1
            !-----

            t0 = cputimer(lcscomm,SYNC_TIMER)

            !
            !Figure out what type of data we will dump, r0,r1,r2 or grid.
            !Only allow one type to be dumped in a single call to this routine
            !to dump multiple datasets to the same file, call with IO_ACTION = APPEND
            !
            if (present(r0)) then
                  NVAR = 1
                  WORK_DATA = R0_DATA
                  ni = r0%ni; nj = r0%nj; nk = r0%nk

                  write(groupname,'(a)') '/'
                  allocate(dataname(1))
                  write(dataname(1),'(a,a)')  trim(groupname),trim(r0%label)
                  if(lcsrank==0) &
                        write(*,'(a,a,a)') 'In structured_io... ',ACTION_STRING(IO_ACTION),trim(r0%label)
            elseif(present(r1)) then
                  NVAR = 3
                  WORK_DATA = R1_DATA
                  ni = r1%ni; nj = r1%nj; nk = r1%nk
                  write(groupname,'(a)')  trim(r1%label)
                  allocate(dataname(3))
                  write(dataname(1),'(a,a,a)')   trim(groupname),'/','-X'
                  write(dataname(2),'(a,a,a)')   trim(groupname),'/','-Y'
                  write(dataname(3),'(a,a,a)')   trim(groupname),'/','-Z'
                  if(lcsrank==0) &
                        write(*,'(a,a,a)') 'In structured_io... ',ACTION_STRING(IO_ACTION),trim(r1%label)
            elseif(present(r2)) then
                  NVAR = 9
                  WORK_DATA = R2_DATA
                  ni = r2%ni; nj = r2%nj; nk = r2%nk
                  write(groupname,'(a)')  trim(r2%label)
                  allocate(dataname(9))
                  write(dataname(1),'(a,a,a)')   trim(groupname),'/','-XX'
                  write(dataname(2),'(a,a,a)')   trim(groupname),'/','-XY'
                  write(dataname(3),'(a,a,a)')   trim(groupname),'/','-XZ'
                  write(dataname(4),'(a,a,a)')   trim(groupname),'/','-YX'
                  write(dataname(5),'(a,a,a)')   trim(groupname),'/','-YY'
                  write(dataname(6),'(a,a,a)')   trim(groupname),'/','-YZ'
                  write(dataname(7),'(a,a,a)')   trim(groupname),'/','-ZX'
                  write(dataname(8),'(a,a,a)')   trim(groupname),'/','-ZY'
                  write(dataname(9),'(a,a,a)')   trim(groupname),'/','-ZZ'
                  if(lcsrank==0) &
                        write(*,'(a,a,a)') 'In structured_io... ',ACTION_STRING(IO_ACTION),trim(r2%label)
            else
                  if(lcsrank==0) &
                        write(*,'(a,a)') 'In structured_io...  No data present'
                  return
            endif


            !Open for read/write
            !No need to write a header here anymore:
            funit = funit_start + lcsrank
            select case(IO_ACTION)
                  case(IO_WRITE)
                        if(TMP_FILE_FORMAT==ASCII) then
                              open(funit,file=fname,status='replace',form='formatted')
                              write(funit,'(a,a,a)') 'TITLE = "',trim(groupname),'"'
                              write(funit,'(a)') 'VARIABLES ='
                              do ivar = 1,NVAR
                                    write(funit,'(a,a,a)') '"',trim(dataname(ivar)) ,'"'
                              enddo
                              write(funit,'(a)') 'ZONE'
                              write(funit,'(a,a,a,i10.10,a)') 'T = "',trim(groupname),'_',lcsrank,'"'
                              write(funit,'(a,i10.10,a,i10.10,a,i10.10)') 'I = ',ni,' J = ',nj, ' K = ',nk
                              write(funit,'(a)') 'DATAPACKING = BLOCK'
                        elseif (TMP_FILE_FORMAT == BINARY) then
                              open(funit,file=fname,status='replace',form='unformatted')
                        endif
                  case(IO_READ)
                        if(TMP_FILE_FORMAT==ASCII) then
                              open(funit,file=fname,action='read',form='formatted')
                              read(funit,*) junk
                              read(funit,*) junk
                              do ivar = 1,NVAR
                                    read(funit,*) junk
                              enddo
                              read(funit,*) junk
                              read(funit,*) junk
                              read(funit,*) junk
                              read(funit,*) junk
                        elseif (TMP_FILE_FORMAT == BINARY) then
                              open(funit,file=fname,action='read',form='unformatted')
                        endif
                  case default
                        if(lcsrank==0) &
                        write(*,*) 'bad IO_ACTION.  Must be IO_WRITE or IO_READ'
                        CFD2LCS_ERROR = 1
            end select


            !
            ! Allocate read/write buffer
            !
            allocate(thedata(1:ni,1:nj,1:nk))

            !Loop through each variable and read/write the data
            do ivar = 1,NVAR

                  !
                  ! Set the data:
                  !
                  !
                  if(IO_ACTION == IO_WRITE) then
                        select case(ivar)
                        case(1)
                              if(WORK_DATA == R0_DATA) then
                                    thedata(1:ni,1:nj,1:nk) = r0%r(1:ni,1:nj,1:nk)
                              elseif(WORK_DATA == R1_DATA) then
                                    thedata(1:ni,1:nj,1:nk) = r1%x(1:ni,1:nj,1:nk)
                              elseif(WORK_DATA == R2_DATA) then
                                    thedata(1:ni,1:nj,1:nk) = r2%xx(1:ni,1:nj,1:nk)
                              endif
                        case(2)
                              if(WORK_DATA == R0_DATA) then
                                    cycle
                              elseif(WORK_DATA == R1_DATA) then
                                    thedata(1:ni,1:nj,1:nk) = r1%y(1:ni,1:nj,1:nk)
                              elseif(WORK_DATA == R2_DATA) then
                                    thedata(1:ni,1:nj,1:nk) = r2%xy(1:ni,1:nj,1:nk)
                              endif
                        case(3)
                              if(WORK_DATA == R0_DATA) then
                                    cycle
                              elseif(WORK_DATA == R1_DATA) then
                                    thedata(1:ni,1:nj,1:nk) = r1%z(1:ni,1:nj,1:nk)
                              elseif(WORK_DATA == R2_DATA) then
                                    thedata(1:ni,1:nj,1:nk) = r2%xz(1:ni,1:nj,1:nk)
                              endif
                        case(4)
                              thedata(1:ni,1:nj,1:nk) = r2%yx(1:ni,1:nj,1:nk)
                        case(5)
                              thedata(1:ni,1:nj,1:nk) = r2%yy(1:ni,1:nj,1:nk)
                        case(6)
                              thedata(1:ni,1:nj,1:nk) = r2%yz(1:ni,1:nj,1:nk)
                        case(7)
                              thedata(1:ni,1:nj,1:nk) = r2%zx(1:ni,1:nj,1:nk)
                        case(8)
                              thedata(1:ni,1:nj,1:nk) = r2%zy(1:ni,1:nj,1:nk)
                        case(9)
                              thedata(1:ni,1:nj,1:nk) = r2%zz(1:ni,1:nj,1:nk)
                        case default
                              cycle
                        end select

                        !write thedata
                        if(TMP_FILE_FORMAT==ASCII) then
                              write(myfmt,'(a,i8,a,a)') '(',ni*nj*nk,trim(ASCII_FMT),')'
                              write(funit,myfmt) thedata
                        elseif(TMP_FILE_FORMAT==BINARY) then
                              write(funit) thedata
                        endif
                  endif

                  !
                  ! Copy from read buffer
                  !
                  if(IO_ACTION == IO_READ) then
                        !read thedata
                        if(TMP_FILE_FORMAT==ASCII) then
                              write(myfmt,'(a,i8,a,a)') '(',ni*nj*nk,trim(ASCII_FMT),')'
                              read(funit,myfmt) thedata
                        elseif(TMP_FILE_FORMAT==BINARY) then
                              read(funit) thedata
                        endif

                        select case(ivar)
                        case(1)
                              if(WORK_DATA == R0_DATA) then
                                    r0%r(1:ni,1:nj,1:nk) = thedata(1:ni,1:nj,1:nk)
                              elseif(WORK_DATA == R1_DATA) then
                                    r1%x(1:ni,1:nj,1:nk) = thedata(1:ni,1:nj,1:nk)
                              elseif(WORK_DATA == R2_DATA) then
                                    r2%xx(1:ni,1:nj,1:nk) = thedata(1:ni,1:nj,1:nk)
                              endif
                        case(2)
                              if(WORK_DATA == R0_DATA) then
                                    cycle
                              elseif(WORK_DATA == R1_DATA) then
                                    r1%y(1:ni,1:nj,1:nk) = thedata(1:ni,1:nj,1:nk)
                              elseif(WORK_DATA == R2_DATA) then
                                    r2%xy(1:ni,1:nj,1:nk) = thedata(1:ni,1:nj,1:nk)
                              endif
                        case(3)
                              if(WORK_DATA == R0_DATA) then
                                    cycle
                              elseif(WORK_DATA == R1_DATA) then
                                    r1%z(1:ni,1:nj,1:nk) = thedata(1:ni,1:nj,1:nk)
                              elseif(WORK_DATA == R2_DATA) then
                                    r2%xz(1:ni,1:nj,1:nk) = thedata(1:ni,1:nj,1:nk)
                              endif
                        case(4)
                              r2%yx(1:ni,1:nj,1:nk) = thedata(1:ni,1:nj,1:nk)
                        case(5)
                              r2%yy(1:ni,1:nj,1:nk) = thedata(1:ni,1:nj,1:nk)
                        case(6)
                              r2%yz(1:ni,1:nj,1:nk) = thedata(1:ni,1:nj,1:nk)
                        case(7)
                              r2%zx(1:ni,1:nj,1:nk) = thedata(1:ni,1:nj,1:nk)
                        case(8)
                              r2%zy(1:ni,1:nj,1:nk) = thedata(1:ni,1:nj,1:nk)
                        case(9)
                              r2%zz(1:ni,1:nj,1:nk) = thedata(1:ni,1:nj,1:nk)
                        case default
                              cycle
                        end select
                  endif

            enddo

            deallocate(thedata)
            close(funit)

            t1 = cputimer(lcscomm,SYNC_TIMER)
            cpu_io = cpu_io + max(t1-t0,0.0)

      end subroutine structured_io

      subroutine unstructured_io(fname,IO_ACTION,r0,r1,r2)
      IMPLICIT NONE
            !-----
            character(len=*):: fname
            integer:: IO_ACTION
            type(ur0_t),optional:: r0
            type(ur1_t),optional:: r1
            type(ur2_t),optional:: r2
            !-----
            if(lcsrank==0) &
                        write(*,'(a,a,a)') 'In unstructured_io... NO ACTION (not supported):  Need to compile with HDF5 support'
            return

      end subroutine unstructured_io

end module io_m
