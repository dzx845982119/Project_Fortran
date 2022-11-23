!****************************************************
!     3次スプライン補間
!     inputファイル:
!     1行目:グリッドファイル名
!     2行目:入力ファイル名(流速のデータ)
!     3行目:出力ファイル名(補間値の出力データ)
!     4行目:補間したい点のx座標
!     以降:未定
!
!
!　　　簡単な説明：
!     スプライン補間を行うのに必要なのは複数のx,y座標および補間したい点のx座標である。
!     本プログラムは第329行からグリッドファイル、流速のデータ、補間したい点のx座標を読み取り、
!     その計算を第46-322行で行う。
!
!　　 ポイントファイル:
!     点の数:n
!     2グリッドファイル:(x1, y1)
!                      (x2, y2)  
!                         .
!                         .
!                         .
!                      (xn, yn)           
!
! 2021/09/14  Dong Zixu  v0.1
!
!****************************************************
module const
    ! DP: 単精度(4), DP: 倍精度(8)
    integer,     parameter :: SP = kind(1.0)
    integer(SP), parameter :: DP = selected_real_kind(2 * precision(1.0_SP))
  end module const
  
  module spline
    use const
    implicit none
    private
    public :: interpolate
  
  contains
    ! 3次スプライン補間
    !
    ! n: 点の数
    ! p(:, :): 点の配列
    ! aaa2: 補間するために入力された横座標の値
    ! out: 3次スプライン補間から求められた補間値
    subroutine interpolate(n, p, aaa2, out)
      use const
      implicit none
      integer(SP), intent(in) :: n
      real(SP),    intent(in) :: p(2, n)
      real(SP) aaa2, out
      real(SP) :: h(n - 1), w(n - 2), mtx(n - 1, n - 2), v(n)
      real(SP) :: a(n - 1), b(n - 1), c(n - 1), d(n)
      real(SP)    :: xp, xp_2, xp_3, x, y
      integer(SP) :: idx, i
  
      ! 初期計算
      h = calc_h(n, p)
      w = calc_w(n, p, h)
      call gen_mtx(n, h, w, mtx)
      call gauss_jordan(n, mtx, v)
      v = cshift(v, -1)
      b = calc_b(n, v)
      a = calc_a(n, v, p)
      d = calc_d(n, p)
      c = calc_c(n, v, p)
  
      ! 補間
        x = aaa2
        idx = find_idx(n, p, x)
        xp   = x - p(1, idx)
        xp_2 = xp * xp
        xp_3 = xp_2 * xp
        y = a(idx) * xp_3 + b(idx) * xp_2 + c(idx) * xp + d(idx)
        out = y
        write(*,*)'aaa2',x
        write(*,*)'idx',idx
        write(*,*)'xp',xp 
        write(*,*)'xp_2',xp_2
        write(*,*)'xp_3',xp_3 
        write(*,*)'out',out 
    end subroutine interpolate
  
    ! h 計算
    !
    ! :param(in) integer(4)     n: 点の数
    ! :param(in) real(8)  p(2, n): 点の配列
    ! :return    real(8) h(n - 1): h の配列
    function calc_h(n, p) result(h)
      implicit none
      integer(SP), intent(in) :: n
      real(SP),    intent(in) :: p(2, n)
      real(SP) :: h(n - 1)
  
      h = p(1, 2:n) - p(1, 1:n-1)
    end function calc_h
  
    ! w 計算
    !
    ! :param(in) integer(4)     n: 点の数
    ! :param(in) real(8)  p(2, n): 点の配列
    ! :param(in) real(8) h(n - 1): h の配列
    ! :return    real(8) w(n - 2): w の配列
    function calc_w(n, p, h) result(w)
      implicit none
      integer(SP), intent(in) :: n
      real(SP),    intent(in) :: p(2, n), h(n - 1)
      real(SP) :: w(n - 2)
  
      w = 6.0_DP * ((p(2, 3:n) - p(2, 2:n-1)) &
      & / h(2:n-1) - (p(2, 2:n-1) - p(2, 1:n-2)) / h(1:n-2))
    end function calc_w
  
    ! 行列生成
    !
    ! :param(in)  integer(4)             n: 点の数
    ! :param(in)  real(8)         h(n - 1): h の配列
    ! :param(in)  real(8)         w(n - 2): w の配列
    ! :param(out) real(8) mtx(n - 1, n - 2: 行列
    subroutine gen_mtx(n, h, w, mtx)
      implicit none
      integer(SP), intent(in)  :: n
      real(SP),    intent(in)  :: h(n - 1), w(n - 1)
      real(SP),    intent(out) :: mtx(n - 1, n - 2)
      integer(SP) :: i
  
      mtx(:, :) = 0.0_DP
      do i = 1, n - 2
        mtx(i, i)  = 2.0_DP * (h(i) + h(i + 1))
        mtx(n - 1, i) = w(i)
        if (i == 1) cycle
        mtx(i, i - 1) = h(i)
        mtx(i - 1, i) = h(i)
      end do
    end subroutine gen_mtx
  
    ! 連立一次方程式を解く（Gauss-Jordan 法）
    !
    ! :param(in)  integer(4)             n: 点の数
    ! :param(in)  real(8) mtx(n - 1, n - 2: 行列
    ! :param(out) real(8)             v(n): v の配列
    subroutine gauss_jordan(n, mtx, v)
      implicit none
      integer(SP), intent(in)  :: n
      real(SP),    intent(in)  :: mtx(n - 1, n - 2)
      real(SP),    intent(out) :: v(n)
      integer(SP) :: i, j
      real(SP)    :: mtx_w(n - 1, n - 2)  ! 作業用
      real(SP)    :: p, d
  
      mtx_w(:, :) = mtx(:, :)
      v(:) = 0.0_DP
      do j = 1, n - 2
        p = mtx_w(j, j)
        mtx_w(j:n-1, j) = mtx_w(j:n-1, j) / p
        do i = 1, n - 2
          if (i == j) cycle
          d = mtx_w(j, i)
          mtx_w(j:n-1, i) = mtx_w(j:n-1, i) - d * mtx_w(j:n-1, j)
        end do
      end do
      v(1:n-2) = mtx_w(n - 1, 1:n-2)
    end subroutine gauss_jordan
  
    ! a 計算
    !
    ! :param(in) integer(4)     n: 点の数
    ! :param(in) real(8)     v(n): v の配列
    ! :param(in) real(8)  p(2, n): 点の配列
    ! :return    real(8) a(n - 1): a の配列
    function calc_a(n, v, p) result(a)
      implicit none
      integer(SP), intent(in) :: n
      real(SP),    intent(in) :: v(n), p(2, n)
      real(SP) :: a(n - 1)
  
      a = (v(2:n) - v(1:n-1)) / (6.0_DP * (p(1, 2:n) - p(1, 1:n-1)))
    end function calc_a
  
    ! b 計算
    !
    ! :param(in) integer(4)     n: 点の数
    ! :param(in) real(8)     v(n): v の配列
    ! :return    real(8) b(n - 1): b の配列
    function calc_b(n, v) result(b)
      implicit none
      integer(SP), intent(in) :: n
      real(SP),    intent(in) :: v(n)
      real(SP) :: b(n - 1)
  
      b = v(1:n-1) / 2.0_DP
    end function calc_b
  
    ! c 計算
    !
    ! :param(in) integer(4)     n: 点の数
    ! :param(in) real(8)     v(n): v の配列
    ! :param(in) real(8)  p(2, n): 点の配列
    ! :return    real(8)     c(n): c の配列
    function calc_c(n, v, p) result(c)
      implicit none
      integer(SP), intent(in) :: n
      real(SP),    intent(in) :: v(n), p(2, n)
      real(SP) :: c(n - 1)
  
      c = (p(2, 2:n) - p(2, 1:n-1)) / (p(1, 2:n) - p(1, 1:n-1)) &
      & - (p(1, 2:n) - p(1, 1:n-1)) * (2.0_DP * v(1:n-1) + v(2:n)) &
      & / 6.0_DP
    end function calc_c
  
    ! d 計算
    !
    ! :param(in) integer(4)     n: 点の数
    ! :param(in) real(8)  p(2, n): 点の配列
    ! :return    real(8)     d(n): d の配列
    function calc_d(n, p) result(d)
      implicit none
      integer(SP), intent(in) :: n
      real(SP),    intent(in) :: p(2, n)
      real(SP) :: d(n)
  
      d = p(2, :)
    end function calc_d
  
    ! インデックス検索
    !
    ! :param(in) integer(4)     n: 点の数
    ! :param(in) real(8)  p(2, n): 点の配列
    ! :param(in) real(8)        x: x 値
    ! :return    integer(4)   idx: インデックス
    function find_idx(n, p, x) result(idx)
      implicit none
      integer(SP), intent(in) :: n
      real(SP),    intent(in) :: p(2, n), x
      integer(SP) :: idx
      integer(SP) :: i, j, k
  
      i = 1
      j = n
      do while (i < j)
        k = (i + j) / 2
        if (p(1, k) < x) then
          i = k + 1
        else
          j = k
        end if
      end do
      if (i > 1) i = i - 1
      idx = i
    end function find_idx
  end module spline
  
  subroutine spline_interpolation(n,aaa2,aaa,uuu,out)
    use const
    use spline
    implicit none
    integer(SP) :: n 
    real(SP) aaa2,out
    real(SP) aaa(n),uuu(n)
    real(SP), allocatable :: p(:, :)  ! 点の配列
    integer(SP) :: i, ios
  

  
    ! 点の配列メモリ確保
    allocate(p(2, n))
  
    ! 点の配列読み込み
    do i = 1, n
      p(1, i)=aaa(i)
      p(2, i)=uuu(i)
    end do

    do i=1,n
    !write(*,*)'p(1, i),p(2, i)=',p(1, i),p(2, i)
    end do
  
    ! 補間
    call interpolate(n, p, aaa2, out)
    write(*,*)'aaa2,out',aaa2,out
  
    ! 点の配列メモリ解放
    deallocate(p)

  end subroutine spline_interpolation

    program main !///dzx
        use const
        use spline
        !/////////////  v-out  ////////////dzx
        integer*4 nn,mm,ll
        parameter(nn=2073,mm=171,ll=82)
        real*8 xx(nn,mm,ll),yy(nn,mm,ll),zz(nn,mm,ll),re0(4)
        real*8, dimension(0:nn+1,0:mm+1,0:ll+1) :: pp,uu,vv,ww
        !real*8 uu(nn,mm,ll),vv(nn,mm,ll),ww(nn,mm,ll),pp(nn,mm,ll),re0(4)
        real*4 uuu(nn-1),vvv(nn-1)
        real*4 aa(nn-1),aaa(nn-1)
  
        integer it0,ii,jjj,iii,i_stat
        integer n0,m0,l0
        integer i,j,k,l
        integer ip
        real*4 pnt,rod,rin,dyy,dxx,dz
        real*4 dA,dB,uout,vout,aaa2
  
        ! 4 bites !
        real*8 time0,vin0,timer
        character*32 rdfile,wtfile,gdfile
        !  ### input data ###
          print *,'Input GRID file name?'
          read(5,'(a32)')gdfile
          write(*,*)gdfile
  
          print *,'Input FLOW file name?'
          read(5,'(a32)')rdfile
          write(*,*)rdfile
  
          print *,'Output point file name?'
          read(5,'(a32)')wtfile
          write(*,*)wtfile
  
          timer=0.d0
  
  ! #####################
  ! ##### READ DATA #####
  ! #####################
  
          ! ### read GRID data ###
          
          open(10,file=gdfile,status='unknown',form='unformatted')
              read(10) n0,m0,l0,xx,yy,zz
              !read(10) (((xx(j,k,l),j=1,n0),k=1,m0),l=1,l0)
              !read(10) (((yy(j,k,l),j=1,n0),k=1,m0),l=1,l0)
              !read(10) (((zz(j,k,l),j=1,n0),k=1,m0),l=1,l0)
          close(10)
  
          pnt=abs(xx(1,1,1)) 
          rod=xx(n0,1,1)  
          rin=yy((n0+1)/2,1,1)  
  
          print *,'### Read GRID data; n0,m0,l0 =',n0,m0,l0    
          print *,'pnt=',real(pnt)
          print *,'rod=',real(rod)-0.5
          print *,'rin=',real(rin)
          ! ##### FLOW data #####
  
          time0=0.
          i_stat=0
          uout=0
          vout=0
          ! read data from file(s)
          open(20,file=rdfile,status='old',form='unformatted')
          read(20)iii,it0,time0,re0,uu,vv,ww,pp
          print *,'### Read FLOW data; iii,it0,time0,re0'
          print *,iii,it0,time0,re0
          !read(20,iostat=i_stat)jjj,(((uu(j,k,l),j=1,n0),k=1,m0),l=1,l0),jjj
          !read(20,iostat=i_stat)jjj,(((vv(j,k,l),j=1,n0),k=1,m0),l=1,l0),jjj
          !read(20,iostat=i_stat)jjj,(((ww(j,k,l),j=1,n0),k=1,m0),l=1,l0),jjj
          close(20)
          write(*,*)'vin=',vin0
          write(*,*)' '

  ! #####################
  ! ##### MAIN LOOP #####
  ! #####################
  
        print *,'Input xp=?'
        read(5,*) xp
        print *,'xp=',real(xp)
        uout=0
        vout=0
  
      open(unit=30,file=wtfile,status='unknown',form='formatted')
      l = 1
      do j = 1,m0  
        do i = 1,n0-1
          !-------J AT P-------
                  dA=xx(i,j,l)-xp
                  dB=xx(i+1,j,l)-xp
                  if(dA*dB<=0) then
                      ip=i
                  end if
        end do
  
        do i= 1,n0
          uuu(i)=uu(i,j,l)
          vvv(i)=vv(i,j,l)
          !write(*,*)'uuu(i)=u(i,j,l)=',uuu(i),uu(i,j,l)
        end do
  
        write(*,*)'j,ip=',j,ip
  
        aa = 0
        aaa= 0
          do i = 1,n0-1
            dx2=xx(i+1,j,l)-xx(i,j,l)
            dy2=yy(i+1,j,l)-yy(i,j,l)
            aa(i)=sqrt(dx2*dx2+dy2*dy2)
            !write(*,*)'aa(',i,')',aa(i)
          end do
  
          aaa(1)=aa(1)
          do i = 2,n0-1
            aaa(i)=aaa(i-1)+aa(i)
            !write(*,*)'aaa(',i,')',aaa(i)
          end do
  
          dyy =yy(ip+1,j,l)-yy(ip,j,l)
          dxx =xx(ip+1,j,l)-xx(ip,j,l)
          dx  =xp-xx(ip,j,l)
          dy  =(dyy/dxx)*dx
          dz  =sqrt(dx*dx+dy*dy)
          aaa2=dz+aaa(ip)
          write(*,*)'dyy,dxx,dy,dx,dz,aaa2',dyy,dxx,dy,dx,dz,aaa2  
          call spline_interpolation (nn-1,aaa2,aaa,uuu,uout)
          call spline_interpolation (nn-1,aaa2,aaa,vvv,vout)

          write(*,*)'out put --->',j,yy(ip,j,l)+dy,uout,vout
          write(30,*)j,yy(ip,j,l)+dy,uout,vout

        end do
        close(30)
  end program main