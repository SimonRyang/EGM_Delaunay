module delaunay_func

use comtypes

! these are a number of functions and subroutines copied from geompack3 http://people.sc.fsu.edu/~jburkardt/f_src/geompack3/geompack3.html

implicit none

contains

!*****************************************************************************
subroutine dtris2 ( npt, maxst, vcl, ind, ntri, til, tnbr, stack, ierr )
!
!! DTRIS2 constructs the Delaunay triangulation of vertices in 2D.
!
!  Discussion: 
!
!    This routine constructs a Delaunay triangulation of 2D vertices using
!    incremental approach and diagonal edge swaps.  Vertices are
!    first sorted in lexicographically increasing (x,y) order, and
!    then are inserted one at a time from outside the convex hull.
!
!  Author:
!
!    Barry Joe, 
!    Department of Computing Science, 
!    University of Alberta,
!    Edmonton, Alberta, Canada  T6G 2H1
!    Email: barry@cs.ualberta.ca
!
!  Reference:
!
!    Barry Joe,
!    GEOMPACK - a software package for the generation of meshes
!      using geometric algorithms, 
!    Advances in Engineering Software,
!    Volume 13, pages 325-331, 1991.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NPT, the number of 2D points (vertices).
!
!    Input, integer ( kind = 4 ) MAXST, the maximum size available for STACK array; 
!    should be about NPT to be safe, but MAX ( 10, 2*LOG2(NPT) ) usually enough.
!
!    Input, real ( kind = 8 ) VCL(1:2,1:*), coordinates of 2D vertices.
!
!    Input/output, IND(1:NPT), the indices in VCL of vertices to be
!    triangulated.  On output, permuted due to sorting.
!
!    Output, integer ( kind = 4 ) NTRI, the number of triangles in triangulation; equal to
!    2*NPT - NB - 2 where NB = number of boundary vertices.
!
!    Output, integer ( kind = 4 ) TIL(1:3,1:NTRI), the triangle incidence list; elements 
!    are indices of VCL; vertices of triangles are in counterclockwise order.
!
!    Output, integer ( kind = 4 ) TNBR(1:3,1:NTRI), the triangle neighbor list; positive
!    elements are indices of TIL; negative elements are used for links
!    of counterclockwise linked list of boundary edges; 
!    LINK = -(3*I + J-1) where I, J = triangle, edge index; 
!    TNBR(J,I) refers to the neighbor along edge from vertex 
!    J to J+1 (mod 3).
!
!    Workspace, integer STACK(1:MAXST), used for stack of triangles for which
!    circumcircle test must be made.
!
!    Output, integer ( kind = 4 ) IERR, error flag, which is zero unless an error occurred.
!

  integer ( kind = 4 ) maxst
  integer ( kind = 4 ) npt

  real    ( prec ) cmax
  integer ( kind = 4 ) e
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierr
  integer ( kind = 4 ) ind(npt)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) l
  integer ( kind = 4 ) ledg
  integer ( kind = 4 ) lr
  integer ( kind = 4 ) lrline
  integer ( kind = 4 ) ltri
  integer ( kind = 4 ) m
  integer ( kind = 4 ) m1
  integer ( kind = 4 ) m2
  integer ( kind = 4 ), parameter :: msglvl = 0
  integer ( kind = 4 ) n
  integer ( kind = 4 ) ntri
  integer ( kind = 4 ) redg
  integer ( kind = 4 ) rtri
  integer ( kind = 4 ) stack(maxst)
  integer ( kind = 4 ) t
  integer ( kind = 4 ) til(3,npt*2)
  integer ( kind = 4 ) tnbr(3,npt*2)
  real    ( prec ) tol
  integer ( kind = 4 ) top
  real    ( prec ) vcl(2,*)

  ierr = 0
  tol = 100.0D+00 * epsilon ( tol )
!
!  Sort vertices by increasing (x,y).
!
  call dhpsrt ( 2, npt, 2, vcl, ind )
!
!  Check that the data is not degenerate.
!
  m1 = ind(1)

  do i = 2, npt

    m = m1
    m1 = ind(i)

    do j = 1, 2
      cmax = max ( abs(vcl(j,m)), abs(vcl(j,m1)) )
      if ( tol * cmax < abs(vcl(j,m) - vcl(j,m1)) .and. tol < cmax ) then
        go to 20
      end if
    end do

    ierr = 224
    return

20  continue

  end do
!
!  Staring with points M1 and M2, find the first point M that is
!  "reasonably" non-colinear.
!
  m1 = ind(1)
  m2 = ind(2)
  j = 3

  do

    if ( npt < j ) then
      ierr = 225
      return
    end if

    m = ind(j)
    lr = lrline ( vcl(1,m), vcl(2,m), vcl(1,m1), vcl(2,m1), vcl(1,m2), &
      vcl(2,m2), 0.0D+00 )

    if ( lr /= 0 ) then
      exit
    end if

    j = j + 1

  end do
!
!  Set up the initial triangle information for M1, M2 and M (and any
!  in-between points we may have skipped over while searching for M.
!
  ntri = j - 2

  if ( lr == -1 ) then

    til(1,1) = m1
    til(2,1) = m2
    til(3,1) = m
    tnbr(3,1) = -3

    do i = 2, ntri
      m1 = m2
      m2 = ind(i+1)
      til(1,i) = m1
      til(2,i) = m2
      til(3,i) = m
      tnbr(1,i-1) = -3*i
      tnbr(2,i-1) = i
      tnbr(3,i) = i - 1
    end do

    tnbr(1,ntri) = -3*ntri - 1
    tnbr(2,ntri) = -5
    ledg = 2
    ltri = ntri

  else

    til(1,1) = m2
    til(2,1) = m1
    til(3,1) = m
    tnbr(1,1) = -4

    do i = 2, ntri
      m1 = m2
      m2 = ind(i+1)
      til(1,i) = m2
      til(2,i) = m1
      til(3,i) = m
      tnbr(3,i-1) = i
      tnbr(1,i) = -3*i - 3
      tnbr(2,i) = i - 1
    end do

    tnbr(3,ntri) = -3*ntri
    tnbr(2,1) = -3*ntri - 2
    ledg = 2
    ltri = 1

  end if

  if ( msglvl == 4 ) then
    m2 = ind(1)
    write ( *,'(i7,4f15.7)') 1,vcl(1,m2),vcl(2,m2),vcl(1,m),vcl(2,m)
    do i = 2, j-1
      m1 = m2
      m2 = ind(i)
      write ( *,'(i7,4f15.7)') 1,vcl(1,m1),vcl(2,m1),vcl(1,m2),vcl(2,m2)
      write ( *,'(i7,4f15.7)') 1,vcl(1,m2),vcl(2,m2),vcl(1,m),vcl(2,m)
    end do
  end if
!
!  Insert vertices one at a time from outside convex hull, determine
!  visible boundary edges, and apply diagonal edge swaps until
!  Delaunay triangulation of vertices (so far) is obtained.
!
  top = 0

  do i = j+1, npt

    if ( msglvl == 4 ) then
      write ( *,'(i7,4f15.7)') i
    end if

    m = ind(i)
    m1 = til(ledg,ltri)

    if ( ledg <= 2 ) then
      m2 = til(ledg+1,ltri)
    else
      m2 = til(1,ltri)
    end if

    lr = lrline(vcl(1,m),vcl(2,m),vcl(1,m1),vcl(2,m1),vcl(1,m2), &
      vcl(2,m2),0.0D+00)

    if ( 0 < lr ) then
      rtri = ltri
      redg = ledg
      ltri = 0
    else
      l = -tnbr(ledg,ltri)
      rtri = l / 3
      redg = mod(l,3) + 1
    end if

    call vbedg(vcl(1,m),vcl(2,m),vcl,til,tnbr,ltri,ledg,rtri,redg)
    n = ntri + 1
    l = -tnbr(ledg,ltri)

    do

      t = l / 3
      e = mod(l,3) + 1
      l = -tnbr(e,t)
      m2 = til(e,t)

      if ( e <= 2 ) then
        m1 = til(e+1,t)
      else
        m1 = til(1,t)
      end if

      ntri = ntri + 1
      tnbr(e,t) = ntri
      til(1,ntri) = m1
      til(2,ntri) = m2
      til(3,ntri) = m
      tnbr(1,ntri) = t
      tnbr(2,ntri) = ntri - 1
      tnbr(3,ntri) = ntri + 1
      top = top + 1

      if ( maxst < top ) then
        ierr = 8
        return
      end if

      stack(top) = ntri

      if ( msglvl == 4 ) then
        write ( *,'(i7,4f15.7)') 1,vcl(1,m),vcl(2,m), vcl(1,m2),vcl(2,m2)
      end if

      if ( t == rtri .and. e == redg ) then
        exit
      end if

    end do

    if ( msglvl == 4 ) then
      write ( *,'(i7,4f15.7)') 1,vcl(1,m),vcl(2,m), vcl(1,m1),vcl(2,m1)
    end if

    tnbr(ledg,ltri) = -3*n - 1
    tnbr(2,n) = -3*ntri - 2
    tnbr(3,ntri) = -l
    ltri = n
    ledg = 2
    call swapec(m,top,maxst,ltri,ledg,vcl,til,tnbr,stack, ierr )

    if ( ierr /= 0 ) then
      return
    end if

  end do

  if ( msglvl == 4 ) then
    write ( *, '(i7,4f15.7)' ) npt+1
  end if

  return
end subroutine dtris2
!*****************************************************************************

!*****************************************************************************80
subroutine walkt2 ( x, y, ntri, vcl, til, tnbr, itri, iedg, ierr )
!
!! WALKT2 walks through a 2D triangulation searching for a point.
!
!  Discussion: 
!
!    This routine walks through neighboring triangles of a 2D (Delaunay)
!    triangulation until a triangle is found containing point (X,Y)
!    or (X,Y) is found to be outside the convex hull.  The search is
!    guaranteed to terminate for a Delaunay triangulation, else a
!    cycle may occur.
!
!  Author:
!
!    Barry Joe, 
!    Department of Computing Science, 
!    University of Alberta,
!    Edmonton, Alberta, Canada  T6G 2H1
!    Email: barry@cs.ualberta.ca
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, Y, the 2D point.
!
!    Input, integer ( kind = 4 ) NTRI, the number of triangles in triangulation; 
!    used to detect cycle.
!
!    Input, real ( kind = 8 ) VCL(1:2,1:*), the coordinates of 2D vertices.
!
!    Input, integer ( kind = 4 ) TIL(1:3,1:*), the triangle incidence list.
!
!    Input, integer ( kind = 4 ) TNBR(1:3,1:*), the triangle neighbor list.
!
!    Input/output, integer ( kind = 4 ) ITRI.  On input, the index of triangle to begin
!    search at.  On output, the index of triangle that search ends at.
!
!    Output, integer ( kind = 4 ) IEDG, 0 if ( X,Y) is in the interior of triangle ITRI; 
!    I = 1, 2, or 3 if ( X,Y) is on interior of edge I of ITRI;
!    I = 4, 5, or 6 if ( X,Y) is (nearly) vertex I-3 of ITRI;
!    I = -1, -2, or -3 if ( X,Y) is outside convex hull due
!    to walking past edge -I of triangle ITRI.
!
!    Output, integer ( kind = 4 ) IERR, error flag, which is zero unless an error occurred.
use comparam

  integer ( kind = 4 ) a
  real    ( prec ) alfa
  integer ( kind = 4 ) b
  real    ( prec ) beta
  integer ( kind = 4 ) c
  integer ( kind = 4 ) cnt
  real    ( prec ) det
  real    ( prec ) dx
  real    ( prec ) dxa
  real    ( prec ) dxb
  real    ( prec ) dy
  real    ( prec ) dya
  real    ( prec ) dyb
  real    ( prec ) gama
  integer ( kind = 4 ) i
  integer ( kind = 4 ) iedg
  integer ( kind = 4 ) ierr
  integer ( kind = 4 ) itri
  integer ( kind = 4 ) ntri
  integer ( kind = 4 ) til(3,2*npt)
  integer ( kind = 4 ) tnbr(3,2*npt)
  real    ( prec ) tol
  real    ( prec ) vcl(2,npt)
  real    ( prec ) x
  real    ( prec ) y
  
  integer ( kind = 4 ) count
!
!  Use barycentric coordinates to determine where (X,Y) is located.
!
  ierr = 0
  tol = 100.0D+00 * epsilon ( tol )
  cnt = 0
  iedg = 0
  count=0

10 continue

  cnt = cnt + 1

  if ( ntri+1 < cnt ) then
    ierr = 226
    return
  end if

  a = til(1,itri)
  b = til(2,itri)
  c = til(3,itri)
  dxa = vcl(1,a) - vcl(1,c)
  dya = vcl(2,a) - vcl(2,c)
  dxb = vcl(1,b) - vcl(1,c)
  dyb = vcl(2,b) - vcl(2,c)
  dx = x - vcl(1,c)
  dy = y - vcl(2,c)
  det = dxa*dyb - dya*dxb
  alfa = (dx*dyb - dy*dxb) / det
  beta = (dxa*dy - dya*dx) / det
  gama = 1.0D+00 - alfa - beta
  count = count +1
  if ( tol < alfa .and. tol < beta .and. tol < gama ) then

    return

  else if ( alfa < -tol ) then

    i = tnbr(2,itri)
    if ( i <= 0 ) then
      iedg = -2
      return
    end if

  else if ( beta < -tol ) then

    i = tnbr(3,itri)
    if ( i <= 0 ) then
      iedg = -3
      return
    end if

  else if ( gama < -tol ) then

    i = tnbr(1,itri)
    if ( i <= 0 ) then
      iedg = -1
      return
    end if

  else if ( alfa <= tol ) then

    if ( beta <= tol ) then
      iedg = 6
    else if ( gama <= tol ) then
      iedg = 5
    else
      iedg = 2
    end if
    return

  else if ( beta <= tol ) then

    if ( gama <= tol ) then
      iedg = 4
    else
      iedg = 3
    end if
    return

  else

    iedg = 1
    return

  end if

  itri = i
  go to 10

end subroutine walkt2
!*****************************************************************************

!*****************************************************************************80
subroutine dhpsrt ( k, n, lda, a, map )
!
!! DHPSRT sorts a list of double precision points in KD.
!
!  Discussion: 
!
!    This routine uses heapsort to obtain the permutation of N K-dimensional
!    double precision points so that the points are in lexicographic
!    increasing order.
!
!  Modified:
!
!    31 August 2005
!
!  Author:
!
!    Barry Joe, 
!    Department of Computing Science, 
!    University of Alberta,
!    Edmonton, Alberta, Canada  T6G 2H1
!    Email: barry@cs.ualberta.ca
!
!  Reference:
!
!    Barry Joe,
!    GEOMPACK - a software package for the generation of meshes
!      using geometric algorithms, 
!    Advances in Engineering Software,
!    Volume 13, pages 325-331, 1991.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) K, the dimension of points.
!
!    Input, integer ( kind = 4 ) N, the number of points.
!
!    Input, integer ( kind = 4 ) LDA, the leading dimension of array A in calling 
!    routine; K <= LDA.
!
!    Input, real ( kind = 8 ) A(1:K,1:*), array of points.
!
!    Input/output, integer ( kind = 4 ) MAP(1:N).  On input, he points of A with indices 
!    MAP(1), MAP(2), ..., MAP(N) are to be sorted.  On output, the elements 
!    are permuted so that A(*,MAP(1)) <= A(*,MAP(2)) <= ... <= A(*,MAP(N)).
!

  integer ( kind = 4 ) lda
  integer ( kind = 4 ) n

  real    ( prec ) a(lda,*)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) k
  integer ( kind = 4 ) map(n)
  integer ( kind = 4 ) t

  do i = n/2, 1, -1
    call dsftdw ( i, n, k, lda, a, map )
  end do

  do i = n, 2, -1
    t = map(1)
    map(1) = map(i)
    map(i) = t
    call dsftdw ( 1, i-1, k, lda, a, map )
  end do

  return
end subroutine dhpsrt
!*****************************************************************************





!*****************************************************************************80
subroutine dsftdw ( l, u, k, lda, a, map )
!
!! DSFTDW does one step of the heap sort algorithm for double precision data.
!
!  Discussion: 
!
!    This routine sifts A(*,MAP(L)) down a heap of size U.
!
!  Author:
!
!    Barry Joe, 
!    Department of Computing Science, 
!    University of Alberta,
!    Edmonton, Alberta, Canada  T6G 2H1
!    Email: barry@cs.ualberta.ca
!
!  Reference:
!
!    Barry Joe,
!    GEOMPACK - a software package for the generation of meshes
!      using geometric algorithms, 
!    Advances in Engineering Software,
!    Volume 13, pages 325-331, 1991.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) L, U, the lower and upper index of part of heap.
!
!    Input, integer ( kind = 4 ) K, the dimension of points.
!
!    Input, integer ( kind = 4 ) LDA, the leading dimension of array A in calling routine.
!
!    Input, real ( kind = 8 ) A(1:K,1:*), see routine DHPSRT.
!
!    Input/output, integer ( kind = 4 ) MAP(1:*), see routine DHPSRT.
!


  integer ( kind = 4 ) lda

  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  integer ( kind = 4 ) map(*)
  integer ( kind = 4 ) u

  real    ( prec ) a(lda,*)
  logical              dless
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) t

  i = l
  j = 2 * i
  t = map(i)

  do

    if ( u < j ) then
      exit
    end if

    if ( j < u ) then
      if ( dless ( k, a(1,map(j)), a(1,map(j+1)) ) ) then
        j = j + 1
      end if
    end if

    if ( dless ( k, a(1,map(j)), a(1,t)) ) then
      exit
    end if

    map(i) = map(j)
    i = j
    j = 2 * i

  end do

  map(i) = t

  return
end subroutine dsftdw
!*****************************************************************************

!*****************************************************************************80
function htsrck ( k, ind, n, p, fc, ht )
!
!! HTSRCK searches for a record in the hash table.
!
!  Discussion: 
!
!    This routine searches for record FC(1:K+4,POS) containing key IND(1:K)
!    in hash table HT.
!
!  Author:
!
!    Barry Joe, 
!    Department of Computing Science, 
!    University of Alberta,
!    Edmonton, Alberta, Canada  T6G 2H1
!    Email: barry@cs.ualberta.ca
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) K, number of vertices in a face.
!
!    Input/output, integer ( kind = 4 ) IND(1:K), vertex indices of face.  On output, these
!    have been sorted into nondecreasing order.
!
!    Input, integer ( kind = 4 ) N, upper bound on vertex indices.
!
!    Input, integer ( kind = 4 ) P, size of hash table.
!
!    Input, integer ( kind = 4 ) FC(1:K+4,1:*), array of face records; see routine DTRISK.
!
!    Input, integer ( kind = 4 ) HT(0:P-1), hash table using direct chaining.
!
!    Output, integer ( kind = 4 ) HTSRCK, position of FC record with key IND(1:K) if found,
!    or 0 if not found.
!


  integer ( kind = 4 ) k
  integer ( kind = 4 ) p

  integer ( kind = 4 ) fc(k+4,*)
  integer ( kind = 4 ) h
  integer ( kind = 4 ) ht(0:p-1)
  integer ( kind = 4 ) htsrck
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ind(k)
  integer ( kind = 4 ) kp3
  integer ( kind = 4 ) n
  integer ( kind = 4 ) pos

  kp3 = k + 3
  call orderk ( k, ind )

  h = ind(1)
  do i = 2, k
    h = mod ( h * n + ind(i), p )
  end do

  pos = ht(h)

20 continue

  if ( pos /= 0 ) then

    i = 1

30  continue

    if ( fc(i,pos) /= ind(i) ) then
      pos = fc(kp3,pos)
      go to 20
    end if

    i = i + 1

    if ( i <= k ) then
      go to 30
    end if

  end if

  htsrck = pos

  return
end function htsrck
!*****************************************************************************



!*****************************************************************************80
subroutine orderk ( k, ind )
!
!! ORDERK reorders K elements of an array in nondecreasing order.
!
!  Discussion: 
!
!    This routine reorders K elements of array IND in nondecreasing order.
!
!    It is assumed that K is small, say <= 15, so that insertion sort
!    is used. If K is larger, a faster sort such as heapsort should
!    be used.
!
!  Author:
!
!    Barry Joe, 
!    Department of Computing Science, 
!    University of Alberta,
!    Edmonton, Alberta, Canada  T6G 2H1
!    Email: barry@cs.ualberta.ca
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) K, the size of array IND.
!
!    Input/output, integer ( kind = 4 ) IND(1:K), an array, which is sorted on output.
!


  integer ( kind = 4 ) k

  integer ( kind = 4 ) i
  integer ( kind = 4 ) ind(k)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) s
  integer ( kind = 4 ) t

  do i = 2, k

    t = ind(i)
    j = i

    do

      s = ind(j-1)

      if ( s <= t ) then
        exit
      end if

      ind(j) = s
      j = j - 1
 
      if ( j <= 1 ) then
        exit
      end if

    end do

    ind(j) = t

  end do

  return
end subroutine orderk
!*****************************************************************************

!*****************************************************************************80
subroutine swapec ( i, top, maxst, btri, bedg, vcl, til, tnbr, stack, ierr )
!
!! SWAPEC swaps diagonal edges in a 2D triangulation
!
!  Discussion: 
!
!    This routine swaps diagonal edges in a 2D triangulation based on empty
!    circumcircle criterion until all triangles are Delaunay, given
!    that I is index of new vertex added to triangulation.
!
!  Author:
!
!    Barry Joe, 
!    Department of Computing Science, 
!    University of Alberta,
!    Edmonton, Alberta, Canada  T6G 2H1
!    Email: barry@cs.ualberta.ca
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) I, the index in VCL of new vertex.
!
!    Input/output, integer ( kind = 4 ) TOP, the index of top of stack, must be greater
!    than or equal to 0.
!
!    Input, integer ( kind = 4 ) MAXST, the maximum size available for STACK array.
!
!    Input/output, integer ( kind = 4 ) BTRI, BEDG, if positive, these are triangle and
!    edge index of a boundary edge whose updated indices must be recorded.
!
!    Input, real ( kind = 8 ) VCL(1:2,1:*), the coordinates of 2D vertices.
!
!    Input/output, integer ( kind = 4 ) TIL(1:3,1:*), the triangle incidence list.
!
!    Input/output, integer ( kind = 4 ) TNBR(1:3,1:*), the triangle neighbor list; negative
!    values are used for links of counterclockwise linked list of boundary
!    edges; LINK = -(3*I + J-1) where I, J = triangle, edge index.
!
!    Input, integer ( kind = 4 ) STACK(1:TOP), the index of initial triangles (involving
!    vertex I) put in stack; the edges opposite I should be in interior.
!
!    Workspace, integer STACK(TOP+1:MAXST), the used as stack.
!
!    Output, integer ( kind = 4 ) IERR, error flag, which is zero unless an error occurred.
!


  integer ( kind = 4 ) maxst

  integer ( kind = 4 ) a
  integer ( kind = 4 ) b
  integer ( kind = 4 ) bedg
  integer ( kind = 4 ) btri
  integer ( kind = 4 ) c
  integer ( kind = 4 ) diaedg
  integer ( kind = 4 ) e
  integer ( kind = 4 ) ee
  integer ( kind = 4 ) em1
  integer ( kind = 4 ) ep1
  integer ( kind = 4 ) f
  integer ( kind = 4 ) fm1
  integer ( kind = 4 ) fp1
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierr
  integer ( kind = 4 ) l
  integer ( kind = 4 ), parameter :: msglvl = 0
  integer ( kind = 4 ) r
  integer ( kind = 4 ) s
  integer ( kind = 4 ) stack(maxst)
  integer ( kind = 4 ) swap
  integer ( kind = 4 ) t
  integer ( kind = 4 ) til(3,*)
  integer ( kind = 4 ) tnbr(3,*)
  integer ( kind = 4 ) top
  integer ( kind = 4 ) tt
  integer ( kind = 4 ) u
  real    ( prec ) x
  real    ( prec ) y
  real    ( prec ) vcl(2,*)
!
!  Determine whether triangles in stack are Delaunay, and swap
!  diagonal edge of convex quadrilateral if not.
!
  ierr = 0
  x = vcl(1,i)
  y = vcl(2,i)

10 continue

  if ( top <= 0 ) then
    return
  end if

  t = stack(top)
  top = top - 1

  if ( til(1,t) == i ) then
    e = 2
    b = til(3,t)
  else if ( til(2,t) == i ) then
    e = 3
    b = til(1,t)
  else
    e = 1
    b = til(2,t)
  end if

  a = til(e,t)
  u = tnbr(e,t)

  if ( tnbr(1,u) == t ) then
    f = 1
    c = til(3,u)
  else if ( tnbr(2,u) == t ) then
    f = 2
    c = til(1,u)
  else
    f = 3
    c = til(2,u)
  end if

  swap = diaedg(x,y,vcl(1,a),vcl(2,a),vcl(1,c),vcl(2,c),vcl(1,b), &
    vcl(2,b))

  if ( swap == 1 ) then

    em1 = e - 1
    if ( em1 == 0) em1 = 3
    ep1 = e + 1
    if ( ep1 == 4) ep1 = 1
    fm1 = f - 1
    if ( fm1 == 0) fm1 = 3
    fp1 = f + 1
    if ( fp1 == 4) fp1 = 1
    til(ep1,t) = c
    til(fp1,u) = i
    r = tnbr(ep1,t)
    s = tnbr(fp1,u)
    tnbr(ep1,t) = u
    tnbr(fp1,u) = t
    tnbr(e,t) = s
    tnbr(f,u) = r

    if ( 0 < tnbr(fm1,u) ) then
      top = top + 1
      stack(top) = u
    end if

    if ( 0 < s ) then

      if ( tnbr(1,s) == u ) then
        tnbr(1,s) = t
      else if ( tnbr(2,s) == u ) then
        tnbr(2,s) = t
      else
        tnbr(3,s) = t
      end if

      top = top + 1

      if ( maxst < top ) then
        ierr = 8
        return
      end if

      stack(top) = t

    else

      if ( u == btri .and. fp1 == bedg ) then
        btri = t
        bedg = e
      end if

      l = -( 3 * t + e - 1 )
      tt = t
      ee = em1

20    continue

      if ( 0 < tnbr(ee,tt) ) then

        tt = tnbr(ee,tt)

        if ( til(1,tt) == a ) then
          ee = 3
        else if ( til(2,tt) == a ) then
          ee = 1
        else
          ee = 2
        end if

        go to 20

      end if

      tnbr(ee,tt) = l

    end if

    if ( 0 < r ) then

      if ( tnbr(1,r) == t ) then
        tnbr(1,r) = u
      else if ( tnbr(2,r) == t ) then
        tnbr(2,r) = u
      else
        tnbr(3,r) = u
      end if

    else

      if ( t == btri .and. ep1 == bedg ) then
        btri = u
        bedg = f
      end if

      l = -(3*u + f-1)
      tt = u
      ee = fm1

30    continue

      if ( 0 < tnbr(ee,tt) ) then

        tt = tnbr(ee,tt)

        if ( til(1,tt) == b ) then
          ee = 3
        else if ( til(2,tt) == b ) then
          ee = 1
        else
          ee = 2
        end if

        go to 30

      end if

      tnbr(ee,tt) = l

    end if

    if ( msglvl == 4 ) then
      write ( *,600) 2,vcl(1,a),vcl(2,a), &
           vcl(1,b),vcl(2,b),x,y,vcl(1,c),vcl(2,c)
    end if

  end if

  go to 10

  600 format (1x,i7,4f15.7/8x,4f15.7)

end subroutine swapec 
!*****************************************************************************

!*****************************************************************************80
subroutine vbedg ( x, y, vcl, til, tnbr, ltri, ledg, rtri, redg )
!
!! VBEDG determines the boundary edges of a 2D triangulation.
!
!  Discussion: 
!
!    This routine determines boundary edges of a 2D triangulation which are
!    visible from point (X,Y) outside convex hull.
!
!  Author:
!
!    Barry Joe, 
!    Department of Computing Science, 
!    University of Alberta,
!    Edmonton, Alberta, Canada  T6G 2H1
!    Email: barry@cs.ualberta.ca
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, Y, the 2D point outside convex hull.
!
!    Input, real ( kind = 8 ) VCL(1:2,1:*), the coordinates of 2D vertices.
!
!    Input, integer ( kind = 4 ) TIL(1:3,1:*), the triangle incidence list.
!
!    Input, integer ( kind = 4 ) TNBR(1:3,1:*), the triangle neighbor list; negative 
!    values are used for links of counterclockwise linked list of boundary
!    edges; LINK = -(3*I + J-1) where I, J = triangle, edge index.
!
!    Input, integer ( kind = 4 ) LTRI, LEDG, if LTRI /= 0 then they are assumed to be as
!    defined below and are not changed, else they are updated.  LTRI is
!    the index of boundary triangle to left of leftmost boundary
!    triangle visible from (X,Y).  LEDG is the boundary edge of triangle
!    LTRI to left of leftmost boundary edge visible from (X,Y)
!
!    Input/output, integer ( kind = 4 ) RTRI.  On input, the index of boundary triangle 
!    to begin search at.  On output, the index of rightmost boundary triangle 
!    visible from (X,Y)
!
!    Input/output, integer ( kind = 4 ) REDG.  On input, the edge of triangle RTRI that is 
!    visible from (X,Y).  On output, the edge of triangle RTRI that is 
!    visible from (X,Y)
!


  integer ( kind = 4 ) a
  integer ( kind = 4 ) b
  integer ( kind = 4 ) e
  integer ( kind = 4 ) l
  logical              ldone
  integer ( kind = 4 ) ledg
  integer ( kind = 4 ) lr
  integer ( kind = 4 ) lrline
  integer ( kind = 4 ) ltri
  integer ( kind = 4 ) redg
  integer ( kind = 4 ) rtri
  integer ( kind = 4 ) t
  integer ( kind = 4 ) til(3,*)
  integer ( kind = 4 ) tnbr(3,*)
  real    ( prec ) vcl(2,*)
  real    ( prec ) x
  real    ( prec ) y
!
!  Find rightmost visible boundary edge using links, then possibly
!  leftmost visible boundary edge using triangle neighbor info.
!
  if ( ltri == 0 ) then
    ldone = .false.
    ltri = rtri
    ledg = redg
  else
    ldone = .true.
  end if

  do

    l = -tnbr(redg,rtri)
    t = l / 3
    e = mod ( l, 3 ) + 1
    a = til(e,t)

    if ( e <= 2 ) then
      b = til(e+1,t)
    else
      b = til(1,t)
    end if

    lr = lrline(x,y,vcl(1,a),vcl(2,a),vcl(1,b),vcl(2,b),0.0D+00)

    if ( lr <= 0 ) then
      exit
    end if

    rtri = t
    redg = e

  end do

  if ( ldone ) then
    return
  end if

  t = ltri
  e = ledg

  do

    b = til(e,t)

    if ( 2 <= e ) then
      e = e - 1
    else
      e = 3
    end if

    do while ( 0 < tnbr(e,t) ) 

      t = tnbr(e,t)

      if ( til(1,t) == b ) then
        e = 3
      else if ( til(2,t) == b ) then
        e = 1
      else
        e = 2
      end if

    end do

    a = til(e,t)
    lr = lrline(x,y,vcl(1,a),vcl(2,a),vcl(1,b),vcl(2,b),0.0D+00)

    if ( lr <= 0 ) then
      exit
    end if

  end do

  ltri = t
  ledg = e

  return
end subroutine vbedg

subroutine sub_basis

use comvar

implicit none

real(prec):: x1,x2,x3,h1,h2,h3, det
real(prec),dimension(2,2):: A
real(prec),dimension(2,1):: X_0

x1 = grid(1,til(1,itri))
x2 = grid(1,til(2,itri))
x3 = grid(1,til(3,itri))
h1 = grid(2,til(1,itri))
h2 = grid(2,til(2,itri))
h3 = grid(2,til(3,itri))

! 2x2 Matrix
A(1,1) = x2 - x1
A(1,2) = x3 - x2
A(2,1) = h2 - h1
A(2,2) = h3 - h2

! Invert Matrix
det = A(1,1)*A(2,2)-A(1,2)*A(2,1)

A_inv(1,1,itri) = A(2,2) / det
A_inv(1,2,itri) = - A(1,2) / det
A_inv(2,1,itri) = - A(2,1) / det
A_inv(2,2,itri) = A(1,1) / det

X_0(1,1) = x1
X_0(2,1) = h1

! Compute constant
CONST(:,:,itri) = matmul(A_inv(:,:,itri),X_0)

end subroutine sub_basis

end module delaunay_func