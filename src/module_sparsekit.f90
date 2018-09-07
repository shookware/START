!------------------------------------------------------------------------------
! TJU/Department of Mechanics, Fluid Mechanics, Code START
!------------------------------------------------------------------------------
!
!  File: module_sparsekit.f90
!> @file
!> @breif sparsekit文件.
!  DESCRIPTION:
!>
!!
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
! TJU/Department of Mechanics, Fluid Mechanics, Code START
!------------------------------------------------------------------------------
!
!  module: sparsekit
!> @breif sparsekit模块.
!  DESCRIPTION:
!>
!!
!  REVISION HISTORY:
!  2017-08-02 - Initial Version
!  TODO_yyyy-mm-dd - TODO_describe_appropriate_changes - TODO_name
!> @author
!> Liu Jianxin
!> @date 2017-08-02
!------------------------------------------------------------------------------
module mod_sparsekit

    use penf, only: R_P
    implicit none

    contains

    !> coo格式稀疏矩阵转换到csr格式
    !!
    !!   This routine converts a matrix that is stored in COO coordinate format
    !!   a, ir, jc into a CSR row general sparse ao, jao, iao format.
    !! @param[in] nrow   the row dimension of the matrix.
    !! @param[in] nnz    the number of nonzero elements.
    !! @note  COO: matrix in coordinate format. a(k), ir(k), jc(k) store the nnz
    !!         nonzero elements of the matrix with a(k) = actual real value of
    !!         the elements, ir(k) = its row number and jc(k) = its column
    !!        number. The order of the elements is arbitrary.
    !! @param[in] a   actual real value of the elements.
    !! @param[in] ir  its row number.
    !! @param[in] jc  its column number.
    !! @note  ao(*), jao(*), iao(nrow+1), the matrix in CSR
    !!              Compressed Sparse Row format.
    !! @param[out] ao   CSR矩阵非零元素
    !! @param[out] jao  CSR矩阵列信息
    !! @param[out] iao  CSR矩阵行信息
    !! @author Youcef Saad
    !! @date 2004-01-07
     subroutine coocsr ( nrow, nnz, a, ir, jc, ao, jao, iao )

    !*****************************************************************************80
    !
    !! COOCSR converts COO to CSR.
    !
    !  Discussion:
    !
    !
    !    Output, real AO(*), JAO(*), IAO(NROW+1), the matrix in CSR
    !    Compressed Sparse Row format.
    !
      implicit none

      integer ( kind = 4 ) nrow

      complex(R_P) a(:)
      complex(R_P) ao(:)
      integer ( kind = 4 ) i
      integer ( kind = 4 ) iad
      integer ( kind = 4 ) iao(nrow+1)
      integer ( kind = 4 ) ir(:)
      integer ( kind = 4 ) j
      integer ( kind = 4 ) jao(:)
      integer ( kind = 4 ) jc(:)
      integer ( kind = 4 ) k
      integer ( kind = 4 ) k0
      integer ( kind = 4 ) nnz
      complex(R_P) x


      iao(1:nrow+1) = 0
    !
    !  Determine the row lengths.
    !
      do k = 1, nnz
        iao(ir(k)) = iao(ir(k)) + 1
      end do
    !
    !  The starting position of each row.
    !
      k = 1
      do j = 1, nrow+1
         k0 = iao(j)
         iao(j) = k
         k = k + k0
      end do
    !
    !  Go through the structure once more.  Fill in output matrix.
    !
      do k = 1, nnz
         i = ir(k)
         j = jc(k)
         x = a(k)
         iad = iao(i)
         ao(iad) = x
         jao(iad) = j
         iao(i) = iad + 1
      end do
    !
    !  Shift back IAO.
    !
      do j = nrow, 1, -1
        iao(j+1) = iao(j)
      end do
      iao(1) = 1

      return
    end subroutine

    subroutine zilu( n, nz_num, ia, ja, l )

    !*****************************************************************************80
    !
    !! ILU_CR computes the incomplete LU factorization of a matrix.
    !
    !  Discussion:
    !
    !    The matrix A is assumed to be stored in compressed row format.  Only
    !    the nonzero entries of A are stored.  The vector JA stores the
    !    column index of the nonzero value.  The nonzero values are sorted
    !    by row, and the compressed row vector IA then has the property that
    !    the entries in A and JA that correspond to row I occur in indices
    !    IA(I) through IA(I+1)-1.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    27 July 2007
    !
    !  Author:
    !
    !    Original C version by Lili Ju.
    !    FORTRAN90 version by John Burkardt.
    !
    !  Parameters:
    !
    !    Input, integer ( kind = 4 ) N, the order of the system.
    !
    !    Input, integer ( kind = 4 ) NZ_NUM, the number of nonzeros.
    !
    !    Input, integer ( kind = 4 ) IA(N+1), JA(NZ_NUM), the row and column
    !    indices of the matrix values.  The row vector has been compressed.
    !
    !    Input, real ( kind = 8 ) A(NZ_NUM), the matrix values.
    !
    !    Input, integer ( kind = 4 ) UA(N), the index of the diagonal element
    !    of each row.
    !
    !    Output, real ( kind = 8 ) L(NZ_NUM), the ILU factorization of A.
    !
      implicit none

      integer ( kind = 4 ) n
      integer ( kind = 4 ) nz_num

      !complex ( kind = 8 ) a(nz_num)
      integer ( kind = 4 ) i
      integer ( kind = 4 ) ia(n+1)
      integer ( kind = 4 ) iw(n)
      integer ( kind = 4 ) j
      integer ( kind = 4 ) ja(nz_num)
      integer ( kind = 4 ) jj
      integer ( kind = 4 ) jrow
      integer ( kind = 4 ) jw
      integer ( kind = 4 ) k
      complex ( kind = 8 ) l(nz_num)
      complex ( kind = 8 ) tl
      integer ( kind = 4 ) ua(n)
    !
    !  Copy A.
    !
     ! l(1:nz_num) = a(1:nz_num)

    call diagonal_pointer_cr ( n, nz_num, ia, ja, ua )

      do i = 1, n
    !
    !  IW points to the nonzero entries in row I.
    !
        iw(1:n) = -1

        do k = ia(i), ia(i+1) - 1
          iw(ja(k)) = k
        end do

        do j = ia(i), ia(i+1) - 1
          jrow = ja(j)
          if ( i <= jrow ) then
            exit
          end if
          tl = l(j) * l(ua(jrow))
          l(j) = tl
          do jj = ua(jrow) + 1, ia(jrow+1) - 1
            jw = iw(ja(jj))
            if ( jw /= -1 ) then
              l(jw) = l(jw) - tl * l(jj)
            end if
          end do
        end do

        ua(i) = j

        if ( jrow /= i ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'ILU_CR - Fatal error!'
          write ( *, '(a)' ) '  JROW ~= I'
          write ( *, '(a,i8)' ) '  JROW = ', jrow
          write ( *, '(a,i8)' ) '  I    = ', i
          stop
        end if

        if ( l(j) == 0.0D+00 ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'ILU_CR - Fatal error!'
          write ( *, '(a,i8)' ) '  Zero pivot on step I = ', i
          write ( *, '(a,i8,a)' ) '  L(', j, ') = 0.0'
          stop
        end if

        l(j) = 1.0D+00 / l(j)

      end do

      l(ua(1:n)) = 1.0D+00 / l(ua(1:n))

      return
    end subroutine zilu

     subroutine diagonal_pointer_cr ( n, nz_num, ia, ja, ua )

    !*****************************************************************************80
    !
    !! DIAGONAL_POINTER_CR finds diagonal entries in a sparse compressed row matrix.
    !
    !  Discussion:
    !
    !    The matrix A is assumed to be stored in compressed row format.  Only
    !    the nonzero entries of A are stored.  The vector JA stores the
    !    column index of the nonzero value.  The nonzero values are sorted
    !    by row, and the compressed row vector IA then has the property that
    !    the entries in A and JA that correspond to row I occur in indices
    !    IA[I] through IA[I+1]-1.
    !
    !    The array UA can be used to locate the diagonal elements of the matrix.
    !
    !    It is assumed that every row of the matrix includes a diagonal element,
    !    and that the elements of each row have been ascending sorted.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    18 July 2007
    !
    !  Author:
    !
    !    Original C version by Lili Ju.
    !    FORTRAN90 version by John Burkardt.
    !
    !  Parameters:
    !
    !    Input, integer ( kind = 4 ) N, the order of the system.
    !
    !    Input, integer ( kind = 4 ) NZ_NUM, the number of nonzeros.
    !
    !    Input, integer ( kind = 4 ) IA(N+1), JA(NZ_NUM), the row and column
    !    indices of the matrix values.  The row vector has been compressed.
    !    On output, the order of the entries of JA may have changed because of
    !    the sorting.
    !
    !    Output, integer ( kind = 4 ) UA(N), the index of the diagonal element
    !    of each row.
    !
      implicit none

      integer ( kind = 4 ) n
      integer ( kind = 4 ) nz_num

      integer ( kind = 4 ) i
      integer ( kind = 4 ) ia(n+1)
      integer ( kind = 4 ) k
      integer ( kind = 4 ) ja(nz_num)
      integer ( kind = 4 ) ua(n)

      ua(1:n) = -1

      do i = 1, n
        do k = ia(i), ia(i+1) - 1
          if ( ja(k) == i ) then
            ua(i) = k
          end if
        end do
      end do

      return
     end subroutine
     
subroutine bsrcsr ( n, nblk, na, a, ja, ia, ao, jao, iao )

!*****************************************************************************80
!
!! BSRCSR converts Block Sparse Row to Compressed Sparse Row (CSR) format.
!
!  Discussion:
!
!    This routine converts a matrix stored in block-reduced
!    a, ja, ia format to the general sparse row a, ja, ia format.
!    A matrix that has a block structure is a matrix whose entries
!    are blocks of the same size nblk (e.g. 3 x 3).  Then it is often
!    preferred to work with the reduced graph of the matrix, i.e.,
!    Instead of storing one element at a time one can store the whole
!    block.  In this storage scheme a row of the array a will
!    hold the nblk**2 entries of a block.
!
!    This code is not "in place".
!
! general picture: (nblk = 2)
!     --- A ---                                --- JA --  -- IA --
! A=  x x x x   1st block in block row 1           x         x
!     x x x x  2-nd block in block row 1           x
!     . . . .                                      .
!     x x x x  last block in block row 1           x
!     -------                                     ---
!     x x x x   1st block in block row 2           x          x
!     x x x x  2-nd block in block row 2           x
!     . . . .                                      x
!     x x x x   last block in block row 2          x
!     -------                                     ---
!     .......                                     ...         .
!     -------                                     ---
!     x x x x   1st block in block row n/nblk      x          x
!     x x x x  2-nd block in block row n/nblk      x
!     . . . .                                      x
!     x x x x  last block in block row n/nblk      x
!     -------                                     ---
!                                               end + 1       x
!
!
! example  with nblk = 2:
!
!
!             1   2   0   0   3   4
!             5   6   0   0   7   8
!             0   0   9  10  11  12
!             0   0  13  14  15  16
!             17 18   0   0   0   0
!             22 23   0   0   0   0
! THEN:
!
!  ---- A ----                                     -- JA --   -- IA --
!-
!  1   5   2  6  Block row 1 (2 block matrices)      | 1  <--- | 1
!  3   7   4  8                                      | 5       |
!  ------------                                      | --      |
!  9  13  10 14  block row 2 (2 block matrices)      | 3  <--- | 3
! 11  15  12 16                                      | 5       |
!  ------------                                      | --      |
! 17  22  18 23  Block row 3 (1 block matrix)        | 1  <--- | 5
!  ------------                                      | --      |
!                                                   end+1 <--- | 6
!
! JA  =  1  5 | 3  5 | 1       column numbers of (1,1) entries of blocks
! IA  =  1      3      5  6    pointers to beginnings of BLOCK-rows
!
!  Modified:
!
!    07 January 2004
!
!  Author:
!
!    Youcef Saad
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NROW, the row dimension of the matrix.
!
! nblk  = integer ( kind = 4 ) equal to the dimension of each block.
!         nblk must divide n.
!
! na      = first dimension of array a as declared in calling program
!
! a      = real array containing the values of the matrix. For details
!         on the format see below. Each row of a contains the nblk x nblk
!         block matrix unpacked column-wise (this allows the user to
!         declare the array a as a(na,nblk,nblk) on entry if desired).
!         the block rows are stored in sequence just as for the compressed
!         sparse row format.
!
! ja      = integer ( kind = 4 ) array of length n/nblk. ja(k) contains the 
!         column index of the leading element, i.e., the element (1,1) of the
!         block that is held in the row a(k,*) of the value array.
!
! ia    = integer ( kind = 4 ) array of length n/nblk+1. ia(i) points to the 
!        beginning of block row number i in the arrays a and ja.
!
!    Output, real AO(*), JAO(*), IAO(N+1), the matrix in CSR
!    Compressed Sparse Row format.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) na

  real ( kind = 8 ) a(na,*)
  real ( kind = 8 ) ao(*)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i1
  integer ( kind = 4 ) i2
  integer ( kind = 4 ) ia(*)
  integer ( kind = 4 ) iao(n+1)
  integer ( kind = 4 ) ii
  integer ( kind = 4 ) ij
  integer ( kind = 4 ) irow
  integer ( kind = 4 ) j
  integer ( kind = 4 ) ja(*)
  integer ( kind = 4 ) jao(*)
  integer ( kind = 4 ) jj
  integer ( kind = 4 ) jst
  integer ( kind = 4 ) k
  integer ( kind = 4 ) krow
  integer ( kind = 4 ) nblk
  integer ( kind = 4 ) nr
!
!  Get the IA, JA data structure for output matrix
!
  nr = n / nblk

  iao(1:n+1) = 0

  irow = 0
  krow = 1

  do ii = 1, nr
!
!  NR is the dimension of the reduced matrix.
!
     i1 = ia(ii)
     i2 = ia(ii+1) - 1
!
!  Create NBLK rows for each K.
!
     do i = 1, nblk
        do k = i1, i2
           jst = ja(k) - 1
           do j = 1, nblk
              ij = ( j - 1 ) * nblk + i
              ao(krow) = a(k,ij)
              jao(krow) = jst + j
              krow = krow + 1
           end do
         end do
      iao(irow+i) = krow
     end do

     irow = irow + nblk

  end do

  do jj = 1, n
     j = n - jj + 1
     iao(j+1) = iao(j)
  end do

  iao(1) = 1

  return
end subroutine     

end module mod_sparsekit
