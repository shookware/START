!------------------------------------------------------------------------------
! TJU/Department of Mechanics, Fluid Mechanics, Code START
!------------------------------------------------------------------------------
!
!  File: module_sparsematrix.f90
!> @file
!> @breif 稀疏矩阵模块文件.
!  DESCRIPTION:
!>
!!
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
! TJU/Department of Mechanics, Fluid Mechanics, Code START
!------------------------------------------------------------------------------
!
!  module: sparsematrix
!> @breif 稀疏矩阵模块.
!  DESCRIPTION:
!>
!!
!  REVISION HISTORY:
!  2017-08-02 - Initial Version
!  TODO_yyyy-mm-dd - TODO_describe_appropriate_changes - TODO_name
!> @author
!> @date 2017-08-02
!------------------------------------------------------------------------------
module mod_sparse_matrix

    use penf, only: R_P
    implicit none

    private

    !> COO系数矩阵的一个非零元素类型
    type, public :: coo_element_type

        private

        integer, private :: ir   !< 行序号
        integer, private :: jc   !< 列序号
        complex(R_P), private :: Avalue !< 非零值

        Contains

        procedure :: Set => set_coo_element   !设置该元素
        procedure :: Get => get_coo_element   !获得该元素

    end type coo_element_type

    !> COO形式系数矩阵类型
    type, public :: mat_coo_type

        private
        integer, private :: nrow  !< 矩阵维数
        integer, private :: nnz   !< 非零元素总个数
        type(coo_element_type), allocatable, private :: item(:) !< 非零元素
        integer, private :: Ncount=0                             !< 当前非零元素数量

        Contains

          Procedure :: Create => create_coo_mat       !< 分配该类型的存储空间
          procedure :: GetNRow => get_coo_nrow        !< 获得矩阵维数
          procedure :: GetNnz => get_coo_nnz          !< 获得矩阵总的非零元素个数
          procedure :: zeros => zeros_coo             !< 矩阵全部赋0值
          procedure :: Print => print_coo            !< 输出矩阵
          procedure :: ISCreate => ISCreate_coo       !< 判断矩阵是否分配内存
          generic :: Set => set_coo_mat_item, set_coo_mat_value, &
                          set_coo_mat_value_block, set_coo_mat_item_group  !<设置矩阵的值
          generic :: Get => get_coo_mat_item, get_coo_mat_value           !<获得矩阵的值


          procedure, private :: set_coo_mat_item
          procedure, private :: set_coo_mat_value
          procedure, private :: get_coo_mat_item
          procedure :: get_coo_mat_value
          procedure, private :: set_coo_mat_value_block
          procedure, private :: set_coo_mat_item_group

    end type mat_coo_type

    !> CSR矩阵类型
    type, public :: mat_csr_type

        private
        integer, private :: nrow  !< 矩阵维数
        integer, private :: nnz   !< 矩阵非零元素个数
        complex(R_P), private, allocatable :: AN(:) !< 矩阵非零元素
        integer, private, allocatable :: IA(:)      !< 矩阵行信息
        integer, private, allocatable :: JA(:)      !< 矩阵列信息

        Contains

          Procedure :: Create => create_csr_mat       !< 分配CSR矩阵内存
          procedure :: set => set_csr_mat_from_coo    !< 从COO类型设置CSR类型
          procedure :: get => get_csr_mat_value       !< 获得CSR类型数组
          procedure :: zeros => zeros_csr             !< 给所有的值赋零值
          procedure :: clear => clear_csr             !< 清除矩阵
          procedure :: IsCreate => IsCreate_csr       !< 判断矩阵是否分配内存

          generic :: operator(.mx.) => csr_matmul_type_array        !< 矩阵乘向量
          generic :: operator(.tmx.) => csr_matmul_type_array_trans !< 矩阵转置乘向量
          generic :: operator(*) => csr_scale_mult_cmplx            !< 标量乘矩阵
          generic :: operator(+) => csr_add_cmplx                   !< 矩阵加法
          generic :: operator(-) => csr_substract_cmplx             !< 矩阵减法

          procedure, private :: csr_matmul_type_array               !< 矩阵乘向量
          procedure, private :: csr_matmul_type_array_trans         !< 矩阵转置乘向量
          procedure, pass(this), private :: csr_scale_mult_cmplx    !< 标量乘矩阵
          procedure, pass(this), private :: csr_add_cmplx           !< 矩阵加法
          procedure, pass(this), private :: csr_substract_cmplx     !< 矩阵减法

    end type mat_csr_type

    type, public :: mat_bsr_type

        private
        integer, private :: nrow !<块行数
        integer, private :: nnz !<块个数
        integer, private :: blockSize !<块大小
        complex(R_P), private, allocatable :: AN(:) !值
        integer, private, allocatable :: IA(:) !<行信息
        integer, private, allocatable :: JA(:) !<列信息
        integer, private :: iterator !<当前坐标迭代器

        Contains

          Procedure :: Create=>create_bsr_mat
          procedure :: Set=> Set_bsr_item_singler
          procedure :: Get=> get_bsr_mat_value
          procedure :: IsCreate=> ISCreate_bsr
          procedure :: Zeros=>zeros_bsr

    end type mat_bsr_type


    contains

    !> 为COO矩阵分配内存,初始化
    !! @param[in] nrow 矩阵维数
    !! @param[in] nnz 矩阵非零元素个数
    subroutine create_coo_mat(this, nrow, nnz)

        implicit none
        integer , intent(in), value :: nrow, nnz
        class(mat_coo_type), intent(out) :: this

        this%nrow=nrow
        this%nnz=nnz
        if(.not. allocated(this%item)) allocate(this%item(nnz))

    end subroutine create_coo_mat

    !> 判断COO矩阵内存是否已经分配
    !! @return 矩阵是否分配
    !! @retval TRUE 已分配内存
    !! @retval FALSE 未分配内存
    function ISCreate_coo(this)

        implicit none
        class(mat_coo_type), intent(inout) :: this
        logical :: iSCreate_coo

        iSCreate_coo=allocated(this%item)

    end function ISCreate_coo


    !> 判断CSR矩阵内存是否已经分配
    !! @return 矩阵是否分配
    !! @retval TRUE 已分配内存
    !! @retval FALSE 未分配内存
    function ISCreate_csr(this)

        implicit none
        class(mat_csr_type), intent(inout) :: this
        logical :: iSCreate_csr

        iSCreate_csr=allocated(this%AN)

    end function ISCreate_csr

    !> 判断CSR矩阵内存是否已经分配
    !! @return 矩阵是否分配
    !! @retval TRUE 已分配内存
    !! @retval FALSE 未分配内存
    function ISCreate_bsr(this)

        implicit none
        class(mat_bsr_type), intent(inout) :: this
        logical :: iSCreate_bsr

        iSCreate_bsr=allocated(this%AN)

    end function ISCreate_bsr

    !> 为CSR矩阵分配内存,初始化
    !! @param[in] nrow 矩阵维数
    !! @param[in] nnz 矩阵非零元素个数
    subroutine create_csr_mat(this, nrow, nnz)

        implicit none
        integer , intent(in) :: nrow, nnz
        class(mat_csr_type), intent(out) :: this

        this%nrow=nrow; this%nnz=nnz
        if(.not. allocated(this%an))then
          allocate(this%AN(nnz))
          allocate(this%IA(nrow+1))
          allocate(this%JA(nnz))
        endif

    end subroutine create_csr_mat

    subroutine Create_bsr_mat(this, nrow, nnz, blocksize)

        implicit none
        class(mat_bsr_type),intent(inout) :: this
        integer, intent(in) :: nrow
        integer, intent(in) :: nnz
        integer, intent(in) :: blocksize

        this%nrow=nrow
        this%nnz=nnz
        this%blockSize=blocksize
        if( .not. allocated(this%an)) then
            allocate(this%an(nnz*blocksize**2))
            allocate(this%ia(this%nrow+1))
            allocate(this%ja(this%nnz))
        end if
        this%an=0.0d0
        this%ia=1
        this%ja=0
        this%iterator=1

    end subroutine Create_bsr_mat

    subroutine Set_bsr_item_singler(this, mat, i_index, j_index)

        implicit none
        class(mat_bsr_type),intent(inout) :: this
        complex(R_P), intent(in) :: mat(this%blockSize, this%blockSize)
        integer, intent(in) :: i_index
        integer, intent(in) :: j_index
        integer :: ibegin, iend

        this%JA(this%iterator)=j_index
        this%ia(i_index+1:this%nrow+1)=this%ia(i_index+1:this%nrow+1)+1
        ibegin= (this%iterator-1)*this%blocksize**2+1
        iend= (this%iterator)*this%blocksize**2
        this%an(ibegin:iend)=reshape(mat, [size(mat)])
        this%iterator=this%iterator+1

    end subroutine Set_bsr_item_singler


    !> 用COO矩阵给CSR矩阵赋值
    !! @param[in] coomat COO形式矩阵
    subroutine set_csr_mat_from_coo(this, coomat)

        use mod_sparsekit
        implicit none
        class(mat_csr_type), intent(inout) :: this
        type(mat_coo_type), intent(inout) :: coomat
        integer :: nrow, nnz
        integer :: ir(coomat%nnz), jc(coomat%nnz)
        complex(R_P) :: Avalue(coomat%nnz)

        call coomat%get(ir, jc, Avalue)
        if(.not. allocated(this%an)) then
            call this%Create(coomat%nrow, coomat%nnz)
        endif

        call coocsr (this%nrow, this%nnz, avalue, ir, jc, this%an, this%ja, this%ia)

    end subroutine set_csr_mat_from_coo

    !> 给coo矩阵非零元素赋值
    !! @param[in] ir 行坐标
    !! @param[in] jc 列坐标
    !! @param[in] Avalue 元素的值
    elemental subroutine set_coo_element(this, ir, jc, Avalue)

        implicit none
        class(coo_element_type), intent(inout) :: this
        integer, intent(in) :: ir, jc
        complex(R_P), intent(in) :: Avalue

        this%ir=ir
        this%jc=jc
        this%Avalue=Avalue

    end subroutine set_coo_element

    !> 获得coo矩阵的一个非零元素
    !! @param[out] ir 行坐标
    !! @param[out] jc 列坐标
    !! @param[out] Avalue 元素的值
    subroutine get_coo_element(this, ir, jc, Avalue)

        implicit none
        class(coo_element_type), intent(in) :: this
        integer , intent(out) :: ir, jc
        complex(R_P), intent(out) :: Avalue

        ir=this%ir
        jc=this%jc
        Avalue=this%Avalue

    end subroutine get_coo_element

    !> 添加n个元素到coo矩阵
    !! @param[in] n 添加元素个数
    !! @param[in] item coo元素
    subroutine set_coo_mat_item_group(this, n, item)

        implicit none
        integer, intent(in) :: n
        type(coo_element_type), intent(in) :: item(n)
        class(mat_coo_type), intent(inout) :: this

        this%item(this%ncount+1:this%ncount+n)=item
        this%ncount=this%ncount+n

    end subroutine set_coo_mat_item_group

    !> 添加1个元素到coo矩阵
    !! @param[in] item coo元素
    elemental subroutine set_coo_mat_item(this, item)

        implicit none
        type(coo_element_type), intent(in) :: item
        class(mat_coo_type), intent(inout) :: this

        this%ncount=this%ncount+1
        this%item(this%ncount)=item

    end subroutine set_coo_mat_item

    !> 获得所有的coo矩阵的非零元素
    !! @param[out] item coo矩阵非零元素
    subroutine get_coo_mat_item(this, item)

        implicit none
        class(mat_coo_type), intent(in) :: this
        type(coo_element_type), intent(out) :: item(this%nnz)

        item=this%item

    end subroutine get_coo_mat_item

    !> 给coo矩阵添加一个非零元素
    !! @param[in] ir 元素的行坐标
    !! @param[in] jc 元素的列坐标
    !! @param[in] Avalue  非零元素的值
    elemental subroutine set_coo_mat_value(this, ir, jc, Avalue)

      implicit none
      class(mat_coo_type), intent(inout) :: this
      integer, intent(in) :: ir, jc
      complex(R_P), intent(in) :: Avalue
      type(coo_element_type) :: item

      call item%set(ir, jc, Avalue)
      call this%set(item)

    end subroutine set_coo_mat_value

    !> 给coo矩阵添加一块非零元素(块方式)
    !! @param[in] 块的自由度(dof*dof的块)
    !! @param[in] iindex 块在整个矩阵中的行序号
    !! @param[in] jindex 块在整个矩阵中的列序号
    !! @param[in] Avalue 非零元素块
    subroutine set_coo_mat_value_block(this, dof, iindex, jindex, Avalue)

        implicit none
        class(mat_coo_type), intent(inout) :: this
        integer , intent(in) :: dof, iindex, jindex
        complex(R_P), intent(in) :: Avalue(dof, dof)
        integer :: i, j
        integer :: n
        type(coo_element_type) :: item(dof*dof)

        n=dof*dof

        do j=1, dof
          do i=1, dof
            call item(i+(j-1)*dof)%Set((iindex-1)*dof+i, (jindex-1)*dof+j, Avalue(i, j))
          end do
        end do
        call this%set(n, item)

    end subroutine set_coo_mat_value_block

    !> 获得COO矩阵全部的非零值
    !! @param[out] ir 非零元素的行坐标
    !! @param[out] jc 非零元素的列坐标
    !! @param[out] Avalue 非零元素
    subroutine get_coo_mat_value(this, ir, jc, Avalue)

      implicit none
      class(mat_coo_type), intent(inout) :: this
      integer, intent(out) :: ir(this%nnz), jc(this%nnz)
      complex(R_P), intent(out) :: Avalue(this%nnz)
      type(coo_element_type) :: item(this%nnz)
      integer :: I

          call this%get(item)
          do I=1, this%nnz
            call item(i)%get(ir(i), jc(i), Avalue(i))
          enddo

    end subroutine get_coo_mat_value

    !> 获得CSR矩阵全部的非零值
    !! @param[out] ia 非零元素的行信息
    !! @param[out] ja 非零元素的列信息
    !! @param[out] An 非零元素
    subroutine get_csr_mat_value(this, ia, ja, an)

        implicit none
        class(mat_csr_type), target, intent(in) :: this
        integer, pointer, intent(out) :: ia(:), ja(:)
        complex(R_P), pointer, intent(out) :: AN(:)

        IA => this%IA; JA => this%JA; AN => this%AN

    end subroutine get_csr_mat_value


    !> 获得BSR矩阵全部的非零值
    !! @param[out] ia 非零元素的行信息
    !! @param[out] ja 非零元素的列信息
    !! @param[out] An 非零元素
    subroutine get_bsr_mat_value(this, ia, ja, an)

        implicit none
        class(mat_bsr_type), target, intent(in) :: this
        integer, pointer, intent(out) :: ia(:), ja(:)
        complex(R_P), pointer, intent(out) :: AN(:)

        IA => this%IA; JA => this%JA; AN => this%AN

    end subroutine get_bsr_mat_value

    !> 获得COO矩阵的维数
    !! @param[out] nrow COO矩阵的维数
    subroutine get_coo_nrow(this, nrow)

        implicit none
        class(mat_coo_type), intent(in) :: this
        integer, intent(out) :: nrow

        nrow=this%nrow

    end subroutine get_coo_nrow

    !> 获得CSR矩阵的维数
    !! @param[out] nrow CSR矩阵的维数
    subroutine get_coo_nnz(this, nnz)

        implicit none
        class(mat_coo_type), intent(in) :: this
        integer, intent(out) :: nnz

        nnz=this%nnz

    end subroutine get_coo_nnz

    !> COO矩阵全部赋零
    subroutine zeros_coo(this)

        implicit none
        class(mat_coo_type), intent(inout) :: this
        integer :: i

        !call this%item%Set(0, 0, (0.0d0, 0.0d0))
        do i=1, this%nnz
          this%item(i)%ir=0
          this%item(i)%jc=0
          this%item(i)%Avalue=0.0d0
        enddo

        this%Ncount=0

    end subroutine zeros_coo

    !> CSR矩阵全部赋零
    subroutine zeros_csr(this)

        implicit none
        class(mat_csr_type), intent(inout) :: this

        this%IA=0; this%JA=0; this%AN=(0.0d0, 0.0d0)

    end subroutine zeros_csr

    !> CSR矩阵清除
    subroutine clear_csr(this)

        implicit none
        class(mat_csr_type), intent(inout) :: this

        this%nnz=0
        this%nrow=0
        deallocate(this%an, this%ia, this%ja)

    end subroutine clear_csr

    !> BSR矩阵全部赋零
    subroutine zeros_bsr(this)

        implicit none
        class(mat_bsr_type), intent(inout) :: this

        this%an=0.0d0
        this%ia=1
        this%ja=0
        this%iterator=1

    end subroutine zeros_bsr

    !> 输出COO矩阵
    subroutine print_coo(this)

      implicit none
      class(mat_coo_type), intent(inout) :: this
      integer :: i

      do i=1, this%nnz
        write(10000, *) this%item(i)%ir, this%item(i)%jc, this%item(i)%Avalue
      end do

    end subroutine print_coo

    !> 矩阵乘向量.
    !!
    !! \f[ \mathbf{B}=\mathbf{A}X \f]
    !! @param[in] array 向量X
    !! @return  矩阵\f$ \mathbf{B} \f$
    function csr_matmul_type_array(this, array) result(outvalue)

       implicit none
       class(mat_csr_type), intent(in) :: this
       complex(R_P), intent(in) :: array(:)
       complex(R_P), allocatable :: outvalue(:)

       if(.not. allocated(outvalue)) allocate(outvalue(size(array)))
       call mkl_zcsrgemv ('n', this%nrow , this%AN, this%IA, this%JA, &
                    &   array, outvalue)

    end function csr_matmul_type_array

    !> 矩阵乘向量.
    !!
    !! \f[ \mathbf{B}=\mathbf{A}^{T}X \f]
    !! @param[in] array 向量X
    !! @return  矩阵\f$ \mathbf{B} \f$
   function csr_matmul_type_array_trans(this, array) result(outvalue)

       implicit none
       class(mat_csr_type), intent(in) :: this
       complex(R_P), intent(in) :: array(:)
       complex(R_P), allocatable :: outvalue(:)

       if(.not. allocated(outvalue)) allocate(outvalue(size(array)))
       call mkl_zcsrgemv ('t', this%nrow , this%AN, this%IA, this%JA, &
                    &   array, outvalue)

   end function csr_matmul_type_array_trans

    !> 标量乘矩阵.
    !!
    !! \f[ \mathbf{B}=\alpha\mathbf{A} \f]
    !! @param[in] scalar 向量\f$\alpha\f$
    !! @return  矩阵\f$ \mathbf{B} \f$
   function csr_scale_mult_cmplx(scalar, this) result(outvalue)

       implicit none
       complex(R_P), intent(in) :: scalar
       class(mat_csr_type), intent(in) :: this
       type(mat_csr_type), allocatable :: outvalue

       if(.not. allocated(outvalue)) allocate(outvalue, source=this)
       outvalue%AN=scalar*this%an

   end function csr_scale_mult_cmplx

    !> CSR矩阵相加.
    !! @return CSR矩阵的和
    function csr_add_cmplx(this, obj1) result(outvalue)

       implicit none
       class(mat_csr_type), intent(in) :: this
       type(mat_csr_type), intent(in) :: obj1
       type(mat_csr_type) :: outvalue
       integer :: ierr, nnz
       complex(R_P), allocatable :: tmp(:)
       integer, allocatable :: tmpI(:)
       integer :: ia(this%nrow+1)

       call  mkl_zcsradd ('n', 1, 3, this%nrow, this%nrow, this%AN, this%JA, &
                &        this%IA, (1.0d0, 0.0d0), obj1%AN, obj1%JA, obj1%IA, &
                &        tmp, tmpI, ia, nnz, ierr)
       call outvalue%Create(this%nrow, ia(this%nrow+1)-1)

       outvalue%ia=ia

       call  mkl_zcsradd ('n', 2, 3, this%nrow, this%nrow, this%AN, this%JA, &
                &   this%IA, (1.0d0, 0.0d0), obj1%AN, obj1%JA, obj1%IA, &
                &   outvalue%AN, outvalue%JA, outvalue%IA, outvalue%nnz, ierr)

    end function csr_add_cmplx

    !> CSR矩阵相减
    !! @return CSR矩阵的差
    function csr_substract_cmplx(this, obj1) result(outvalue)

       implicit none
       class(mat_csr_type), intent(in) :: this
       type(mat_csr_type), intent(in) :: obj1
       type(mat_csr_type) :: outvalue
       integer :: ierr, nnz
       complex(R_P), allocatable :: tmp(:)
       integer , allocatable :: tmpI(:)
       integer :: ia(this%nrow+1)

       call  mkl_zcsradd ('n', 1, 3, this%nrow, this%nrow, this%AN, this%JA, &
                &   this%IA, (1.0d0, 0.0d0), obj1%AN, obj1%JA, obj1%IA, &
                &   tmp, tmpI, ia, nnz, ierr)
       call outvalue%Create(this%nrow, ia(this%nrow+1)-1)

       outvalue%ia=ia

       call  mkl_zcsradd ('n', 2, 3, this%nrow, this%nrow, this%AN, this%JA, &
                &   this%IA, (-1.0d0, 0.0d0), obj1%AN, obj1%JA, obj1%IA, &
                &   outvalue%AN, outvalue%JA, outvalue%IA, outvalue%nnz, ierr)

    end function csr_substract_cmplx

end module mod_sparse_matrix
