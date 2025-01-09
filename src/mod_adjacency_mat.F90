module mod_adjacency_mat
    !! This module is used to efficiently handle the topological problem of
    !! finding all the atom pairs in the molecule separated by exactly n bonds.
    !! To do so the molecule is tought as an undirected unweighted graph, 
    !! represented in memory with a boolean adjacency matrix. Since this matrix
    !! is sparse and symmetric, it is effeciently represented in the Yale format
    !! (sometime referred as compressed sparse row, CSR) omitting the value
    !! array.
    !!
    !! From the adjacency matrix \(\mathbb{C}_1\), the matrix containing the atom
    !! pairs separated by exactly n bonds are computed according the following
    !! recoursive formula:
    !! \[\mathbb{C}_n := \mathbb{C}_{n-1} \cdot \mathbb{C}_1 - (\mathbb{C}_{n-1}
    !! + \mathbb{C}_{n-2});\]
    !! Note that in this context since we are manipulating boolean matrices, 
    !! the product is replaced by logical and, sum by logical or, and - by  
    !! logical and not. To apply this recursive definition, we also assume that
    !! \(\mathbb C_0 := \mathbb 1\).
    !! 
    !! The module contains all the routines needed to perform such operations,
    !! including the linear scaling multiplication between sparse matrices.

    use mod_memory, only: ip

    implicit none
    private

    type yale_sparse
        !! Describe a square boolean sparse matrix n x n.  
        !! Indices of non-zero elements in the i-th row are contained in
        !! the array ci(ri(i):ri(i+1)).
        integer(ip) :: n 
        !! Dimension of the matrix
        integer(ip), allocatable :: ri(:) 
        !! Row index vector, always dimensioned n+1
        integer(ip), allocatable :: ci(:) 
        !! Column index vector, dimension number of non-zero elements
    end type yale_sparse

    public :: yale_sparse
    public :: adj_mat_from_conn, build_conn_upto_n, free_yale_sparse, copy_yale_sparse, &
              reallocate_mat, reverse_grp_tab, &
              compress_list, compress_data, allocate_yale_sparse

    contains
        subroutine print_yale_sparse(m, header, level, logpre, u)
            use mod_io, only: ommp_message
            use mod_constants, only: OMMP_STR_CHAR_MAX

            implicit none 

            type(yale_sparse), intent(in) :: m 
            !! Matrix to print
            character(len=*), intent(in), optional :: header
            !! Header of the matrix printout
            character(len=*), intent(in), optional :: logpre
            !! String that explains the type of message, if missing a default type
            !! is assigned based on requested verbosity level
            integer(ip), intent(in) :: level
            !! Requested verbosity level
            integer(ip), intent(in), optional :: u
            !! Output unit for the message, if missing, [[iof_mmpol]] is used.

            integer(ip), parameter :: max_print_size = 20
            integer(ip) :: i, j
            character(len=OMMP_STR_CHAR_MAX) :: s

            write(s, "(A, A)") "Yale Sparse Matrix :: ", header
            call ommp_message(s, level, logpre, u)

            if(m%n < max_print_size) then
                do i=1, m%n
                    write(s, "(I2, ':', 20I2)") i, m%ci(m%ri(i):m%ri(i+1)-1)
                    call ommp_message(s, level, logpre, u)
                end do
                write(s, "(A, A)") "End Yale Sparse Matrix :: ", header
                call ommp_message(s, level, logpre, u)
            else
                write(s, "(A, A, A)") "****** Matrix too big to be printed (", header, ")"
                call ommp_message(s, level, logpre, u)
            end if
        end subroutine

        subroutine allocate_yale_sparse(m, n, nnz)
            implicit none

            type(yale_sparse), intent(inout) :: m
            !! Matrix to be initialized
            integer(ip), intent(in) :: n
            !! Major dimension for the matrix
            integer(ip), intent(in) :: nnz
            !! Number of nonzero element for the matrix

            m%n = n
            allocate(m%ri(n+1))
            allocate(m%ci(nnz))
        end subroutine

        subroutine copy_yale_sparse(f, t)
            !! Copy boolean sparse matrix in yale format f to t.
            implicit none
            type(yale_sparse), intent(in) :: f
            !! Matrix to copy
            type(yale_sparse), intent(out) :: t
            !! Destination matrix

            allocate(t%ri(size(f%ri)))
            allocate(t%ci(size(f%ci)))
            
            t%n = f%n
            t%ri = f%ri
            t%ci = f%ci
        end subroutine copy_yale_sparse

        subroutine free_yale_sparse(m)
            !! Deallocate a boolean sparse matrix (in Yale format) 
            implicit none
            
            type(yale_sparse), intent(inout) :: m !! Matrix to free

            if(allocated(m%ri)) deallocate(m%ri)
            if(allocated(m%ci)) deallocate(m%ci)
            m%n=0_ip
        end subroutine free_yale_sparse

        subroutine adj_mat_from_conn(i12, sparse)
            !! Create adjacency matrix \(\mathbb C_1\) from connectivity lists.   
            !! Array i12 and n12 contain the connectivity list in the following
            !! format: i12(0:n(j),j) contains the index of all the atoms connected
            !! to atom with index j.

            use mod_utils, only: sort_ivec_inplace
            use mod_memory, only: mallocate, mfree
            implicit none

            integer(ip), intent(inout) :: i12(:,:)
            !! Indices of connected atoms for each atom in the molecule
            type(yale_sparse), intent(out) :: sparse
            !! Adjacency matrix in Yale format (\(\mathbb C_1\))

            integer(ip) :: i, j, n, m
            integer(ip), allocatable :: nnz(:)

            n = size(i12, 2)
            m = size(i12, 1)
            call mallocate('adj_mat_from_conn [nnz]', n, nnz)

            !$omp parallel do default(shared) &
            !$omp private(j,i)
            do i=1, n
                ! Count the number of non-zero elements and move them to the left
                nnz(i) = 0
                do j=1, m
                    if(i12(j,i) /= 0) then
                        nnz(i) = nnz(i) + 1
                        i12(nnz(i),i) = i12(j,i)
                    end if
                end do
                call sort_ivec_inplace(i12(1:nnz(i),i))
            end do

            call compress_list(n, i12, nnz, sparse)
            
            call mfree('adj_mat_from_conn [nnz]', nnz)

        end subroutine adj_mat_from_conn

        subroutine reallocate_mat(m, nnz)
            !! Reshape a boolean sparse matrix in Yale format to accomodate a 
            !! larger number of non-zero elements or to trim unused non-zero
            !! elements after a guess allocation.

            implicit none

            type(yale_sparse), intent(inout) :: m
            !! Matrix to be reshaped
            integer(ip), intent(in) :: nnz
            !! New shape for the matrix
            
            integer(ip), allocatable :: tmp(:)

            if(nnz == 0) then
                ! The matrix is empty, just allocate empty ci
                deallocate(m%ci) 
                allocate(m%ci(0))
            else
                allocate(tmp(nnz))
                if(size(m%ci) > nnz) then
                    ! We are shrinking the matrix
                    tmp = m%ci(1:nnz)
                else
                    ! We are enlarging the matrix
                    tmp(1:size(m%ci)) = m%ci
                end if

                deallocate(m%ci)
                allocate(m%ci(nnz))
                m%ci = tmp
                deallocate(tmp)
            end if
        end subroutine reallocate_mat

        subroutine mat_mult2(sp1, sp2, res)
            !! Performs the operation \(res := sp1 \cdot sp2\) on boolean
            !! sparse matrices; product correspond to logical and while sum
            !! correspond to logical or.    
            !! This subroutine is just for test pourposes as it performs the 
            !! matrix product without exploiting the sparsity, and therefore 
            !! with an unfavorable scaling \(\mathcal O(n^2)\).
            implicit none
            
            type(yale_sparse), intent(in) :: sp1, sp2
            !! Input matrices (sparse boolean matrices in Yale format)
            type(yale_sparse), intent(out) :: res
            !! Output matrix (sparse boolean matrices in Yale format)

            integer(ip) :: ic, ir, m1, n1, k, nnz, res_nnz
            
            res%n = sp1%n
            allocate(res%ri(res%n+1))
            
            ! This is just a guess, it could be increased later
            nnz = max(size(sp1%ci), size(sp2%ci))
            res_nnz = nnz
            allocate(res%ci(res_nnz))
            
            res%ri(1) = 1
            do ir = 1, sp1%n
                m1 = sp1%ri(ir)
                n1 = sp1%ri(ir+1) - 1

                res%ri(ir+1) = res%ri(ir)

                do ic = 1, sp2%n
                    do k = sp2%ri(ic), sp2%ri(ic+1)-1
                        if(any(sp1%ci(m1:n1) == sp2%ci(k))) then
                            res%ci(res%ri(ir+1)) = ic
                            res%ri(ir+1) = res%ri(ir+1) + 1

                            if(res%ri(ir+1) > res_nnz) then
                                res_nnz = res_nnz + nnz
                                call reallocate_mat(res, res_nnz)
                            endif
                            
                            exit
                        end if
                    end do
                end do
            end do
            
            ! Trim the output vector
            nnz = res%ri(res%n+1) - 1
            call reallocate_mat(res, nnz)    
        end subroutine
        
        subroutine mat_mult(sp1, sp2, res)
            !! Performs the operation \(res := sp1 \cdot sp2\) on boolean
            !! sparse matrices; product correspond to logical and while sum
            !! correspond to logical or.    
            !! This subroutine is the one actually used in the code; It performs 
            !! matrix product exploiting the sparsity, and therefore 
            !! with a scaling \(\mathcal O(n)\).
            implicit none
            
            type(yale_sparse), intent(in) :: sp1, sp2
            !! Input matrices (sparse boolean matrices in Yale format)
            type(yale_sparse), intent(out) :: res
            !! Output matrix (sparse boolean matrices in Yale format)

            integer(ip) :: ic, ir, k, nnz, res_nnz
            
            res%n = sp1%n
            allocate(res%ri(res%n+1))
            
            ! This is just a guess, it could be increased later
            nnz = max(size(sp1%ci), size(sp2%ci))
            res_nnz = nnz
            allocate(res%ci(res_nnz))
            
            res%ri(1) = 1
            do ir = 1, sp1%n
                res%ri(ir+1) = res%ri(ir)
                do k = sp1%ri(ir), sp1%ri(ir+1) - 1
                    do ic = sp2%ri(sp1%ci(k)), sp2%ri(sp1%ci(k)+1)-1
                        ! ir -> sp2%ci(ic)
                         if(.not. any(res%ci(res%ri(ir):res%ri(ir+1)-1) == sp2%ci(ic))) then
                             res%ci(res%ri(ir+1)) = sp2%ci(ic)
                             res%ri(ir+1) = res%ri(ir+1) + 1

                             if(res%ri(ir+1) > res_nnz) then
                                 res_nnz = res_nnz + nnz
                                 call reallocate_mat(res, res_nnz)
                             endif
                         end if
                    end do
                end do
            end do
            
            ! Trim the output vector
            nnz = res%ri(res%n+1) - 1
            call reallocate_mat(res, nnz)    

        end subroutine

        subroutine mat_andnot(sp1, sp2, res)
            !! Performs the operation \(res := sp1 \land \neg sp2\) on
            !! boolean sparse matrices.

            implicit none
            
            type(yale_sparse), intent(in) :: sp1, sp2
            !! Input matrices (sparse boolean matrices in Yale format)
            type(yale_sparse), intent(out) :: res
            !! Output matrix (sparse boolean matrices in Yale format)

            integer(ip) :: ic, ir, ic2_1, ic2_2, nnz

            res%n = sp1%n
            allocate(res%ri(res%n+1))
            !! Worst case scenario, nnz(res) will be nnz(sp1)
            allocate(res%ci(sp1%ri(sp1%n+1))) 
            
            res%ri(1) = 1
            do ir = 1, sp1%n
                res%ri(ir+1) = res%ri(ir)
                do ic = sp1%ri(ir), sp1%ri(ir+1) - 1
                    ic2_1 = sp2%ri(ir)
                    ic2_2 = sp2%ri(ir+1) - 1
                    if(.not. any(sp2%ci(ic2_1:ic2_2) == sp1%ci(ic))) then
                        res%ci(res%ri(ir+1)) = sp1%ci(ic)
                        res%ri(ir+1) = res%ri(ir+1) + 1
                    end if
                end do
            end do
            
            ! Trim the output vector
            nnz = res%ri(res%n+1) - 1
            call reallocate_mat(res, nnz)    

        end subroutine

        subroutine sparse_identity(n, res)
            !! Create an identity matrix (boolean sparse, represented in
            !! Yale format) of dimension \(n\).
            implicit none

            integer(ip), intent(in) :: n
            !! Rank of the output matrix
            type(yale_sparse), intent(out) :: res
            !! Output matrix

            integer(ip) :: i

            res%n = n
            allocate(res%ci(res%n))
            allocate(res%ri(res%n+1))
            
            !$omp parallel do
            do i = 1, n
                res%ci(i) = i
                res%ri(i) = i
            end do
            res%ri(n+1) = n+1

        end subroutine sparse_identity

        subroutine build_conn_upto_n(adj, n, res, start_id)
            !! Build connectivity matrices up to \(\mathbb C_n\) 
            !! starting from \(\mathbb C_1\). Results are stored in an
            !! array of boolean sparse matrix in Yale format in such a way that
            !! \(res(n) := \mathbb C_n\); since FORTRAN is 1-based the useless
            !! \(\mathbb C_0\) is not stored.
            implicit none
            
            type(yale_sparse), intent(in) :: adj
            !! Adjacency matrix in Yale format
            integer(ip), intent(in) :: n
            !! Maximum level of connectivity that should be computed
            type(yale_sparse), intent(out), allocatable :: res(:)
            !! Results connectivity matrices
            logical :: start_id
            !! Specifies if the first matrix allocated res(1) should be the 
            !! identity (true) or the adjacency (false).

            integer(ip) :: i, adj_idx
            type(yale_sparse) :: tmp, id

            if(start_id) then
                allocate(res(n+1))
                call sparse_identity(adj%n, res(1))
                adj_idx = 2
            else
                allocate(res(n))
                adj_idx = 1
            end if
            call copy_yale_sparse(adj, res(adj_idx))
            
            do i=adj_idx+1, adj_idx+n-1
                if(size(res(i-1)%ci) == 0) then
                    ! Create a null matrix
                    res(i)%n = res(i-1)%n
                    allocate(res(i)%ri(res(i)%n+1))
                    res(i)%ri = 1
                    allocate(res(i)%ci(0))
                else
                    call mat_mult(res(i-1), res(adj_idx), res(i))
                    !call matcpy(res(i-1), res(i))
                    call mat_andnot(res(i), res(i-1), tmp)
                    if(i == adj_idx+1) then
                        call sparse_identity(adj%n, id)
                        call mat_andnot(tmp, id, res(i))
                        if(start_id) then
                            call copy_yale_sparse(id, res(1))
                        end if
                        call free_yale_sparse(id)
                    else
                        call mat_andnot(tmp, res(i-2), res(i))
                    end if
                end if
            end do
            call free_yale_sparse(tmp)
        end subroutine build_conn_upto_n
        
        subroutine reverse_grp_tab(a2g, g2a, ng_in)
            use mod_memory, only: mallocate, mfree
            !! Takes as argument an array of  group index for each
            !! atom, and create a list of atms in each group using the
            !! sparse matrix format (saved as Yale format).
            !! This is used by cell lists, polarization group etc.
            
            implicit none

            integer(ip), intent(in) :: a2g(:)
            !! Index of polarization group for each MM atom
            type(yale_sparse), intent(out) :: g2a
            !! Indices of atoms included in each polarization group;
            !! Atom indeces for the n-th group are found at 
            !! pg2mm%ci(pg2mm%ri(n):pg2mm%ri(n+1)-1)
            integer(ip), intent(in), optional :: ng_in
            !! Number of groups if it is not provided in input it is
            !! assumed that the number of group equals the largest group
            !! index, that is no empty groups are present after the one
            !! with the largest index.

            integer(ip) :: i, j, na, ng, ig
            integer(ip), allocatable :: uc_data(:, :), g_dim(:)

            na = size(a2g)
            if(present(ng_in)) then
                ng = ng_in
            else
                ng = maxval(a2g)
            end if

            ! Find largest group
            call mallocate('reverse_grp_tab [g_dim]', ng, g_dim)
            g_dim = 0

            
            do i=1, na
                g_dim(a2g(i)) = g_dim(a2g(i)) + 1
            end do

            ! Struct for uncompressed data
            call mallocate('reverse_grp_tab [uc_data]', maxval(g_dim), ng, uc_data)
            ! First invert in an uncompressed structure
            uc_data = 0
            g_dim = 0

            do i=1, na
                ig = a2g(i)
                g_dim(ig) = g_dim(ig) + 1
                uc_data(g_dim(ig),ig) = i
            end do
            
            ! Compress the list
            !g2a%n = ng
            !if(.not. allocated(g2a%ri)) &
            !    call mallocate('reverse_grp_tab [ri]', ng+1, g2a%ri)
            !if(.not. allocated(g2a%ci)) &
            !    call mallocate('reverse_grp_tab [ci]', na, g2a%ci)
            !g2a%ri(1) = 1
            !do i=1, ng
            !    g2a%ri(i+1) = g2a%ri(i) + g_dim(i) - 1
            !    g2a%ci(g2a%ri(i):g2a%ri(i+1)-1) = uc_data(1:g_dim(i)-1,i)
            !end do
            call compress_list(ng, uc_data, g_dim, g2a)

            ! Free temporary mem
            call mfree('reverse_grp_tab [uc_data]', uc_data)
            call mfree('reverse_grp_tab [g_dim]', g_dim)

        end subroutine reverse_grp_tab

        subroutine compress_list(n, uc_list, nit, s)
            !! This subroutine takes as input a sparse matrix (rank [n]) in an 
            !! uncompressed yale format [uc_list], as a rectangular 
            !! matrix ([n] x max_el_per_row) and the actual number of items 
            !! [nit] for each row (remaining
            !! elements are not considered) and compress in a Yale format
            !! sparse matrix [s].
            !! The task is parallelized to handle large matrices.
            
            use mod_memory, only: mallocate, mfree

            implicit none

            integer(ip), intent(in) :: n
            !! Rank of matrix
            integer(ip), intent(in) :: uc_list(:,:)
            !! Uncompressed list/boolean sparse matrix
            integer(ip), intent(in) :: nit(n)
            !! Number of elements for each row of [uc_list]
            type(yale_sparse), intent(out) :: s
            !! Output sparse matrix

            integer(ip), allocatable :: idx(:)
            !! Indices where a certain row begins (RI)
            integer(ip) :: nnz
            !! Number of non-zero elements

            integer(ip) :: i, j

            if(n == 0) then
                s%n = 0
                if(allocated(s%ri)) call mfree('compress_list [ri]', s%ri)
                call mallocate('compress_list [ri]', 1_ip, s%ri)
                s%ri = 1
                if(allocated(s%ci)) call mfree('compress_list [ci]', s%ci)
                call mallocate('compress_list [ci]', 0_ip, s%ci)
                return
            end if

            call mallocate('compress_list [idx]', n, idx)
            idx(1) = 1
            do i=1, n-1
                idx(i+1) = idx(i) + nit(i)
            end do
            nnz = idx(n) + nit(n)

            s%n = n
            if(allocated(s%ri)) call mfree('compress_list [ri]', s%ri)
            call mallocate('compress_list [ri]', n+1, s%ri)
            if(allocated(s%ci)) call mfree('compress_list [ci]', s%ci)
            call mallocate('compress_list [ci]', nnz-1, s%ci)

            !$omp parallel do default(shared) schedule(dynamic) &
            !$omp private(i,j)
            do i=1, n
                s%ri(i)=idx(i)
                do j=0, nit(i)-1
                    s%ci(idx(i)+j) = uc_list(j+1,i)
                end do
            end do
            s%ri(n+1) = nnz
        end subroutine

        subroutine compress_data(s, uc_data, c_data)
            !! Compress the data in uc_data to the same Yale sparse
            !! format described in s

            use mod_memory, only: mallocate, rp

            implicit none

            type(yale_sparse), intent(in) :: s
            !! Input Yale format binary matrix/sparse list
            real(rp), intent(in) :: uc_data(:,:)
            !! Uncompressed data in input
            real(rp), allocatable, intent(out) :: c_data(:)
            !! Compressed data in output
            
            integer(ip) :: nnz, n, i, j0, j

            n = s%n
            nnz = size(s%ci)
            call mallocate('compress_data [c_data]', nnz, c_data)

            !$omp parallel do default(shared) schedule(dynamic) & 
            !$omp private(i,j,j0)
            do i=1, s%n
                j0 = s%ri(i) - 1
                do j=s%ri(i), s%ri(i+1)-1
                    c_data(j) = uc_data(j-j0, i)
                end do
            end do
        end subroutine

end module mod_adjacency_mat
