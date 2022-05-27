module adjacency_matrix
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
    public :: adj_mat_from_conn, build_conn_upto_n

    contains

        subroutine matcpy(f, t)
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
        end subroutine matcpy

        subroutine matfree(m)
            !! Deallocate a boolean sparse matrix (in Yale format) 
            implicit none
            
            type(yale_sparse), intent(inout) :: m !! Matrix to free

            deallocate(m%ri)
            deallocate(m%ci)
            m%n=0_ip
        end subroutine matfree

        subroutine adj_mat_from_conn(i12, n12, sparse)
            !! Create adjacency matrix \(\mathbb C_1\) from connectivity lists.   
            !! Array i12 and n12 contain the connectivity list in the following
            !! format: i12(0:n(j),j) contains the index of all the atoms connected
            !! to atom with index j.

            implicit none

            integer(ip), intent(in) :: i12(:,:)
            !! Indices of connected atoms for each atom in the molecule
            integer(ip), intent(in) :: n12(:)
            !! Number of connected atoms for each atom in the molecule
            type(yale_sparse), intent(out) :: sparse
            !! Adjacency matrix in Yale format (\(\mathbb C_1\))

            integer(ip) :: i, j, nnz

            nnz = sum(n12)
            
            sparse%n = size(i12, 2)
            allocate(sparse%ri(sparse%n+1))
            allocate(sparse%ci(nnz))

            ! Compress adjacency matrix
            sparse%ri(1) = 1
            do i = 1, sparse%n
                sparse%ri(i+1) = sparse%ri(i)

                do j = 1, size(i12, 1)
                    if(i12(j,i) /= 0) then
                        sparse%ci(sparse%ri(i+1)) = i12(j,i)
                        sparse%ri(i+1) = sparse%ri(i+1) + 1
                    end if
                end do
            end do

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
                                ! write(*, *) "REALLOCATING RES"
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
            

            do i = 1, n
                res%ci(i) = i
                res%ri(i) = i
            end do
            res%ri(n+1) = n+1

        end subroutine sparse_identity

        subroutine build_conn_upto_n(adj, n, res)
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

            integer(ip) :: i
            type(yale_sparse) :: tmp, id

            allocate(res(n))

            call matcpy(adj, res(1))

            do i=2, n
                call mat_mult(res(i-1), res(1), res(i))
                !call matcpy(res(i-1), res(i))
                call mat_andnot(res(i), res(i-1), tmp)
                if(i == 2) then
                    call sparse_identity(adj%n, id)
                    call mat_andnot(tmp, id, res(i))
                    call matfree(id)
                else
                    call mat_andnot(tmp, res(i-2), res(i))
                end if
            end do
            call matfree(tmp)
        end subroutine build_conn_upto_n
end module adjacency_matrix
