      module vectormodule
         use goverlapmodule
         implicit none
         type :: Vector
            type(GOverlap), allocatable, dimension(:) :: data
            integer :: capacity
            integer :: length
            contains
               procedure, pass(vec) :: create_vector 
               procedure, pass(vec) :: push_back 
               procedure, pass(vec) :: pop_back 
               procedure, pass(vec) :: clear
               procedure, pass(vec) :: vector_size
               procedure, pass(vec) :: empty 
               procedure, pass(vec) :: get_element
               procedure, pass(vec) :: set_element
               procedure, pass(vec) :: sort
         end type Vector

         contains

         subroutine create_vector(vec)
            class(Vector), intent(out) :: vec
            vec%capacity = 10  ! Initial capacity
            vec%length = 0
            allocate(vec%data(0:vec%capacity-1))  ! Start index at 0
         end subroutine create_vector

         subroutine push_back(vec, value)
            class(Vector), intent(inout) :: vec
            type(GOverlap), intent(in) :: value
            type(GOverlap), dimension(:), allocatable :: newData
            integer :: i
            if (vec%length == 0) then
               call vec%create_vector()
            else if (vec%length == vec%capacity) then
               vec%capacity = vec%capacity * 2  ! Double the capacity
               allocate(newData(0:vec%capacity-1))
               do i = 0, vec%length - 1
                  newData(i) = vec%data(i)
               end do
               deallocate(vec%data) ! Deallocate the old array
               vec%data = newData   ! Point the data pointer to the new array
            end if
            vec%data(vec%length) = value
            vec%length = vec%length + 1
         end subroutine push_back

         subroutine pop_back(vec)
            class(Vector), intent(inout) :: vec
            if (vec%length > 0) then
               vec%length = vec%length - 1
            end if
         end subroutine pop_back

         subroutine clear(vec)
            class(Vector), intent(inout) :: vec
            deallocate(vec%data)
            vec%capacity = 0
            vec%length = 0
         end subroutine clear

         function vector_size(vec) result(n)
            class(Vector), intent(in) :: vec
            integer :: n
            n = vec%length
         end function vector_size

         function empty(vec) result(isEmpty)
            class(Vector), intent(in) :: vec
            logical :: isEmpty
            isEmpty = (vec%length == 0)
         end function empty

         function get_element(vec, index) result(value)
            class(Vector), intent(in) :: vec
            integer, intent(in) :: index
            type(GOverlap) :: value
            if (index >= 0 .and. index <= vec%length - 1) then
               value = vec%data(index)
            end if
         end function get_element

         subroutine set_element(vec, index, value)
            class(Vector), intent(inout) :: vec
            integer, intent(in) :: index
            type(GOverlap), intent(in) :: value
            if (index >= 0 .and. index <= vec%length - 1) then
               vec%data(index) = value
            end if
         end subroutine set_element

         subroutine sort(vec)
            class(Vector), intent(inout) :: vec
            integer :: low, high
            low = 0
            high = vec%length - 1
            call quicksort_overlaps(vec, low, high)
         end subroutine sort

         recursive subroutine quicksort_overlaps(vec, low, high)
            class(Vector), intent(inout) :: vec
            integer, intent(in) :: low, high
            integer :: pivot
            if (low < high) then
               pivot = partition_overlaps(vec, low, high)
               call quicksort_overlaps(vec, low, pivot - 1)
               call quicksort_overlaps(vec, pivot + 1, high)
            end if
         end subroutine quicksort_overlaps

         function partition_overlaps(vec, low, high) result(pivot)
            class(Vector), intent(inout) :: vec
            integer, intent(in) :: low, high
            type(GOverlap) :: temp1, temp2
            integer :: i, j, pivot
            integer :: pivot_index
            pivot_index = high
            i = low - 1
            do j = low, high - 1
               temp1 = vec%get_element(j)
               temp2 = vec%get_element(pivot_index)
               if (temp1%volume < temp2%volume) then
                  i = i + 1
                  temp1 = vec%get_element(i)
                  temp2 = vec%get_element(j)
                  call vec%set_element(i, temp2)
                  call vec%set_element(j, temp1)
               end if
            end do
            temp1 = vec%get_element(i + 1)
            call vec%set_element(i + 1, vec%get_element(pivot_index))
            call vec%set_element(pivot_index, temp1)
            pivot = i + 1
         end function partition_overlaps

      end module vectormodule
