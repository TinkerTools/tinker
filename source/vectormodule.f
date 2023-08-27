c
c
c     ##################################################################
c     ##  COPYRIGHT (C) 2023 by MKJ Chung, MJ Schnieders & JW Ponder  ##
c     ##                     All Rights Reserved                      ##
c     ##################################################################
c
c     ##############################################################
c     ##                                                          ##
c     ##  module vectormodule  --  dynamic vector implementation  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     data       vector of GOverlap objects
c     capacity   capacity of vector
c     length     length of vector
c
c
      module vectormodule
      use goverlapmodule
      implicit none
      type :: Vector
      type(GOverlap), allocatable :: data(:)
      integer capacity
      integer length
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
c
c
c     ############################################################
c     ##                                                        ##
c     ##  subroutine create_vector  --  allocate inital vector  ##
c     ##                                                        ##
c     ############################################################
c
c
c     "create_vector" initializes a vector of length 0 and capactiy 10
c
c
      subroutine create_vector(vec)
      class(Vector) vec
c
c
c     set initial capacity and length
c
      vec%capacity = 10
      vec%length = 0
c
c     allocate vector whose index starts at 0
c
      allocate(vec%data(0:vec%capacity-1))
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine push_back  --  push element to end of vector  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "push_back" pushes back an object to the end of the vector
c
c
      subroutine push_back(vec, value)
      class(Vector) vec
      type(GOverlap) value
      type(GOverlap), allocatable :: newData(:)
      integer i
c
c
c     initialize vector if not initialized
c
      if (vec%capacity == 0) then
         call vec%create_vector()
c
c     double capacity if length = capacity
c
      else if (vec%length == vec%capacity) then
         vec%capacity = vec%capacity * 2
         allocate(newData(0:vec%capacity-1))
         do i = 0, vec%length - 1
            newData(i) = vec%data(i)
         end do
c
c     deallocate old array, point data pointer to new array
c
         deallocate(vec%data)
         vec%data = newData
      end if
c
c     append element to end of vector
c
      vec%data(vec%length) = value
      vec%length = vec%length + 1
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine pop_back  --  removes last element of vector  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "pop_back" removes element at the end of vector
c
c
      subroutine pop_back(vec)
      class(Vector) vec
      if (vec%length > 0) then
         vec%length = vec%length - 1
      end if
      return
      end
c
c
c     ###############################################
c     ##                                           ##
c     ##  subroutine clear  --  deallocate vector  ##
c     ##                                           ##
c     ###############################################
c
c
c     "clear" deallocates vector and sets length and capacity to 0
c
c
      subroutine clear(vec)
      class(Vector) vec
      deallocate(vec%data)
      vec%capacity = 0
      vec%length = 0
      return
      end
c
c
c     ####################################################
c     ##                                                ##
c     ##  function vector_size  --  return vector size  ##
c     ##                                                ##
c     ####################################################
c
c
c     "vector_size" returns the size of vector
c
c
      function vector_size(vec) result(n)
      class(Vector) vec
      integer n
      n = vec%length
      return
      end
c
c
c     ####################################################
c     ##                                                ##
c     ##  function empty  --  check if vector is empty  ##
c     ##                                                ##
c     ####################################################
c
c
c     "empty" returns true if vector has no elements
c
c
      function empty(vec) result(isEmpty)
      class(Vector) vec
      logical isEmpty
      isEmpty = (vec%length == 0)
      return
      end
c
c
c     ##########################################################
c     ##                                                      ##
c     ##  function get_element  --  getter method for vector  ##
c     ##                                                      ##
c     ##########################################################
c
c
c     "get_element" returns the element at the vector index
c
c
      function get_element(vec, index) result(value)
      class(Vector) vec
      type(GOverlap) value
      integer index
      if (index >= 0 .and. index <= vec%length - 1) then
         value = vec%data(index)
      end if
      return
      end
c
c
c     ############################################################
c     ##                                                        ##
c     ##  subroutine set_element  --  setter method for vector  ##
c     ##                                                        ##
c     ############################################################
c
c
c     "set_element" sets the element at the vector index
c
c
      subroutine set_element(vec, index, value)
      class(Vector) vec
      type(GOverlap) value
      integer index
      if (index >= 0 .and. index <= vec%length - 1) then
         vec%data(index) = value
      end if
      return
      end
c
c
c     ###################################################
c     ##                                               ##
c     ##  subroutine sort  --  sort method for vector  ##
c     ##                                               ##
c     ###################################################
c
c
c     "sort" uses the quicksort algorithm to sort the vector
c
c
      subroutine sort(vec)
      class(Vector) vec
      integer low, high
      low = 0
      high = vec%length - 1
      call quicksort_overlaps(vec, low, high)
      return
      end
c
c
c     ###########################################################
c     ##                                                       ##
c     ##  subroutine quicksort_overlaps  --  quicksort method  ##
c     ##                                                       ##
c     ###########################################################
c
c
c     "quicksort_overlaps" recursively sorts vector using quicksort
c
c
      recursive subroutine quicksort_overlaps(vec, low, high)
      class(Vector) vec
      integer low, high
      integer pivot
      if (low < high) then
         pivot = partition_overlaps(vec, low, high)
         call quicksort_overlaps(vec, low, pivot - 1)
         call quicksort_overlaps(vec, pivot + 1, high)
      end if
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  function partition_overlaps  --  partition in quicksort  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "partition_overlaps" is the partitioning step of the quicksort
c     algorithm, where elements are rearranged so that elements with
c     smaller volumes are placed before elements with larger volumes
c
c
      function partition_overlaps(vec, low, high) result(pivot)
      class(Vector) vec
      type(GOverlap) temp1, temp2
      integer low, high
      integer i, j, pivot
      integer pivot_index
      pivot_index = high
      i = low - 1
      do j = low, high - 1
         temp1 = vec%get_element(j)
         temp2 = vec%get_element(pivot_index)
         if (temp1%volume <= temp2%volume) then
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
      return
      end
      end module vectormodule
