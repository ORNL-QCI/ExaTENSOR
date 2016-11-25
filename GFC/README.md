GFC: Generic Fortran Containers: Implementation based on OOP Fortran-2003
     that heavily relies on dynamic polymorphism.
Author: Dmitry I. Lyakh (Liakh): liakhdi@ornl.gov

Copyright (C) 2014-2016 Dmitry I. Lyakh (Liakh)
Copyright (C) 2014-2016 Oak Ridge National Laboratory (UT-Battelle)

LICENSE: GNU Lesser General Public License v.3

GENERAL INFORMATION:
 GFC provides a number of heterogeneous containers for higher-level
 algorithm expression in object-oriented Fortran (2003 and later).
 It is still under development, but certain containers have already
 been implemented and tested. The GFC source also contains non-GFC
 based containers developed by me in the past prior to GFC. In future
 they will be fully replaced by their GFC-based reimplementations.
 All GFC-based containers are prefixed with "gfc_".
 Provided (and to be provided) Fortran containers:
 1. Bi-directional Linked List: gfc_list (GFC), lists (legacy);
 2. Tree: gfc_tree (GFC);
 3  Stack: gfc_stack (under development), stack (legacy);
 4. Dictionary (ordered map): gfc_dictionary (under development), dictionary (legacy);
 5. Queue: gfc_queue (under development);
 6. Priority queue: gfc_pri_queue (under development);
 Vector is not implemented due to the power of Fortran arrays,
 except resizability feature.

DESIGN AND FEATURES:
 # A GFC container is a structured collection of objects of any class;
   Objects of different classes can be stored in the same container.
 # Objects can be stored in the container either by value or by reference;
   The same container may have elements stored both ways.
   When storing objects in the container by value, only allocatable components
   of derived types are recursively cloned whereas the pointer components
   are just pointer associated. To change this default behavior, one can
   supply a user-defined generic copy constructor (see interfaces in gfc_base.F90).
 # A GFC subcontainer is a container linked as a part of another container.
   As a consequence, its boundary elements may have outside links.
   In this case, the larger container will be called a composite container.
 # Each container has an associated iterator for scanning over its elements.
   The structure of the container determines the scanning sequence, that is,
   the way the elements of the container are traversed over.
   There are four major ways of scanning a container:
    1) Unpredicated passive scan: The iterator returns each new encountered
       element of the container.
    2) Predicated passive scan: The iterator returns each new encountered
       element of the container that satisfies a certain condition.
    3) Unpredicated active scan: The iterator traverses the entire container
       or its part and applies a user-defined action to each element.
    4) Predicated active scan: The iterator traverses the entire container
       or its part and applies a user-defined action to each element that
       satisfies a certain condition.
   Additionally, active scans allow for a time-limited execution, in which
   the scan is interrupted after some time interval specified by a user.
   Each specific class of containers has its own iterator class because
   of different linkage between the elements of different containers.
   All insertion, deletion, and search operations are done via iterators,
   that is, the iterator methods provide the only possible way of accessing,
   updating, and performing other actions on the associated container.
   All relevant iterator methods are thread-safe, thus enabling
   a parallel execution of container scans (via OpenMP threads).
   However, it is the user responsibility to avoid race conditions
   when updating the value of container elements, that is, the structure
   of the container is protected from races by GFC, but the values of
   container elements are not protected. In the latter case, a thread
   is supposed to acquire an exclusive access to the container element
   if necessary for avoiding race conditions on container values.
   An iterator must not be shared among two or more threads!
 # The container element deletion operation may require a user-defined
   destructor which will release all resources occupied by the object
   stored in that element, unless the object has FINAL methods defined
   (the interface for a user-defined generic destructor is provided in
    gfc_base.F90).

NOTES:
 # This implementation of Generic Fortran Containers heavily relies on
   dynamic polymorphism, thus introducing certain overhead! Also, dynamic
   type inferrence may decrease the efficiency when storing and operating
   on small objects. Thus, this implementation aims at providing a set
   of higher-level abstractions for control logic and other high-level
   operations for which ultimate Flop/s efficiency is not required.
   The use of GFC containers in the inner loop of compute intensive kernels
   is highly discouraged (please resort to plain data, like arrays).
 # Due to the limitations of Fortran class inheritence, public methods
   with a trailing underscore shall NEVER be used by the end user!
