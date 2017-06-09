GFC: Generic Fortran Containers: Implementation based on OOP Fortran-2003
     that heavily relies on dynamic polymorphism.
Author: Dmitry I. Lyakh (Liakh): liakhdi@ornl.gov

Copyright (C) 2014-2017 Dmitry I. Lyakh (Liakh)
Copyright (C) 2014-2017 Oak Ridge National Laboratory (UT-Battelle)

LICENSE: GNU Lesser General Public License v.3 (LGPLv3)

GENERAL INFORMATION:
 GFC provides a number of heterogeneous containers for higher-level
 algorithm expression in object-oriented Fortran (2003 and later).
 It is still under development, but certain containers have already
 been implemented and tested. The GFC source also contains non-GFC
 legacy containers developed by me in the past prior to GFC.
 All GFC-based containers are prefixed with "gfc_".
 Provided (and to be provided) Fortran containers:
 1. Bi-directional Linked List: gfc_list (GFC), lists (legacy);
 2. Tree: gfc_tree (GFC);
 3. Dictionary (ordered map): gfc_dictionary (GFC), dictionary (legacy);
 4. Vector: gfc_vector (GFC);
 5. Vector tree: gfc_vec_tree (GFC): Vector with an imposed tree structure;
 6. Graph: gfc_graph (GFC): Generic graph/hypergraph, a bit slow;
 7. Stack: gfc_stack (GFC, under development), stack (legacy);
 8. Queue: gfc_queue (GFC, under development);
 9. Priority queue: gfc_pri_queue (GFC, under development);
 10. Hash map (unordered map): gfc_hash_map (GFC, under development);

DESIGN AND FEATURES:
 1. A GFC container is a structured collection of scalar objects of any class;
    Objects of different classes can be stored in the same container.
 2. Objects can be stored in the container either by value or by reference;
    The same container may have elements stored both ways.
    When storing objects in the container by value, only allocatable components
    of derived types are recursively cloned whereas the pointer components
    are just pointer associated. To change this default behavior, one can
    supply a user-defined non-member generic copy constructor (see interfaces
    in gfc_base.F90). As a general rule, objects with allocated pointer components
    will have to provide an explicit non-member copy constructor to make sure
    those pointer components are cloned by value instead of pointer association.
 3. A GFC subcontainer is a container linked as a part of another container.
    As a consequence, its boundary elements may have outside links.
    In this case, the larger container will be called the host container.
 4. Each container has an associated iterator for scanning over its elements.
    The structure of the container determines the scanning sequence, that is,
    the path the elements of the container are traversed over.
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
    Active scans require a user-specified action defined either as a stateless
    procedure with a specific generic interface or a Fortran function object,
    both defined in file gfc_base.F90.
    Additionally, active scans allow for a time-limited execution, in which
    the scan is interrupted after some time interval specified by a user.
    Each specific class of containers has its own iterator class because
    of different linkage between the elements of different containers.
    All insertion, deletion, and search operations are done via iterators,
    that is, the iterator methods provide the only possible way of accessing,
    updating, and performing other actions on the associated container.
    The only exception is the duplication of containers as a whole,
    either by value or by reference (e.g., list-to-list, list-to-vector).
    In general multiple iterators may simultaneously be associated with the
    same container and care is required when performing concurrent insertions
    or deletions via different iterators: If the container structure is updated
    via some iterator, all other iterators become invalid until they are
    explicitly reset via the .reset() method.
 5. The container element deletion operation may require a user-defined
    non-member destructor which will properly release all resources occupied
    by the object stored in that element, unless the object has FINAL methods
    defined (the interface for a user-defined non-member generic destructor
    is provided in gfc_base.F90).

 FEATURES UNDER DEVELOPMENT:
 1. Currently SCAN is the only global operation defined on any container and
    provided by the corresponding iterator class. Equipped with a proper user-defined
    function object, it can deliver broad functionality (transforms, reductions, etc.).
    In future, however, other predefined global methods are expected to be added.
 2. Subcontainer feature is not available for all containers, and even when it is
    available it has not been carefully tested yet.
 3. Multi-threading support in future: All relevant iterator methods will be thread-safe,
    thus enabling a parallel execution of container operations (via OpenMP threads).
    However, it is the user responsibility to avoid race conditions
    when updating the values of container elements, that is, the structure
    of the container will be protected from races by GFC, but the values of
    container elements are not protected automatically. A thread is supposed to
    acquire an exclusive access to the container element via GFC lock API,
    if necessary for avoiding race conditions on container values.
    An iterator must never be shared between multiple OpenMP threads!

NOTES:
 1. This implementation of Generic Fortran Containers heavily relies on
    dynamic polymorphism (and RTTI), thus introducing certain overhead! Also,
    dynamic type inferrence may decrease the efficiency when storing and
    operating on small objects. Thus, this implementation aims at providing a set
    of rather expensive higher-level abstractions for control logic and other
    high-level operations for which ultimate numeric Flop/s efficiency is not
    required. The use of GFC containers in the inner loop of compute intensive
    kernels is highly discouraged (please resort to plain data, like arrays).
 2. Due to limitations of Fortran class inheritence, public methods
    with a trailing underscore shall NEVER be used by the end user!
 3. Currently an insertion of an allocated pointer variable into a GFC container
    by reference will result in the loss of the allocated status of the stored
    pointer in case it was allocated. Consequently, the stored pointer cannot
    later be used for deallocating the corresponding target, although some
    compilers ignore this.
