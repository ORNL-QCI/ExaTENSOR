GFC: Generic Fortran Containers: Implementation based on OOP Fortran-2003
     that heavily relies on dynamic polymorphism and RTTI.
Author: Dmitry I. Lyakh (Liakh): liakhdi@ornl.gov

Copyright (C) 2014-2022 Dmitry I. Lyakh (Liakh)
Copyright (C) 2014-2022 Oak Ridge National Laboratory (UT-Battelle)

LICENSE: BSD 3-Clause

ISSUES:
 GFC heavily uses modern Fortran-2003 features which are
 implemented by a number of Fortran compilers, including
 GCC 5.3 and higher, Intel 17 and higher, IBM XL 15 and higher,
 PGI 17 and higher, and Cray. Unfortunately almost all above
 compilers may still have bugs in the implementation of Fortran-2003
 features. GCC/8.0.0 has fixed a number of critical bugs (special thanks
 to Paul Thomas and other GCC developers). Other compilers to try are:
 Intel 18+, IBM XL 15+, PGI 17.10+, Cray 8.6.x+.

GENERAL INFORMATION:
 GFC provides a number of heterogeneous containers for higher-level
 algorithm expression in object-oriented Fortran (2003 and later).
 It is still under development, but certain containers have already
 been implemented and tested. The GFC source also contains non-GFC
 legacy containers developed by me in the past prior to GFC.
 All GFC-based containers are prefixed with "gfc_".
 Provided (and to be provided) Fortran containers:
 0. Bank: gfc_bank (GFC): A pool of reusable objects of any type;
 1. Bi-directional Linked List: gfc_list (GFC), lists (legacy);
 2. Tree: gfc_tree (GFC);
 3. Dictionary (ordered map): gfc_dictionary (GFC), dictionary (legacy);
 4. Vector: gfc_vector (GFC);
 5. Vector tree: gfc_vec_tree (GFC): Vector with an imposed tree structure;
 6. Graph: gfc_graph (GFC): Generic graph/hypergraph (slow);
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
    either define a user-defined assignment operator in the relevant classes or
    supply a user-defined non-member generic GFC copy constructor (see interfaces
    in gfc_base.F90). As a general rule, objects with allocated pointer components
    will have to provide either of the above to make sure those pointer components
    are cloned by value instead of mere pointer association.
 3. A GFC subcontainer is a container linked as a part of another container.
    As a consequence, its boundary elements may have outside links.
    In this case, the encompassing container will be called the host container.
 4. Each container has an associated iterator for traversing its elements.
    The structure of the container determines the scanning sequence, that is,
    the path the elements of the container are traversed over.
    There are four major ways of container scanning:
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
    procedure with a specific generic interface or a GFC function object,
    both defined in file gfc_base.F90.
    Additionally, active scans allow for a time-limited execution, in which
    the scan is interrupted after some time interval specified by a user.
    Each specific class of containers has its own iterator class because
    of different linkage between the elements of different containers.
    All access, insertion, deletion, and search operations are done via iterators,
    that is, the iterator methods provide the only possible way of accessing,
    updating, and performing other actions on the associated container.
    The only exception is the duplication of containers as a whole,
    either by value or by reference (e.g., list-to-list, list-to-vector).
    In general multiple iterators may simultaneously be associated with the
    same container and care is required when performing concurrent insertions
    or deletions via different iterators: If the container structure is updated
    via some iterator (insert/delete), all other iterators become invalid until
    they are explicitly reset via the .reset() method (TO BE RELAXED IN FUTURE).
 5. The container element deletion operation may require a user-defined
    non-member GFC destructor which will properly release all resources occupied
    by the object stored in that element, unless the object has FINAL methods
    defined (the interface for a user-defined non-member generic GFC destructor
    is provided in gfc_base.F90).

FEATURES UNDER DEVELOPMENT:
 1. Currently SCAN is the only global operation defined on any container and
    provided by the corresponding iterator class. Equipped with a proper user-defined
    function object, it can deliver broad functionality (transforms, reductions, etc.).
    In future, however, other predefined global methods are expected to be added.
 2. Subcontainer feature is not available for all containers and, even when it is
    available, it has not been carefully tested yet.
 3. Multi-threading support in future: All relevant iterator methods will be thread-safe,
    thus enabling a parallel execution of container operations (via OpenMP threads).
    However, it is the user responsibility to avoid race conditions
    when updating the values of container elements, that is, the structure
    of the container will be protected from races by GFC, but the values of
    container elements are not protected automatically. A thread is supposed to
    acquire an exclusive access to the container element via GFC lock API
    if necessary for avoiding race conditions on container values.
    An iterator must never be shared between multiple OpenMP threads!

NOTES:
 1. This implementation of Generic Fortran Containers heavily relies on
    dynamic polymorphism (and RTTI), thus introducing certain overhead! Also,
    dynamic type inferrence may decrease the efficiency when storing and
    operating on small objects. Thus, this implementation aims at providing a set
    of rather EXPENSIVE higher-level abstractions for control logic and other
    high-level operations for which ultimate performance is not required.
    The use of GFC containers in the inner loop of compute intensive
    kernels is highly discouraged (please resort to plain data, like arrays).
 2. Due to the limitations of Fortran class inheritence, public methods
    with a trailing underscore shall NEVER be used by the end user!
 3. Currently an insertion of an allocated pointer variable into a GFC container
    by reference will result in the loss of the allocated status of the stored
    pointer in case it was allocated. Consequently, the stored pointer cannot
    later be used for deallocating the corresponding target, although some
    compilers ignore this. In other words, storing allocated pointers by
    reference assumes no resource ownership by GFC.

TO DO (developers only):
 1. Thread-safety and multi-threading:
    a) All iterator moves must be done via explicit iterator move methods,
       including the most general .jump() method. This will allow reference
       counting for the number of iterator instances currently associated
       with each container element. Iterator initialization/release and
       all iterator moves must be atomic, where the atomicity is ensured
       by atomicity of the reference count update in each container element.
       Container elements with a non-zero reference count cannot be deleted.
    b) The default (suboptimal) race protection for container structure
       updates should be implemented via a per-container-instance lock.
       While a container is locked, no other action can be performed on it,
       except accessing/mutating the values of its elements via existing
       iterators (here the race conditions are user responsibility).
       In particular, a locked container will not accept iterator moves.
    c) All iterators should provide a .sync() method to synchronize their
       current state and accomodate recent actions performed on the same
       container via other iterators.
 2. Type-safety:
    a) Provide a generic method for each container class called
       .register_type(character_type_name,empty_object_of_some_type):
       It will clone empty_object_of_some_type and store it inside the
       container to test all subsequent insertions to be of the same
       type/class, if needed. The character_type_name can later be used
       to retrieve the information on the stored type.
