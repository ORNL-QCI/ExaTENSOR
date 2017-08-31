!Domain-specific virtual processor (DSVP): Abstract base module.
!AUTHOR: Dmitry I. Lyakh (Liakh): quant4me@gmail.com
!REVISION: 2017/08/31

!Copyright (C) 2014-2017 Dmitry I. Lyakh (Liakh)
!Copyright (C) 2014-2017 Oak Ridge National Laboratory (UT-Battelle)

!This file is part of ExaTensor.

!ExaTensor is free software: you can redistribute it and/or modify
!it under the terms of the GNU Lesser General Public License as published
!by the Free Software Foundation, either version 3 of the License, or
!(at your option) any later version.

!ExaTensor is distributed in the hope that it will be useful,
!but WITHOUT ANY WARRANTY; without even the implied warranty of
!MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!GNU Lesser General Public License for more details.

!You should have received a copy of the GNU Lesser General Public License
!along with ExaTensor. If not, see <http://www.gnu.org/licenses/>.

!NOTES:
! # A domain-specific virtual processor (DSVP) processes domain-specific instructions
!   operating on domain-specific operands. Domain-specific operands encapsulate
!   domain-specific data, making it processible by domain-specific instructions
!   which encapsulate domain-specific operations. These abstract classes are intended
!   for a further specialization to concrete domain-specific classes via inheritance.
! # The abstract domain-specific processor classes provide basic instruction processing
!   primitives, relevant deferred interfaces, and other bookkeeping:
!   a) Domain-specific operand class:
!      * Acquire local resources;
!      * Prefetch operand;
!      * Upload operand;
!      * Sync operand prefetching/uploading;
!      * Release local resources.
!   b) Domain-specific instruction:
!      * Acquire resource: Acquires local resources for instruction operands;
!      * Prefetch input: Starts prefetching (remote) input operands;
!      * Sync prefetch: Synchronizes input prefetch;
!      * Execute: Starts execution of the domain-specific instruction;
!      * Sync execution: Synchronizes instruction execution (completion or error);
!      * Upload output: Starts uploading (remote) output operands;
!      * Sync upload: Synchronizes output upload;
!      * Release resource: Releases local resources occupied by instruction operands;
!      * Encode: Encodes a domain-specific instruction into the bytecode packet.
!   c) Domain-specific virtual processor:
!      * Start;
!      * Shutdown;
!      * Decode: Unpacks the instruction bytecode into a domain-specific instruction.
! # A concrete domain-specific virtual processor implementation derives from the abstract classes
!   specified in this module. The derivative concrete classes have to provide implementation
!   of all deferred type-bound procedures. The concrete domain-specific virtual processor class
!   must define the microcode for each concrete domain-specific instruction and override
!   the .decode() method. The .decode() method is supposed to unpack the instruction code
!   from a raw bytecode packet (incoming instruction) and, based on that code, execute the instruction
!   setup microcode which will associate all dynamic bindings (dynamic methods) of the instruction object
!   as well as unpack the instruction control field and instruction operands. After a domain-specific
!   instruction has been succesfully decoded, it is ready for use in the domain-specific virtual processor
!   pipeline implemented by a concrete domain-specific virtual processor class composed from the concrete
!   domain-specific virtual unit (DSVU) classes: DSVP consists of one or more DSVU, for example,
!   decoder DSVU, encoder DSVU, caching DSVU, CPU executing DSVU, GPU executing DSVU, etc.

       module dsvp_base
        use dil_basic
        use timers
        use pack_prim, only: obj_pack_t
        use gfc_base !contains OpenMP also
        use gfc_list
        implicit none
        private
!PARAMETERS:
 !Basic:
        integer(INTD), private:: CONS_OUT=6 !default output for this module
        integer(INTD), private:: DEBUG=1    !debugging mode
        logical, private:: VERBOSE=.TRUE.   !verbosity for errors
 !Error codes:
        integer(INTD), parameter, public:: DSVP_SUCCESS=SUCCESS        !success
        integer(INTD), parameter, public:: DSVP_ERROR=GENERIC_ERROR    !generic error code
        integer(INTD), parameter, public:: DSVP_ERR_INVALID_ARGS=-1    !invalid arguments passed
        integer(INTD), parameter, public:: DSVP_ERR_INVALID_REQ=-2     !invalid request
        integer(INTD), parameter, public:: DSVP_ERR_MEM_ALLOC_FAIL=-3  !failed memory allocation
        integer(INTD), parameter, public:: DSVP_ERR_MEM_FREE_FAIL=-4   !failed memory deallocation
        integer(INTD), parameter, public:: DSVP_ERR_BROKEN_OBJ=-5      !broken object
        integer(INTD), parameter, public:: DSVP_ERR_UNABLE_COMPLETE=-6 !unable to complete
        integer(INTD), parameter, public:: DSVP_ERR_RSC_EXCEEDED=-7    !resource exceeded
 !DSVP status (negative numbers are specific error codes):
        integer(INTD), parameter, public:: DSVP_STAT_OFF=0          !DSVP is off (either not initialized or turned off)
        integer(INTD), parameter, public:: DSVP_STAT_ON=1           !DSVP has been initialized and is active now
        integer(INTD), parameter, public:: DSVP_STAT_ERR=-1         !DSVP encountered an error (generic)
 !DSVP specific kind (valid specific DSVP kinds must be non-negative):
        integer(INTD), parameter, public:: DSVP_NO_KIND=-1          !no specficic kind
 !Domain-specific operand:
  !Operand status:
        integer(INTD), parameter, public:: DS_OPRND_EMPTY=0         !empty operand
        integer(INTD), parameter, public:: DS_OPRND_DEFINED=1       !defined operand, but the data has not been delivered yet
        integer(INTD), parameter, public:: DS_OPRND_PRESENT=2       !defined operand with the data locally present (delivered)
  !Data communication status:
        integer(INTD), parameter, public:: DS_OPRND_NO_COMM=0       !no pending communication on the domain-specific operand
        integer(INTD), parameter, public:: DS_OPRND_FETCHING=1      !domain-specific operand is being fetched
        integer(INTD), parameter, public:: DS_OPRND_UPLOADING=-1    !domain-specific operand is being uploaded
 !Domain-specific instruction:
  !Instruction code (valid codes must be non-negative):
        integer(INTD), parameter, public:: DS_INSTR_NOOP=-1         !no operation (all valid instruction codes are non-negative)
  !Instruction status (instruction pipeline stages):
        integer(INTD), parameter, public:: DS_INSTR_EMPTY=0         !empty instruction
        integer(INTD), parameter, public:: DS_INSTR_NEW=1           !new (freshly decoded) instruction
        integer(INTD), parameter, public:: DS_INSTR_RSC_WAIT=2      !waiting for resource allocation required by instruction operands (local resources are acquired)
        integer(INTD), parameter, public:: DS_INSTR_INPUT_WAIT=3    !waiting for the (remote) input operands to be delivered
        integer(INTD), parameter, public:: DS_INSTR_READY_TO_EXEC=4 !instruction is ready for execution
        integer(INTD), parameter, public:: DS_INSTR_SCHEDULED=5     !instruction has been scheduled to a specific computing unit (specific queue)
        integer(INTD), parameter, public:: DS_INSTR_ISSUED=6        !instruction has been issued for execution on a specific computing unit
        integer(INTD), parameter, public:: DS_INSTR_COMPLETED=7     !instruction has completed execution (either successfully or with an error)
        integer(INTD), parameter, public:: DS_INSTR_OUTPUT_WAIT=8   !waiting for the (remote) output operands to be uploaded back
        integer(INTD), parameter, public:: DS_INSTR_RETIRED=9       !instruction retired (all temporary resources have been released)
!DERIVED TYPES:
 !Domain-specific resource (normally local):
        type, abstract, public:: ds_resrc_t
         contains
          procedure(ds_resrc_query_i), deferred, public:: is_empty  !returns TRUE if the resource is empty (unacquired), FALSE otherwise
        end type ds_resrc_t
 !Domain-specific operand (will contain domain-specific data to be processed by DSVP):
        type, abstract, public:: ds_oprnd_t
         integer(INTD), private:: stat=DS_OPRND_EMPTY       !current status of the domain-specific operand: {DS_OPRND_EMPTY,DS_OPRND_DEFINED,DS_OPRND_PRESENT}
         integer(INTD), private:: in_route=DS_OPRND_NO_COMM !communication status: {DS_OPRND_NO_COMM,DS_OPRND_FETCHING,DS_OPRND_UPLOADING}
         contains
          procedure(ds_oprnd_query_i), deferred, public:: is_remote  !checks whether the domain-specific operand is local or remote
          procedure(ds_oprnd_self_i), deferred, public:: acquire_rsc !explicitly acquires local resources for the domain-specific operand
          procedure(ds_oprnd_self_i), deferred, public:: prefetch    !starts prefetching a remote domain-specific operand (acquires local resources!)
          procedure(ds_oprnd_self_i), deferred, public:: upload      !starts uploading the domain-specific operand to its remote location
          procedure(ds_oprnd_sync_i), deferred, public:: sync        !synchronizes the currently pending communication on the domain-specific operand (either test or wait)
          procedure(ds_oprnd_self_i), deferred, public:: release     !destroys the present local copy of the domain-specific operand (releases local resources!), but the operand will stay defined
          procedure(ds_oprnd_self_i), deferred, public:: destruct    !performs complete destruction back to an empty (undefined) state
          procedure, public:: is_active=>DSOprndIsActive               !returns TRUE if the domain-specific operand is active (defined and maybe present)
          procedure, public:: is_present=>DSOprndIsPresent             !returns TRUE if the domain-specific operand is locally available (present)
          procedure, public:: mark_active=>DSOprndMarkActive           !marks the domain-specific operand active (defined)
          procedure, public:: mark_empty=>DSOprndMarkEmpty             !marks the domain-specific operand inactive (empty), local resources are released
          procedure, public:: mark_delivered=>DSOprndMarkDelivered     !marks the domain-specific operand locally available (present)
          procedure, public:: mark_undelivered=>DSOprndMarkUndelivered !marks the domain-specific operand locally unavailable (but defined), local resources are released
          procedure, public:: get_comm_stat=>DSOprndGetCommStat        !gets the communication status
          procedure, public:: set_comm_stat=>DSOprndSetCommStat        !sets the communication status
        end type ds_oprnd_t
 !Wrapped reference to a domain-specific operand:
        type, private:: ds_oprnd_ref_t
         class(ds_oprnd_t), pointer, private:: oprnd_ref=>NULL() !non-owning pointer (reference) to a domain-specific operand
        end type ds_oprnd_ref_t
 !Domain-specific instruction control field:
        type, abstract, public:: ds_instr_ctrl_t
         contains
          procedure(ds_instr_ctrl_pack_i), deferred, public:: pack     !packs the domain-specific instruction control into a plain byte packet
          procedure(ds_instr_ctrl_unpack_i), deferred, public:: unpack !unpacks the domain-specific instruction control from a plain byte packet
        end type ds_instr_ctrl_t
 !Domain-specific instruction (realization of a domain-specific operation for DSVP):
        type, abstract, public:: ds_instr_t
         integer(INTD), private:: code=DS_INSTR_NOOP                !all valid instruction codes are non-negative (negative means no operation)
         integer(INTD), private:: num_oprnds=0                      !number of domain-specific operands: Numeration:[0,1,2,3,...]
         integer(INTD), private:: stat=DS_INSTR_EMPTY               !status of the domain-specific instruction in the processing pipeline
         integer(INTD), private:: error_code=DSVP_SUCCESS           !error code (success:DSVP_SUCCESS)
         class(ds_instr_ctrl_t), pointer, private:: control=>NULL() !instruction control field: set up by the DECODE procedure
         type(ds_oprnd_ref_t), allocatable, private:: operand(:)    !domain-specific operands (wrapped pointers): set up by the DECODE procedure
         !dynamic methods (dynamically bound by the DECODE procedure):
         procedure(ds_instr_self_i), pass(this), pointer, public:: acquire_resource=>NULL() !acquires local resources for instruction operands: dynamic binding set by decode
         procedure(ds_instr_self_i), pass(this), pointer, public:: prefetch_input=>NULL()   !starts prefetching input operands: dynamic binding set by decode
         procedure(ds_instr_sync_i), pass(this), pointer, public:: sync_prefetch=>NULL()    !synchronizes the input prefetch (either test or wait): dynamic binding set by decode
         procedure(ds_instr_self_i), pass(this), pointer, public:: execute=>NULL()          !executes the domain-specific instruction: dynamic binding set by decode
         procedure(ds_instr_sync_i), pass(this), pointer, public:: sync_execution=>NULL()   !synchronizes the execution (either test or wait): dynamic binding set by decode
         procedure(ds_instr_self_i), pass(this), pointer, public:: upload_output=>NULL()    !starts uploading the output: dynamic binding set by decode
         procedure(ds_instr_sync_i), pass(this), pointer, public:: sync_upload=>NULL()      !synchronizes the output upload (either test or wait): dynamic binding set by decode
         procedure(ds_instr_self_i), pass(this), pointer, public:: release_resource=>NULL() !releases local resources occupied by instruction operands: dynamic binding set by decode
         contains
          procedure(ds_instr_encode_i), deferred, public:: encode       !encoding procedure: Packs the domain-specific instruction into a raw byte packet (bytecode)
          procedure, public:: is_empty=>DSInstrIsEmpty                  !returns TRUE if the domain-specific instruction is empty
          procedure, public:: is_retired=>DSInstrIsRetired              !returns TRUE if the domain-specific instruction is retired
          procedure, public:: is_active=>DSInstrIsActive                !returns TRUE if the domain-specific instruction is active (defined)
          procedure, public:: get_code=>DSInstrGetCode                  !returns the instruction code
          procedure, public:: set_code=>DSInstrSetCode                  !sets the instruction code
          procedure, public:: get_status=>DSInstrGetStatus              !returns the current status of the domain-specific instruction and error code
          procedure, public:: set_status=>DSInstrSetStatus              !sets the status of the domain-specific instruction and error code
          procedure, public:: get_control=>DSInstrGetControl            !returns the pointer to the instruction control field
          procedure, public:: set_control=>DSInstrSetControl            !associates the control field pointer
          procedure, public:: free_control=>DSInstrFreeControl          !frees the control field pointer
          procedure, public:: alloc_operands=>DSInstrAllocOperands      !allocates instruction operand placeholders
          procedure, public:: dealloc_operands=>DSInstrDeallocOperands  !deallocates instruction operands
          procedure, public:: get_operand=>DSInstrGetOperand            !returns a pointer to a specific operand of the domain-specific instruction
          procedure, public:: set_operand=>DSInstrSetOperand            !associates a specific operand of the domain-specific instruction with its target
          procedure, public:: free_operand=>DSInstrFreeOperand          !frees a specific instruction operand
          procedure, public:: get_num_operands=>DSInstrGetNumOperands   !returns the number of operands in the domain-specific instruction
          procedure, public:: all_set=>DSInstrAllSet                    !returns TRUE if all operands and control are set
          procedure, public:: activate=>DSInstrActivate                 !activates the domain-specific instruction in the final stage of decoding (sets up dynamic microcode bindings, opcode and status)
          procedure, public:: terminate=>DSInstrTerminate               !terminates the normal instruction execution workflow, but leaves instruction defined (retired)
          procedure, public:: clean=>DSInstrClean                       !resets the domain-specific instruction to an empty state (after it has been retired)
        end type ds_instr_t
 !Domain-specific microcode binding (for a DS instruction, see above):
        type, public:: ds_microcode_t
         procedure(ds_instr_self_i), nopass, pointer, public:: acquire_resource=>NULL() !acquires local resources for instruction operands
         procedure(ds_instr_self_i), nopass, pointer, public:: prefetch_input=>NULL()   !starts prefetching input operands
         procedure(ds_instr_sync_i), nopass, pointer, public:: sync_prefetch=>NULL()    !synchronizes the input prefetch (either test or wait)
         procedure(ds_instr_self_i), nopass, pointer, public:: execute=>NULL()          !executes the domain-specific instruction
         procedure(ds_instr_sync_i), nopass, pointer, public:: sync_execution=>NULL()   !synchronizes the execution (either test or wait)
         procedure(ds_instr_self_i), nopass, pointer, public:: upload_output=>NULL()    !starts uploading the output operand
         procedure(ds_instr_sync_i), nopass, pointer, public:: sync_upload=>NULL()      !synchronizes the output upload (either test or wait)
         procedure(ds_instr_self_i), nopass, pointer, public:: release_resource=>NULL() !releases local resources occupied by instruction operands
         contains
          procedure, public:: reset=>DSMicrocodeReset !resets microcode binding to NULL
        end type ds_microcode_t
 !DSVP/DSVU configuration:
        type, abstract, public:: dsv_conf_t
        end type dsv_conf_t
 !Domain-specific virtual unit port:
        type, public:: ds_unit_port_t
         logical, private:: locked=.FALSE.              !lock for updates
         type(list_bi_t), private:: queue               !queue of incoming DS instructions stored by reference
         type(list_iter_t), private:: iqueue            !queue iterator
         contains
          procedure, public:: accept=>DSUnitPortAccept  !accepts new DS instructions from other DS units
          procedure, public:: absorb=>DSUnitPortAbsorb  !absorbs new DS instructions from the DS unit port into the DS unit queue
          procedure, public:: free=>DSUnitPortFree      !releases all DS instructions and resets everything
          final:: ds_unit_port_dtor                     !dtor
        end type ds_unit_port_t
 !Domain-specific virtual unit (DSVU):
        type, abstract, public:: ds_unit_t
         integer(INTD), private:: id=-1                     !unique ID of the DS unit: [0..max]
         integer(INTD), private:: error_code=DSVP_SUCCESS   !error code
         class(dsvp_t), pointer, private:: dsvp_p=>NULL()   !non-owning pointer to the DSVP the DS unit is part of
         type(ds_unit_port_t), private:: port               !DS unit port (for incoming DS instructions from other DS units)
         type(list_bi_t), private:: queue                   !queue of the currently processed DS instructions stored by reference
         type(list_iter_t), private:: iqueue                !queue iterator
         contains
          procedure(ds_unit_ctor_i), deferred, public:: configure !configures the DS unit
          procedure(ds_unit_self_i), deferred, public:: start     !starts, lives and stops the DS unit (the corresponding thread will run here until termination)
          procedure(ds_unit_self_i), deferred, public:: shutdown  !stops the DS unit (called from .start)
          procedure, public:: load_port=>DSUnitLoadPort           !loads the DS unit port with new DS instructions (called by other DS units)
          procedure, public:: flush_port=>DSUnitFlushPort         !flushes the port content into the DS unit queue (called by the current DS unit)
          procedure, public:: get_dsvp=>DSUnitGetDSVP             !returns a pointer to the DSVP the DS unit is part of
          procedure, public:: get_id=>DSUnitGetId                 !returns the DS unit id
          procedure, private:: set_id=>DSUnitSetId                !sets the DS unit id (when the DSVU table is constructed)
        end type ds_unit_t
 !Wrapped reference to a DSVU:
        type, private:: ds_unit_ref_t
         class(ds_unit_t), pointer, private:: unit_ref=>NULL() !non-owning pointer (reference) to a domain-specific unit
        end type ds_unit_ref_t
 !Domain-specific decoder unit:
        type, abstract, extends(ds_unit_t), public:: ds_decoder_t
         contains
          procedure(ds_decoder_decode_i), deferred, public:: decode !decoding procedure: Unpacks the raw byte packet (instruction bytecode) and constructs a domain-specific instruction
        end type ds_decoder_t
 !Domain-specific encoder unit:
        type, abstract, extends(ds_unit_t), public:: ds_encoder_t
         contains
          procedure(ds_encoder_encode_i), deferred, public:: encode !encoding procedure: Packs a domain-specific instruction into the raw byte packet (instruction bytecode)
        end type ds_encoder_t
 !Domain-specific virtual processor (DSVP):
        type, abstract, public:: dsvp_t
         integer(INTD), private:: stat=DSVP_STAT_OFF       !current DSVP status: {DSVP_STAT_OFF, DSVP_STAT_ON, negative integers = errors}
         integer(INTD), private:: spec_kind=DSVP_NO_KIND   !specific kind of DSVP (to be specilized): Non-negative
         integer(INTL), private:: id=-1                    !DSVP unique ID: [0..max]
         integer(INTL), private:: instr_received=0_INTL    !total number of received domain-specific instructions
         integer(INTL), private:: instr_processed=0_INTL   !total number of processed (retired) instructions (both successful and failed)
         integer(INTL), private:: instr_failed=0_INTL      !total number of retired failed instructions
         type(ds_microcode_t), pointer, private:: microcode(:) !pointer to DS microcode bindings for each DS instruction code: [0..ISA_size-1]: Set by .configure()
         character(:), allocatable, private:: description      !symbolic description of the DSVP: Set by .configure()
         type(ds_unit_ref_t), allocatable, private:: units(:)  !DSVU table (enumerated references to DSVU the DSVP is composed of): Set by .configure()
         integer(INTD), private:: num_units=0                  !number of set DSVU in the DSVU table: [0..num_units-1]: Set by .configure()
         real(8), private:: time_start                         !start time stamp (sec): Set by .start()
         contains
          procedure(dsvp_ctor_i), deferred, public:: configure                   !configures DSVP: Allocates/configures DS units, allocates DSVU table, sets up microcode, sets description, etc.
          procedure, public:: start=>DSVPStart                                   !launches configured DSVP to its life cycle
          procedure, public:: shutdown=>DSVPShutdown                             !shuts down DSVP but keeps it configured
          procedure, public:: destroy=>DSVPDestroy                               !destroys DSVP completely
          procedure, public:: alloc_units=>DSVPAllocUnits                        !allocates the DSVU table
          procedure, public:: free_units=>DSVPFreeUnits                          !frees the DSVU table
          procedure, public:: set_unit=>DSVPSetUnit                              !sets up a DSVU entry in the DSVU table
          procedure, public:: get_unit=>DSVPGetUnit                              !returns a pointer to the DSVU from the DSVU table
          procedure, public:: set_description=>DSVPSetDescription                !sets DSVP ID, kind, and symbolic description
          procedure, public:: get_description=>DSVPGetDescription                !gets DSVP ID, kind, and symbolic description
          procedure, public:: incr_recv_instr_counter=>DSVPIncrRecvInstrCounter  !increments the receieved instruction counter
          procedure, public:: incr_rtrd_instr_counter=>DSVPIncrRtrdInstrCounter  !increments the processed (retired) instruction counter
          procedure, public:: incr_fail_instr_counter=>DSVPIncrFailInstrCounter  !increments the failed instruction counter
          procedure, public:: time_active=>DSVPTimeActive                        !returns the time DSVP is active in seconds
          procedure, public:: get_status=>DSVPGetStatus                          !returns the current status of the DSVP
          procedure, private:: set_status=>DSVPSetStatus                         !sets the DSVP status
          procedure, private:: start_time=>DSVPStartTime                         !starts the time when DSVP is initialized (status set to DSVP_STAT_ON)
        end type dsvp_t
!INTERFACES:
 !Abstract:
        abstract interface
  !ds_resrc_t:
         function ds_resrc_query_i(this,ierr) result(ans)
          import:: ds_resrc_t,INTD
          implicit none
          logical:: ans                               !out: answer
          class(ds_resrc_t), intent(in):: this        !in: domain-specific resource
          integer(INTD), intent(out), optional:: ierr !out: error code
         end function ds_resrc_query_i
  !ds_oprnd_t:
   !self:
         subroutine ds_oprnd_self_i(this,ierr)
          import:: ds_oprnd_t,INTD
          implicit none
          class(ds_oprnd_t), intent(inout):: this     !inout: domain-specific operand
          integer(INTD), intent(out), optional:: ierr !out: error code
         end subroutine ds_oprnd_self_i
   !query:
         function ds_oprnd_query_i(this,ierr) result(res)
          import:: ds_oprnd_t,INTD
          implicit none
          logical:: res                               !out: result
          class(ds_oprnd_t), intent(in):: this        !in: domain-specific operand
          integer(INTD), intent(out), optional:: ierr !out: error code
         end function ds_oprnd_query_i
   !sync:
         function ds_oprnd_sync_i(this,ierr,wait) result(res)
          import:: ds_oprnd_t,INTD
          implicit none
          logical:: res                               !out: result
          class(ds_oprnd_t), intent(inout):: this     !in: domain-specific operand
          integer(INTD), intent(out), optional:: ierr !out: error code
          logical, intent(in), optional:: wait        !in: TRUE activates WAIT instead of TEST synchronization
         end function ds_oprnd_sync_i
  !ds_instr_ctrl_t:
   !pack:
         subroutine ds_instr_ctrl_pack_i(this,packet,ierr)
          import:: ds_instr_ctrl_t,obj_pack_t,INTD
          class(ds_instr_ctrl_t), intent(in):: this   !in: domain-specific instruction control
          class(obj_pack_t), intent(inout):: packet   !out: packet
          integer(INTD), intent(out), optional:: ierr !out: error code
         end subroutine ds_instr_ctrl_pack_i
   !unpack:
         subroutine ds_instr_ctrl_unpack_i(this,packet,ierr)
          import:: ds_instr_ctrl_t,obj_pack_t,INTD
          class(ds_instr_ctrl_t), intent(inout):: this !out: domain-specific instruction control
          class(obj_pack_t), intent(inout):: packet    !in: packet
          integer(INTD), intent(out), optional:: ierr  !out: error code
         end subroutine ds_instr_ctrl_unpack_i
  !ds_instr_t:
   !self:
         subroutine ds_instr_self_i(this,ierr)
          import:: ds_instr_t,INTD
          implicit none
          class(ds_instr_t), intent(inout):: this     !inout: domain-specific instruction
          integer(INTD), intent(out), optional:: ierr !out: error code
         end subroutine ds_instr_self_i
   !sync:
         function ds_instr_sync_i(this,ierr,wait) result(res)
          import:: ds_instr_t,INTD
          implicit none
          logical:: res                               !out: TRUE if synchronized
          class(ds_instr_t), intent(inout):: this     !inout: domain-specific instruction
          integer(INTD), intent(out), optional:: ierr !out: error code
          logical, intent(in), optional:: wait        !in: WAIT versus TEST (defaults to WAIT)
         end function ds_instr_sync_i
   !encode:
         subroutine ds_instr_encode_i(this,instr_packet,ierr)
          import:: ds_instr_t,obj_pack_t,INTD
          class(ds_instr_t), intent(in):: this            !in: domain-specific instruction to be encoded
          class(obj_pack_t), intent(inout):: instr_packet !out: instruction byte packet (bytecode)
          integer(INTD), intent(out), optional:: ierr     !out: error code
         end subroutine ds_instr_encode_i
  !ds_unit_t:
   !ctor:
         subroutine ds_unit_ctor_i(this,conf,ierr)
          import:: ds_unit_t,dsv_conf_t,INTD
          implicit none
          class(ds_unit_t), intent(inout):: this      !out: configured (constructed) DSVP
          class(dsv_conf_t), intent(in):: conf        !in: DSVU configuration
          integer(INTD), intent(out), optional:: ierr !out: error code
         end subroutine ds_unit_ctor_i
   !self:
         subroutine ds_unit_self_i(this,ierr)
          import:: ds_unit_t,INTD
          class(ds_unit_t), intent(inout):: this      !inout: DS unit
          integer(INTD), intent(out), optional:: ierr !out: error code
         end subroutine ds_unit_self_i
  !ds_decoder_t:
         subroutine ds_decoder_decode_i(this,ds_instr,instr_packet,ierr)
          import:: ds_decoder_t,ds_instr_t,obj_pack_t,INTD
          class(ds_decoder_t), intent(inout):: this               !inout: DS decoder unit
          class(ds_instr_t), intent(inout), target:: ds_instr     !out: decoded domain-specific instruction ready for the DS pipeline
          class(obj_pack_t), intent(inout):: instr_packet         !in: instruction byte packet (bytecode)
          integer(INTD), intent(out), optional:: ierr             !out: error code
         end subroutine ds_decoder_decode_i
  !ds_encoder_t:
         subroutine ds_encoder_encode_i(this,ds_instr,instr_packet,ierr)
          import:: ds_encoder_t,ds_instr_t,obj_pack_t,INTD
          class(ds_encoder_t), intent(inout):: this               !inout: DS encoder unit
          class(ds_instr_t), intent(in), target:: ds_instr        !in: domain-specific instruction
          class(obj_pack_t), intent(inout):: instr_packet         !out: instruction byte packet (bytecode)
          integer(INTD), intent(out), optional:: ierr             !out: error code
         end subroutine ds_encoder_encode_i
  !dsvp_t:
   !ctor:
         subroutine dsvp_ctor_i(this,conf,ierr)
          import:: dsvp_t,dsv_conf_t,INTD
          implicit none
          class(dsvp_t), intent(inout):: this         !out: configured (constructed) DSVP
          class(dsv_conf_t), intent(in):: conf        !in: DSVP configuration
          integer(INTD), intent(out), optional:: ierr !out: error code
         end subroutine dsvp_ctor_i
        end interface
!VISIBILITY:
 !ds_resrc_t:
        public ds_resrc_query_i
 !ds_oprnd_t:
        private DSOprndIsActive
        private DSOprndIsPresent
        private DSOprndMarkActive
        private DSOprndMarkEmpty
        private DSOprndMarkDelivered
        private DSOprndMarkUndelivered
        private DsOprndGetCommStat
        private DSOprndSetCommStat
        public ds_oprnd_self_i
        public ds_oprnd_query_i
        public ds_oprnd_sync_i
 !ds_instr_ctrl_t:
        public ds_instr_ctrl_pack_i
        public ds_instr_ctrl_unpack_i
 !ds_instr_t:
        private DSInstrIsEmpty
        private DSInstrIsRetired
        private DSInstrIsActive
        private DSInstrGetCode
        private DSInstrSetCode
        private DSInstrGetStatus
        private DSInstrSetStatus
        private DSInstrGetControl
        private DSInstrSetControl
        private DSInstrFreeControl
        private DSInstrAllocOperands
        private DSInstrDeallocOperands
        private DSInstrGetOperand
        private DSInstrSetOperand
        private DSInstrFreeOperand
        private DSInstrGetNumOperands
        private DSInstrAllSet
        private DSInstrActivate
        private DSInstrTerminate
        private DSInstrClean
        public ds_instr_self_i
        public ds_instr_sync_i
        public ds_instr_encode_i
 !ds_microcode_t:
        private DSMicrocodeReset
 !ds_unit_port_t:
        private DSUnitPortAccept
        private DSUnitPortAbsorb
        private DSUnitPortFree
        public ds_unit_port_dtor
 !ds_unit_t:
        private DSUnitLoadPort
        private DSUnitFlushPort
        private DSUnitGetDSVP
        private DSUnitGetId
        private DSUnitSetId
        public ds_unit_ctor_i
        public ds_unit_self_i
 !ds_decoder_t:
        public ds_decoder_decode_i
 !ds_encoder_t:
        public ds_encoder_encode_i
 !dsvp_t:
        private DSVPStart
        private DSVPShutdown
        private DSVPDestroy
        private DSVPAllocUnits
        private DSVPFreeUnits
        private DSVPSetUnit
        private DSVPGetUnit
        private DSVPSetDescription
        private DSVPGetDescription
        private DSVPIncrRecvInstrCounter
        private DSVPIncrRtrdInstrCounter
        private DSVPIncrFailInstrCounter
        private DSVPTimeActive
        private DSVPGetStatus
        private DSVPSetStatus
        private DSVPStartTime
        public dsvp_ctor_i
!IMPLEMENTATION:
       contains
![ds_oprnd_t]==========================================
        function DSOprndIsActive(this,ierr) result(res)
!Returns TRUE if the domain-specific operand is active (defined).
         implicit none
         logical:: res                               !out: answer
         class(ds_oprnd_t), intent(in):: this        !in: domain-specific operand
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc

         errc=DSVP_SUCCESS
         res=(this%stat.gt.DS_OPRND_EMPTY)
         if(present(ierr)) ierr=errc
         return
        end function DSOprndIsActive
!-------------------------------------------------------
        function DSOprndIsPresent(this,ierr) result(res)
!Returns TRUE if the domain-specific operand is locally available (present).
!The input domain-specific operand must have already been defined.
         implicit none
         logical:: res                               !out: answer
         class(ds_oprnd_t), intent(in):: this        !in: defined domain-specific operand
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc

         errc=DSVP_SUCCESS; res=.FALSE.
         if(this%is_active(errc)) then
          if(errc.eq.DSVP_SUCCESS) res=(this%stat.gt.DS_OPRND_DEFINED)
         else
          errc=DSVP_ERR_INVALID_REQ
         endif
         if(present(ierr)) ierr=errc
         return
        end function DSOprndIsPresent
!----------------------------------------------
        subroutine DSOprndMarkActive(this,ierr)
!Marks the domain-specific operand active (defined).
!Trying to mark an already active operand active again will result in an error.
         implicit none
         class(ds_oprnd_t), intent(inout):: this     !inout: domain-specific operand
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc

         errc=DSVP_SUCCESS
         if(.not.this%is_active(errc)) then
          if(errc.eq.DSVP_SUCCESS) this%stat=DS_OPRND_DEFINED
         else
          errc=DSVP_ERR_INVALID_REQ
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine DSOprndMarkActive
!---------------------------------------------
        subroutine DSOprndMarkEmpty(this,ierr)
!Marks the domain-specific operand as empty (undefined).
!The local resources will automatically be released.
!It is allowed to call this procedure on an empty operand.
         implicit none
         class(ds_oprnd_t), intent(inout):: this     !inout: domain-specific operand
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc,ier

         errc=DSVP_SUCCESS
         if(this%is_active(errc)) then
          if(errc.eq.DSVP_SUCCESS) then
           if(this%in_route.eq.DS_OPRND_NO_COMM) then
            if(this%is_present(ier)) call this%mark_undelivered(errc) !will release local resources
            if(errc.eq.DSVP_SUCCESS) errc=ier
            this%stat=DS_OPRND_EMPTY
           else
            errc=DSVP_ERR_INVALID_REQ
           endif
          endif
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine DSOprndMarkEmpty
!------------------------------------------------
        subroutine DSOprndMarkDelivered(this,ierr)
!Marks the domain-specific operand as delivered (present, locally available).
!Trying to mark an already delivered operand delivered again will cause an error.
         implicit none
         class(ds_oprnd_t), intent(inout):: this     !inout: domain-specific operand
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc

         errc=DSVP_SUCCESS
         if(this%is_active(errc)) then
          if(errc.eq.DSVP_SUCCESS) then
           if(.not.this%is_present(errc)) then
            call this%set_comm_stat(DS_OPRND_NO_COMM,errc)
            if(errc.eq.DSVP_SUCCESS) this%stat=DS_OPRND_PRESENT
           else
            errc=DSVP_ERR_INVALID_REQ
           endif
          endif
         else
          errc=DSVP_ERR_INVALID_REQ
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine DSOprndMarkDelivered
!-----------------------------------------------------------
        subroutine DSOprndMarkUndelivered(this,ierr,sync_it)
!Marks the domain-specific operand as undelivered (locally unavailable).
!The local resources are automatically released, but the operand stays defined.
!It is allowed to call this procedure on an undelivered operand. However, trying to
!mark undelivered an operand with a pending communication will cause an error,
!unless the optional parameter <sync_it> is set to TRUE (will enforce synchronization).
         implicit none
         class(ds_oprnd_t), intent(inout):: this     !inout: domain-specific operand
         integer(INTD), intent(out), optional:: ierr !out: error code
         logical, intent(in), optional:: sync_it     !in: if TRUE, a possible pending communication will be completed before resource release
         integer(INTD):: errc,ier
         logical:: compl_comm,dirty

         errc=DSVP_SUCCESS; dirty=.FALSE.
         if(this%is_active(errc)) then
          if(errc.eq.DSVP_SUCCESS) then
           compl_comm=.FALSE.; if(present(sync_it)) compl_comm=sync_it
           if(compl_comm) then
            compl_comm=this%sync(errc,wait=.TRUE.); dirty=(errc.ne.DSVP_SUCCESS)
            errc=DSVP_SUCCESS
           endif
           if((this%in_route.eq.DS_OPRND_NO_COMM).or.dirty) then !no pending communication check
            if(this%is_present(ier)) call this%release(errc) !will release local resources
            if(errc.eq.DSVP_SUCCESS) errc=ier
            this%stat=DS_OPRND_DEFINED !status will be changed regardless the success of resource release
           else
            errc=DSVP_ERR_INVALID_REQ
           endif
          endif
         else
          if(errc.eq.DSVP_SUCCESS) errc=DSVP_ERR_INVALID_REQ
         endif
         if(dirty) errc=NOT_CLEAN
         if(present(ierr)) ierr=errc
         return
        end subroutine DSOprndMarkUndelivered
!---------------------------------------------------------
        function DSOprndGetCommStat(this,ierr) result(sts)
!Gets the current communication status on the domain-specific operand.
         implicit none
         integer(INTD):: sts                         !out: current communication status
         class(ds_oprnd_t), intent(in):: this        !in: domain-specific operand
         integer(INTD), intent(out), optional:: ierr !out: error code

         sts=this%stat
         if(present(ierr)) ierr=DSVP_SUCCESS
         return
        end function DSOprndGetCommStat
!---------------------------------------------------
        subroutine DSOprndSetCommStat(this,sts,ierr)
!Sets the communication status on the domain-specific operand.
         implicit none
         class(ds_oprnd_t), intent(inout):: this     !inout: domain-specific operand
         integer(INTD), intent(in):: sts             !in: communication status to be set
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc

         if(this%is_active(errc)) then
          if(errc.eq.DSVP_SUCCESS) then
           select case(sts)
           case(DS_OPRND_NO_COMM)
            this%in_route=sts
           case(DS_OPRND_FETCHING,DS_OPRND_UPLOADING)
            if(this%in_route.eq.DS_OPRND_NO_COMM) then
             this%in_route=sts
            else
             errc=DSVP_ERR_INVALID_REQ
            endif
           case default
            errc=DSVP_ERR_INVALID_ARGS
           end select
          endif
         else
          errc=DSVP_ERR_INVALID_REQ
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine DSOprndSetCommStat
![ds_instr_t]=========================================
        function DSInstrIsEmpty(this,ierr) result(res)
!Returns TRUE if the domain-specific instruction is empty (undefined).
         implicit none
         logical:: res                               !out: answer
         class(ds_instr_t), intent(in):: this        !in: domain-specific instruction
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc

         errc=DSVP_SUCCESS; res=.TRUE.
         if(this%stat.ne.DS_INSTR_EMPTY) then
          if(this%code.ne.DS_INSTR_NOOP) then
           res=.FALSE.
          else
           errc=DSVP_ERR_BROKEN_OBJ
          endif
         endif
         if(present(ierr)) ierr=errc
         return
        end function DSInstrIsEmpty
!-------------------------------------------------------
        function DSInstrIsRetired(this,ierr) result(res)
!Returns TRUE if the domain-specific instruction is retired.
         implicit none
         logical:: res                               !out: answer
         class(ds_instr_t), intent(in):: this        !in: domain-specific instruction
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc

         errc=DSVP_SUCCESS; res=.FALSE.
         if(this%stat.eq.DS_INSTR_RETIRED) then
          if(this%code.ne.DS_INSTR_NOOP) then
           res=.TRUE.
          else
           errc=DSVP_ERR_BROKEN_OBJ
          endif
         endif
         if(present(ierr)) ierr=errc
         return
        end function DSInstrIsRetired
!------------------------------------------------------
        function DSInstrIsActive(this,ierr) result(res)
!Returns TRUE if the domain-specific instruction is active (defined).
         implicit none
         logical:: res                               !out: answer
         class(ds_instr_t), intent(in):: this        !in: domain-specific instruction
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc

         res=.not.this%is_empty(errc)
         if(present(ierr)) ierr=errc
         return
        end function DSInstrIsActive
!------------------------------------------------------
        function DSInstrGetCode(this,ierr) result(code)
!Returns the instruction code. All valid codes are non-negative.
!A negative code, e.g., DS_INSTR_NOOP=-1, means an empty instruction.
         implicit none
         integer(INTD):: code                        !out: instruction code
         class(ds_instr_t), intent(in):: this        !in: domain-specific instruction
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc

         errc=DSVP_SUCCESS
         code=this%code
         if(present(ierr)) ierr=errc
         return
        end function DSInstrGetCode
!------------------------------------------------
        subroutine DSInstrSetCode(this,code,ierr)
!Sets the instruction code. All valid codes are non-negative.
!A negative code, e.g., DS_INSTR_NOOP=-1, means an empty instruction.
         implicit none
         class(ds_instr_t), intent(inout):: this     !inout: domain-specific instruction
         integer(INTD), intent(in):: code            !in: instruction code
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc

         errc=DSVP_SUCCESS
         this%code=code
         if(present(ierr)) ierr=errc
         return
        end subroutine DSInstrSetCode
!-------------------------------------------------------
        function DSInstrGetStatus(this,ierr) result(sts)
!Returns the instruction status.
         implicit none
         integer(INTD):: sts                         !out: instruction status
         class(ds_instr_t), intent(in):: this        !in: domain-specific instruction
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc

         errc=DSVP_SUCCESS
         sts=this%stat
         if(present(ierr)) ierr=errc
         return
        end function DSInstrGetStatus
!------------------------------------------------------------
        subroutine DSInstrSetStatus(this,sts,ierr,error_code)
!Sets the instruction status and (optionally) error code.
         implicit none
         class(ds_instr_t), intent(inout):: this          !inout: domain-specific instruction
         integer(INTD), intent(in):: sts                  !in: instruction status to be set
         integer(INTD), intent(out), optional:: ierr      !out: error code
         integer(INTD), intent(in), optional:: error_code !in: instruction error code to set
         integer(INTD):: errc

         errc=DSVP_SUCCESS
         this%stat=sts
         if(present(error_code)) this%error_code=error_code
         if(present(ierr)) ierr=errc
         return
        end subroutine DSInstrSetStatus
!---------------------------------------------------------
        function DSInstrGetControl(this,ierr) result(ctrl)
!Returns a pointer to the instruction control field.
         implicit none
         class(ds_instr_ctrl_t), pointer:: ctrl      !out: pointer to the instruction control field
         class(ds_instr_t), intent(in):: this        !in: domain-specific instruction
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc

         errc=DSVP_SUCCESS
         if(associated(this%control)) then
          ctrl=>this%control
         else
          errc=DSVP_ERR_INVALID_REQ
         endif
         if(present(ierr)) ierr=errc
         return
        end function DSInstrGetControl
!---------------------------------------------------
        subroutine DSInstrSetControl(this,ctrl,ierr)
!Associates the instruction control field.
         implicit none
         class(ds_instr_t), intent(inout):: this            !inout: domain-specific instruction
         class(ds_instr_ctrl_t), pointer, intent(in):: ctrl !in: pointer to a control field (target)
         integer(INTD), intent(out), optional:: ierr        !out: error code
         integer(INTD):: errc

         errc=DSVP_SUCCESS
         if(associated(ctrl)) then
          if(.not.associated(this%control)) then
           this%control=>ctrl
          else
           errc=DSVP_ERR_INVALID_REQ
          endif
         else
          errc=DSVP_ERR_INVALID_ARGS
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine DSInstrSetControl
!-----------------------------------------------------------
        subroutine DSInstrFreeControl(this,ierr,dissoc_only)
!Frees the instructon control field. By default it will be deallocated,
!but if <dissoc_only>=TRUE, only dissociation will happen.
         implicit none
         class(ds_instr_t), intent(inout):: this     !inout: domain-specific instruction
         integer(INTD), intent(out), optional:: ierr !out: error code
         logical, intent(in), optional:: dissoc_only !in: if TRUE, only dissociation, no deallocation
         integer(INTD):: errc
         integer:: ier

         errc=DSVP_SUCCESS
         if(associated(this%control)) then
          ier=1; if(present(dissoc_only)) then; if(dissoc_only) ier=0; endif
          if(ier.ne.0) then
           deallocate(this%control,STAT=ier); if(ier.ne.0) errc=DSVP_ERR_MEM_FREE_FAIL
          endif
          this%control=>NULL()
         else
          errc=DSVP_ERR_INVALID_REQ
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine DSInstrFreeControl
!--------------------------------------------------------------
        subroutine DSInstrAllocOperands(this,num_operands,ierr)
!Allocates the operands reference storage (to an empty state).
         implicit none
         class(ds_instr_t), intent(inout):: this     !inout: domain-specific instruction
         integer(INTD), intent(in):: num_operands    !in: number of operands
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc
         integer:: ier

         errc=DSVP_SUCCESS
         if(num_operands.gt.0) then
          if(.not.allocated(this%operand)) then
           allocate(this%operand(0:num_operands-1),STAT=ier)
           if(ier.ne.0) errc=DSVP_ERR_MEM_ALLOC_FAIL
           if(errc.eq.DSVP_SUCCESS) then
            if(this%num_oprnds.ne.0) errc=DSVP_ERR_BROKEN_OBJ
            this%num_oprnds=num_operands
           endif
          else
           errc=DSVP_ERR_INVALID_REQ
          endif
         else
          errc=DSVP_ERR_INVALID_ARGS
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine DSInstrAllocOperands
!---------------------------------------------------------------
        subroutine DSInstrDeallocOperands(this,ierr,dissoc_only)
!Deallocates the operands reference storage. By default, a deallocation
!will be performed on each operand, unless <dissoc_only>=TRUE.
         implicit none
         class(ds_instr_t), intent(inout):: this     !inout: domain-specific instruction
         integer(INTD), intent(out), optional:: ierr !out: error code
         logical, intent(in), optional:: dissoc_only !in: if TRUE, only dissociation will be performed on operands
         integer(INTD):: errc,ier,i
         integer:: j
         logical:: dis

         errc=DSVP_SUCCESS
         if(allocated(this%operand)) then
          if(this%num_oprnds.gt.0) then
           if(present(dissoc_only)) then; dis=dissoc_only; else; dis=.FALSE.; endif
           do i=this%num_oprnds-1,0,-1
            call this%free_operand(i,ier,dis); if(ier.ne.DSVP_SUCCESS) errc=ier
           enddo
           if(errc.eq.DSVP_SUCCESS) then
            deallocate(this%operand,STAT=j); if(j.ne.0) errc=DSVP_ERR_MEM_FREE_FAIL
            this%num_oprnds=0
           endif
          else
           errc=DSVP_ERR_BROKEN_OBJ
          endif
         else
          errc=DSVP_ERR_INVALID_REQ
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine DSInstrDeallocOperands
!-------------------------------------------------------------------
        function DSInstrGetOperand(this,op_num,ierr) result(oprnd_p)
!Provides an access to a specific operand of the domain-specific instruction.
         implicit none
         class(ds_oprnd_t), pointer:: oprnd_p         !out: pointer to a specific operand
         class(ds_instr_t), intent(in):: this         !in: domain-specific instruction
         integer(INTD), intent(in):: op_num           !in: operand number (0,1,2,...)
         integer(INTD), intent(out), optional:: ierr  !out: error code
         integer(INTD):: errc

         errc=DSVP_SUCCESS; oprnd_p=>NULL()
         if(this%is_active(errc)) then
          if(errc.eq.DSVP_SUCCESS) then
           if(this%num_oprnds.gt.0) then
            if(op_num.ge.0.and.op_num.lt.this%num_oprnds) then
#ifdef DSVP_DEBUG
             if(DEBUG.gt.0) then
              if(allocated(this%operand)) then
               if(lbound(this%operand,1).ne.0.or.ubound(this%operand,1)+1.ne.this%num_oprnds) errc=DSVP_ERR_BROKEN_OBJ
              else
               errc=DSVP_ERR_BROKEN_OBJ
              endif
             endif
             if(errc.eq.DSVP_SUCCESS) then
#endif
              oprnd_p=>this%operand(op_num)%oprnd_ref
#ifdef DSVP_DEBUG
             endif
#endif
            else
             errc=DSVP_ERR_INVALID_ARGS
            endif
           else
            errc=DSVP_ERR_INVALID_REQ
           endif
          endif
         else
          errc=DSVP_ERR_INVALID_REQ
         endif
         if(present(ierr)) ierr=errc
         return
        end function DSInstrGetOperand
!-----------------------------------------------------------
        subroutine DSInstrSetOperand(this,op_num,oprnd,ierr)
!Associates a specific operand of the domain-specific instruction with its target.
         implicit none
         class(ds_instr_t), intent(inout):: this        !inout: domain-specific instruction
         integer(INTD), intent(in):: op_num             !in: operand number (0,1,2,...)
         class(ds_oprnd_t), pointer, intent(in):: oprnd !in: domain-specific operand (target)
         integer(INTD), intent(out), optional:: ierr    !out: error code
         integer(INTD):: errc

         errc=DSVP_SUCCESS
         if(associated(oprnd)) then
          if(this%num_oprnds.gt.0) then
           if(op_num.ge.0.and.op_num.lt.this%num_oprnds) then
#ifdef DSVP_DEBUG
            if(DEBUG.gt.0) then
             if(allocated(this%operand)) then
              if(lbound(this%operand,1).ne.0.or.ubound(this%operand,1)+1.ne.this%num_oprnds) errc=DSVP_ERR_BROKEN_OBJ
             else
              errc=DSVP_ERR_BROKEN_OBJ
             endif
            endif
            if(errc.eq.DSVP_SUCCESS) then
#endif
             if(.not.associated(this%operand(op_num)%oprnd_ref)) then
              this%operand(op_num)%oprnd_ref=>oprnd
             else
              errc=DSVP_ERR_INVALID_REQ
             endif
#ifdef DSVP_DEBUG
            endif
#endif
           else
            errc=DSVP_ERR_INVALID_ARGS
           endif
          else
           errc=DSVP_ERR_INVALID_REQ
          endif
         else
          errc=DSVP_ERR_INVALID_ARGS
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine DSInstrSetOperand
!------------------------------------------------------------------
        subroutine DSInstrFreeOperand(this,op_num,ierr,dissoc_only)
!Frees a specific instruction operand. By default the operand poiter will be
!deallocated, unless <dissoc_only>=TRUE, which will only dissociate it.
         implicit none
         class(ds_instr_t), intent(inout):: this     !inout: domain-specific instruction
         integer(INTD), intent(in):: op_num          !in: operand number (0,1,3,...)
         integer(INTD), intent(out), optional:: ierr !out: error code
         logical, intent(in), optional:: dissoc_only !in: if TRUE, no pointer deallocation will be performed
         integer(INTD):: errc
         integer:: ier
         logical:: dis

         errc=DSVP_SUCCESS
         if(present(dissoc_only)) then; dis=dissoc_only; else; dis=.FALSE.; endif
         if(this%num_oprnds.gt.0) then
          if(op_num.ge.0.and.op_num.lt.this%num_oprnds) then
           call this%operand(op_num)%oprnd_ref%mark_empty(errc) !will call destructor (release local resources)
           if(.not.dis) then
            deallocate(this%operand(op_num)%oprnd_ref,STAT=ier)
            if(ier.ne.0.and.errc.eq.DSVP_SUCCESS) errc=DSVP_ERR_MEM_FREE_FAIL
           endif
           this%operand(op_num)%oprnd_ref=>NULL()
          else
           errc=DSVP_ERR_INVALID_ARGS
          endif
         else
          errc=DSVP_ERR_INVALID_REQ
          if(allocated(this%operand)) errc=DSVP_ERR_BROKEN_OBJ
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine DSInstrFreeOperand
!---------------------------------------------------------------------
        function DSInstrGetNumOperands(this,ierr) result(num_operands)
!Returns the number of instruction operands.
         implicit none
         integer(INTD):: num_operands                !out: number of operands
         class(ds_instr_t), intent(in):: this        !in: domain-specific instruction
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc

         errc=DSVP_SUCCESS
         num_operands=this%num_oprnds
         if(num_operands.lt.0) errc=DSVP_ERR_BROKEN_OBJ
         if(present(ierr)) ierr=errc
         return
        end function DSInstrGetNumOperands
!----------------------------------------------------
        function DSInstrAllSet(this,ierr) result(res)
!Returns TRUE if the domain-specific instruction is fully defined.
         implicit none
         logical:: res                               !out: answer
         class(ds_instr_t), intent(in):: this        !in: domain-specific instruction
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc,i

         errc=DSVP_SUCCESS; res=.FALSE.
         if(this%code.eq.DS_INSTR_NOOP.or.this%stat.eq.DS_INSTR_EMPTY) return
         if(.not.associated(this%control)) return
         if(this%num_oprnds.gt.0) then
          if(allocated(this%operand)) then
           do i=0,this%num_oprnds-1
            if(.not.associated(this%operand(i)%oprnd_ref)) return
           enddo
          else
           errc=DSVP_ERR_BROKEN_OBJ
          endif
         endif
         if(.not.associated(this%acquire_resource)) return
         if(.not.associated(this%prefetch_input)) return
         if(.not.associated(this%sync_prefetch)) return
         if(.not.associated(this%execute)) return
         if(.not.associated(this%sync_execution)) return
         if(.not.associated(this%upload_output)) return
         if(.not.associated(this%sync_upload)) return
         if(.not.associated(this%release_resource)) return
         if(errc.eq.DSVP_SUCCESS) res=.TRUE.
         if(present(ierr)) ierr=errc
         return
        end function DSInstrAllSet
!---------------------------------------------------------------
        subroutine DSInstrActivate(this,op_code,micro_bind,ierr)
!Activates the domain-specific instruction after it has been constructed or decoded.
         implicit none
         class(ds_instr_t), intent(inout):: this        !inout: domain-specific instruction
         integer(INTD), intent(in):: op_code            !in: instruction code
         class(ds_microcode_t), intent(in):: micro_bind !in: concrete microcode binding for the domain-specific instruction
         integer(INTD), intent(out), optional:: ierr    !out: error code
         integer(INTD):: errc

         call this%set_code(op_code,errc)
         if(errc.eq.DSVP_SUCCESS) then
          if(this%get_status().eq.DS_INSTR_EMPTY) then
           call set_microcode(errc)
           if(errc.eq.0) then
            call this%set_status(DS_INSTR_NEW,errc,DSVP_SUCCESS)
            if(errc.ne.DSVP_SUCCESS) then
             call this%set_status(DS_INSTR_RETIRED,errc,DSVP_ERROR)
             errc=-4
            endif
           else
            call this%set_status(DS_INSTR_RETIRED,errc,DSVP_ERR_INVALID_ARGS)
            errc=-3
           endif
          else
           call this%set_status(DS_INSTR_RETIRED,errc,DSVP_ERR_INVALID_REQ)
           errc=-2
          endif
         else
          call this%set_status(DS_INSTR_RETIRED,errc,DSVP_ERR_UNABLE_COMPLETE)
          errc=-1
         endif
         if(present(ierr)) ierr=errc
         return

        contains

         subroutine set_microcode(jerr)
          implicit none
          integer(INTD), intent(out):: jerr

          jerr=0
          this%acquire_resource=>micro_bind%acquire_resource
          this%prefetch_input=>micro_bind%prefetch_input
          this%sync_prefetch=>micro_bind%sync_prefetch
          this%execute=>micro_bind%execute
          this%sync_execution=>micro_bind%sync_execution
          this%upload_output=>micro_bind%upload_output
          this%sync_upload=>micro_bind%sync_upload
          this%release_resource=>micro_bind%release_resource
          return
         end subroutine set_microcode

        end subroutine DSInstrActivate
!--------------------------------------------------------
        subroutine DSInstrTerminate(this,error_code,ierr)
!Terminates the normal instruction execution workflow,
!but leaves the instruction defined (active).
!The instruction status will set to DS_INSTR_RETIRED with
!a given <error_code>.
         implicit none
         class(ds_instr_t), intent(inout):: this     !inout: domain-specific instruction
         integer(INTD), intent(in):: error_code      !in: instruction error code to set
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc,ier,sts
         logical:: nostat,res

         if(this%is_active(errc)) then
          sts=this%get_status(errc)
          if(errc.eq.DSVP_SUCCESS) then
           if(sts.ge.DS_INSTR_INPUT_WAIT) then
            if(sts.ge.DS_INSTR_READY_TO_EXEC) then
             if(sts.ge.DS_INSTR_OUTPUT_WAIT) then
              if(sts.lt.DS_INSTR_RETIRED) then
               res=this%sync_upload(ier); if(ier.ne.DSVP_SUCCESS.and.errc.eq.DSVP_SUCCESS) errc=ier
              endif
             else
              if(sts.ge.DS_INSTR_SCHEDULED) then
               res=this%sync_execution(ier); if(ier.ne.DSVP_SUCCESS.and.errc.eq.DSVP_SUCCESS) errc=ier
              endif
             endif
            else
             res=this%sync_prefetch(ier); if(ier.ne.DSVP_SUCCESS.and.errc.eq.DSVP_SUCCESS) errc=ier
            endif
           endif
           if(sts.ge.DS_INSTR_RSC_WAIT) then
            call this%release_resource(ier); if(ier.ne.DSVP_SUCCESS.and.errc.eq.DSVP_SUCCESS) errc=ier
           endif
          else !instruction status unknown
           res=this%sync_prefetch(ier)
           res=this%sync_execution(ier)
           res=this%sync_upload(ier)
           call this%release_resource(ier)
          endif
          call this%set_status(DS_INSTR_RETIRED,ier,error_code)
          if(ier.ne.DSVP_SUCCESS.and.errc.eq.DSVP_SUCCESS) errc=ier
         else
          if(errc.eq.DSVP_SUCCESS) errc=DSVP_ERR_INVALID_REQ
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine DSInstrTerminate
!-----------------------------------------------------
        subroutine DSInstrClean(this,ierr,dissoc_only)
!Resets the domain-specific instruction to an empty state. By default,
!the instruction control field and instruction operands will be deallocated,
!unless <dissoc_only>=TRUE, which will only dissociate the corresponing pointers.
         implicit none
         class(ds_instr_t), intent(inout):: this     !inout: domain-specific instruction
         integer(INTD), intent(out), optional:: ierr !out: error code
         logical, intent(in), optional:: dissoc_only !in: if TRUE, the control and individual operands will simply be dissociated
         integer(INTD):: errc,ier
         logical:: dis

         errc=DSVP_SUCCESS; if(present(dissoc_only)) then; dis=dissoc_only; else; dis=.FALSE.; endif
         call this%dealloc_operands(ier,dis); if(ier.ne.DSVP_SUCCESS.and.errc.eq.DSVP_SUCCESS) errc=ier
         call this%free_control(ier,dis); if(ier.ne.DSVP_SUCCESS.and.errc.eq.DSVP_SUCCESS) errc=ier
         this%code=DS_INSTR_NOOP; this%stat=DS_INSTR_EMPTY
         this%acquire_resource=>NULL()
         this%prefetch_input=>NULL()
         this%sync_prefetch=>NULL()
         this%execute=>NULL()
         this%sync_execution=>NULL()
         this%upload_output=>NULL()
         this%sync_upload=>NULL()
         this%release_resource=>NULL()
         if(present(ierr)) ierr=errc
         return
        end subroutine DSInstrClean
![ds_microcode_t]========================
        subroutine DSMicrocodeReset(this)
         implicit none
         class(ds_microcode_t), intent(inout):: this !inout: microcode binding

         this%acquire_resource=>NULL()
         this%prefetch_input=>NULL()
         this%sync_prefetch=>NULL()
         this%execute=>NULL()
         this%sync_execution=>NULL()
         this%upload_output=>NULL()
         this%sync_upload=>NULL()
         this%release_resource=>NULL()
         return
        end subroutine DSMicrocodeReset
![ds_unit_port_t]============================================
        function DSUnitPortAccept(this,new_list) result(ierr)
!Accepts new DS instructions from other DS units in the current DS unit port.
!The new DS instructions are assumed stored by reference in the <new_list> and
!they will be moved into the port, thus leaving <new_list> empty at the end.
         implicit none
         integer(INTD):: ierr                        !out: error code
         class(ds_unit_port_t), intent(inout):: this !inout: DS unit port
         type(list_bi_t), intent(inout):: new_list   !inout: list of new DS instructions for the DS unit (from other DS units): List items are stored by reference
         type(list_iter_t):: ilist
         integer(INTD):: errc

         ierr=DSVP_SUCCESS
!$OMP CRITICAL (DSVU_PORT_LOCK)
         errc=ilist%init(new_list)
         if(errc.eq.GFC_SUCCESS) then
          errc=ilist%move_list(this%iqueue); if(errc.ne.GFC_SUCCESS) ierr=DSVP_ERROR
         else
          ierr=DSVP_ERROR
         endif
         errc=ilist%release(); if(errc.ne.GFC_SUCCESS.and.ierr.eq.DSVP_SUCCESS) ierr=DSVP_ERROR
!$OMP END CRITICAL (DSVU_PORT_LOCK)
         return
        end function DSUnitPortAccept
!-----------------------------------------------------------------
        function DSUnitPortAbsorb(this,dsvu_queue_it) result(ierr)
!Absorbs new DS instructions from the DS unit port into the DS unit queue.
         implicit none
         integer(INTD):: ierr                             !out: error code
         class(ds_unit_port_t), intent(inout):: this      !inout: DS unit port
         type(list_iter_t), intent(inout):: dsvu_queue_it !inout: DS unit queue iterator
         integer(INTD):: errc

         ierr=DSVP_SUCCESS
!$OMP CRITICAL (DSVU_PORT_LOCK)
         errc=this%iqueue%reset()
         if(errc.eq.GFC_SUCCESS) then
          errc=this%iqueue%move_list(dsvu_queue_it)
          if(errc.ne.GFC_SUCCESS) ierr=DSVP_ERROR
         else
          ierr=DSVP_ERROR
         endif
!$OMP END CRITICAL (DSVU_PORT_LOCK)
         return
        end function DSUnitPortAbsorb
!-------------------------------------------------
        function DSUnitPortFree(this) result(ierr)
!Deletes all DS instructions from the port and resets everything.
         implicit none
         integer(INTD):: ierr                        !out: error code
         class(ds_unit_port_t), intent(inout):: this !inout: DS unit port

!$OMP CRITICAL (DSVU_PORT_LOCK)
         ierr=this%iqueue%delete_all()
!$OMP END CRITICAL (DSVU_PORT_LOCK)
         return
        end function DSUnitPortFree
!-----------------------------------------
        subroutine ds_unit_port_dtor(this)
         implicit none
         type(ds_unit_port_t):: this
         integer(INTD):: ierr

         ierr=this%free()
         return
        end subroutine ds_unit_port_dtor
![ds_unit_t]===============================================
        function DSUnitLoadPort(this,new_list) result(ierr)
!Called by other DS units in order to put new DS instructions
!into the port of the current DS unit. <new_list> with new
!DS instructions stored by reference will be emptied upon exit.
         implicit none
         integer(INTD):: ierr                      !out: error code
         class(ds_unit_t), intent(inout):: this    !inout: DS unit whose port is being loaded
         type(list_bi_t), intent(inout):: new_list !inout: list of new DS instructions for the DS unit (from other DS units): List items are stored by reference

         ierr=this%port%accept(new_list) !DS instructions (references) will be moved into the port
         return
        end function DSUnitLoadPort
!--------------------------------------------------
        function DSUnitFlushPort(this) result(ierr)
!Flushes the content of the DS unit port into the DS unit queue.
         implicit none
         integer(INTD):: ierr                   !out: error code
         class(ds_unit_t), intent(inout):: this !inout: DS unit

         ierr=this%port%absorb(this%iqueue) !DS instructions (references) will be moved from the port into the DS unit queue
         return
        end function DSUnitFlushPort
!-----------------------------------------------------
        function DSUnitGetDSVP(this,ierr) result(dsvp)
!Returns a pointer to the DSVP the DS unit is part of.
         implicit none
         class(dsvp_t), pointer:: dsvp               !out: pointer to the DSVP
         class(ds_unit_t), intent(in), target:: this !in: configured DS unit
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc

         errc=DSVP_SUCCESS
         dsvp=>this%dsvp_p; if(.not.associated(dsvp)) errc=DSVP_ERR_INVALID_REQ
         if(present(ierr)) ierr=errc
         return
        end function DSUnitGetDSVP
!------------------------------------------------------
        function DSUnitGetId(this,ierr) result(unit_id)
!Returns the DS unit id.
         implicit none
         integer(INTD):: unit_id                     !out: unit id
         class(ds_unit_t), intent(in):: this         !in: DS unit
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc

         errc=DSVP_SUCCESS; unit_id=this%id
         if(unit_id.lt.0) errc=DSVP_ERR_INVALID_REQ !no valid id
         if(present(ierr)) ierr=errc
         return
        end function DSUnitGetId
!-----------------------------------------------------
        subroutine DSUnitSetId(this,dsvp,unit_id,ierr)
!Sets DS unit id (>=0).
         implicit none
         class(ds_unit_t), intent(inout):: this      !inout: DS unit
         class(dsvp_t), intent(in), target:: dsvp    !in: DSVP the DS unit is part of
         integer(INTD), intent(in):: unit_id         !in: DS unit id: >=0
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc

         errc=DSVP_SUCCESS
         if(unit_id.ge.0) then
          this%id=unit_id
          this%dsvp_p=>dsvp
         else
          errc=DSVP_ERR_INVALID_ARGS
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine DSUnitSetId
![dsvp_t]==============================
        subroutine DSVPStart(this,ierr)
!Starts DSVP active life cycle (starts all DS units).
         implicit none
         class(dsvp_t), intent(inout):: this         !inout: configured DSVP becomes active
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: ier,errc,nthreads,mthreads,tid

         errc=DSVP_SUCCESS
         if(this%get_status().eq.DSVP_STAT_OFF.and.this%num_units.gt.0) then
          nthreads=this%num_units !number of threads needed = number of units
          mthreads=omp_get_max_threads()
          if(mthreads.lt.nthreads) then
           call omp_set_num_threads(nthreads)
           mthreads=omp_get_max_threads()
          endif
          if(mthreads.ge.nthreads) then
!$OMP PARALLEL DEFAULT(SHARED) NUM_THREADS(mthreads)
           if(omp_get_num_threads().ge.nthreads) then
!$OMP MASTER
            call this%set_status(DSVP_STAT_ON,errc)
!$OMP END MASTER
!$OMP BARRIER
            if(errc.eq.DSVP_SUCCESS) then
             write(CONS_OUT,'("#MSG(dsvp_base:dsvp_t.start): Spawned ",i5," threads")') omp_get_num_threads() !debug
             tid=omp_get_thread_num()
             call this%units(tid)%unit_ref%start(errc)
            endif
           else
            write(CONS_OUT,'("#FATAL(dsvp_base:dsvp_t.start): Unable to spawn ",i5," threads!")') nthreads
            errc=DSVP_ERR_RSC_EXCEEDED
           endif
!$OMP END PARALLEL
           call this%shutdown(ier); if(ier.ne.DSVP_SUCCESS.and.errc.eq.DSVP_SUCCESS) errc=ier
          else
           write(CONS_OUT,'("#FATAL(dsvp_base:dsvp_t.start): Insufficient number of threads: ",i5," when need ",i5)')&
           &mthreads,nthreads
           errc=DSVP_ERR_RSC_EXCEEDED
          endif
         else
          errc=DSVP_ERR_INVALID_REQ
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine DSVPStart
!-----------------------------------------
        subroutine DSVPShutdown(this,ierr)
!Shuts down DSVP.
         implicit none
         class(dsvp_t), intent(inout):: this         !inout: active DSVP becomes inactive but still configured
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc

         errc=DSVP_SUCCESS
         !`Finish
         call this%set_status(DSVP_STAT_OFF,errc)
         if(present(ierr)) ierr=errc
         return
        end subroutine DSVPShutdown
!----------------------------------------
        subroutine DSVPDestroy(this,ierr)
!Destroys inactive DSVP completely.
         implicit none
         class(dsvp_t), intent(inout):: this         !inout: inactive DSVP
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc
         integer:: ier

         errc=DSVP_SUCCESS
         if(this%get_status(errc).eq.DSVP_STAT_OFF) then
          this%microcode=>NULL()
          if(allocated(this%description)) then
           deallocate(this%description,STAT=ier)
           if(ier.ne.0.and.errc.eq.DSVP_SUCCESS) errc=DSVP_ERR_MEM_FREE_FAIL
          endif
          if(allocated(this%units)) then
           deallocate(this%units,STAT=ier)
           if(ier.ne.0.and.errc.eq.DSVP_SUCCESS) errc=DSVP_ERR_MEM_FREE_FAIL
          endif
          this%num_units=0
          this%instr_received=0
          this%instr_processed=0
          this%instr_failed=0
          this%spec_kind=DSVP_NO_KIND
          this%id=-1
         else
          errc=DSVP_ERR_INVALID_REQ
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine DSVPDestroy
!-----------------------------------------------------
        subroutine DSVPAllocUnits(this,num_units,ierr)
!Allocates the DSVU table in DSVP.
         implicit none
         class(dsvp_t), intent(inout):: this         !inout: DSVP
         integer(INTD), intent(in):: num_units       !in: DSVU table size
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc

         errc=DSVP_SUCCESS; this%num_units=0
         if(.not.allocated(this%units)) then
          if(num_units.gt.0) then
           allocate(this%units(0:num_units-1),STAT=errc)
           if(errc.ne.0) errc=DSVP_ERR_MEM_ALLOC_FAIL
          else
           errc=DSVP_ERR_INVALID_ARGS
          endif
         else
          errc=DSVP_ERR_INVALID_REQ
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine DSVPAllocUnits
!------------------------------------------
        subroutine DSVPFreeUnits(this,ierr)
!Frees the DSVU table.
         implicit none
         class(dsvp_t), intent(inout):: this         !inout: DSVP
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc

         errc=DSVP_SUCCESS
         if(allocated(this%units)) then
          deallocate(this%units,STAT=errc)
          if(errc.ne.0) errc=DSVP_ERR_MEM_FREE_FAIL
         else
          errc=DSVP_ERR_INVALID_REQ
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine DSVPFreeUnits
!--------------------------------------------------
        subroutine DSVPSetUnit(this,virt_unit,ierr)
!Sets up a new DSVU entry in the DSVU table in a sequential order.
         implicit none
         class(dsvp_t), intent(inout), target:: this       !inout: DSVP
         class(ds_unit_t), intent(in), pointer:: virt_unit !in: valid DSVU
         integer(INTD), intent(out), optional:: ierr       !out: error code
         integer(INTD):: errc

         errc=DSVP_SUCCESS
         if(associated(virt_unit)) then
          if(this%num_units.le.ubound(this%units,1)) then
           this%units(this%num_units)%unit_ref=>virt_unit
           call this%units(this%num_units)%unit_ref%set_id(this,this%num_units,errc)
           if(errc.eq.DSVP_SUCCESS) then
            this%num_units=this%num_units+1
           else
            this%units(this%num_units)%unit_ref=>NULL()
           endif
          else
           errc=DSVP_ERR_RSC_EXCEEDED
          endif
         else
          errc=DSVP_ERR_INVALID_ARGS
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine DSVPSetUnit
!-------------------------------------------------------------
        function DSVPGetUnit(this,unit_id,ierr) result(dsvu_p)
!Returns a pointer to the DSVU from the DSVU table.
         implicit none
         class(ds_unit_t), pointer:: dsvu_p          !out: pointer to DSVU
         class(dsvp_t), intent(in):: this            !in: DSVP
         integer(INTD), intent(in):: unit_id         !in: DSVU id in the DSVU table
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc

         errc=DSVP_SUCCESS; dsvu_p=>NULL()
         if(unit_id.ge.0.and.unit_id.lt.this%num_units) then
          dsvu_p=>this%units(unit_id)%unit_ref
         else
          errc=DSVP_ERR_INVALID_ARGS
         endif
         if(present(ierr)) ierr=errc
         return
        end function DSVPGetUnit
!------------------------------------------------------------------
        subroutine DSVPSetDescription(this,id,descr,ierr,spec_kind)
!Sets the DSVP ID, kind, and symbolic description.
         implicit none
         class(dsvp_t), intent(inout):: this             !inout: active DSVP
         integer(INTL), intent(in):: id                  !in: unique ID: >=0
         character(*), intent(in):: descr                !in: symbolic description
         integer(INTD), intent(out), optional:: ierr     !out: error code
         integer(INTD), intent(in), optional:: spec_kind !in: specific kind of DSVP: >=0
         integer(INTD):: errc,l
         integer:: ier

         errc=DSVP_SUCCESS
         if(this%get_status(errc).ne.DSVP_STAT_OFF) then
          if(allocated(this%description)) then
           deallocate(this%description,STAT=ier)
           if(ier.ne.0.and.errc.eq.DSVP_SUCCESS) errc=DSVP_ERR_MEM_FREE_FAIL
          endif
          l=len_trim(descr)
          if(l.gt.0) then
           allocate(character(len=l)::this%description,STAT=ier)
           if(ier.eq.0) then
            this%description(1:l)=descr(1:l)
           else
            if(errc.eq.DSVP_SUCCESS) errc=DSVP_ERR_MEM_ALLOC_FAIL
           endif
          endif
          this%id=id
          if(present(spec_kind)) this%spec_kind=spec_kind
         else
          errc=DSVP_ERR_INVALID_REQ
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine DSVPSetDescription
!----------------------------------------------------------------------------
        subroutine DSVPGetDescription(this,id,descr,descr_len,ierr,spec_kind)
!Gets the DSVP ID, kind, and symbolic description.
         implicit none
         class(dsvp_t), intent(in):: this                 !in: active DSVP
         integer(INTL), intent(out):: id                  !out: DSVP ID
         character(*), intent(inout):: descr              !out: symbolic description
         integer(INTD), intent(out):: descr_len           !out: length of the description string
         integer(INTD), intent(out), optional:: ierr      !out: error code
         integer(INTD), intent(out), optional:: spec_kind !out: specific kind of DSVP
         integer(INTD):: errc,l

         errc=DSVP_SUCCESS
         if(this%get_status(errc).ne.DSVP_STAT_OFF) then
          if(allocated(this%description)) then
           descr_len=len_trim(this%description)
           if(descr_len.gt.0) then
            if(len(descr).ge.descr_len) then
             descr(1:descr_len)=this%description(1:descr_len)
            else
             errc=DSVP_ERR_INVALID_ARGS !insufficient length of <descr>
            endif
           endif
          else
           descr_len=0
          endif
          id=this%id
          if(present(spec_kind)) spec_kind=this%spec_kind
         else
          errc=DSVP_ERR_INVALID_REQ
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine DSVPGetDescription
!----------------------------------------------------------
        subroutine DSVPIncrRecvInstrCounter(this,ierr,incr)
!Increments the received instruction counter.
         implicit none
         class(dsvp_t), intent(inout):: this         !inout: DSVP
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTL), intent(in), optional:: incr  !in: specific increment
         integer:: errc

         errc=DSVP_SUCCESS
         if(this%get_status(errc).ne.DSVP_STAT_OFF) then
          if(errc.eq.DSVP_SUCCESS) then
           if(present(incr)) then
            this%instr_received=this%instr_received+incr
           else
            this%instr_received=this%instr_received+1_INTL
           endif
          endif
         else
          errc=DSVP_ERR_INVALID_REQ
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine DSVPIncrRecvInstrCounter
!----------------------------------------------------------
        subroutine DSVPIncrRtrdInstrCounter(this,ierr,incr)
!Increments the retired instruction counter.
         implicit none
         class(dsvp_t), intent(inout):: this         !inout: DSVP
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTL), intent(in), optional:: incr  !in: specific increment
         integer:: errc

         errc=DSVP_SUCCESS
         if(this%get_status(errc).ne.DSVP_STAT_OFF) then
          if(errc.eq.DSVP_SUCCESS) then
           if(present(incr)) then
            this%instr_processed=this%instr_processed+incr
           else
            this%instr_processed=this%instr_processed+1_INTL
           endif
          endif
         else
          errc=DSVP_ERR_INVALID_REQ
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine DSVPIncrRtrdInstrCounter
!----------------------------------------------------------
        subroutine DSVPIncrFailInstrCounter(this,ierr,incr)
!Increments the failed instruction counter.
         implicit none
         class(dsvp_t), intent(inout):: this         !inout: DSVP
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTL), intent(in), optional:: incr  !in: specific increment
         integer:: errc

         errc=DSVP_SUCCESS
         if(this%get_status(errc).ne.DSVP_STAT_OFF) then
          if(errc.eq.DSVP_SUCCESS) then
           if(present(incr)) then
            this%instr_failed=this%instr_failed+incr
           else
            this%instr_failed=this%instr_failed+1_INTL
           endif
          endif
         else
          errc=DSVP_ERR_INVALID_REQ
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine DSVPIncrFailInstrCounter
!----------------------------------------------------
        function DSVPTimeActive(this,ierr) result(tm)
!Returns the time DSVP is active in seconds.
         implicit none
         real(8):: tm                                !out: time active in seconds
         class(dsvp_t), intent(in):: this            !in: DSVP
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc

         errc=DSVP_SUCCESS; tm=-1d0
         if(this%get_status(errc).ne.DSVP_STAT_OFF) then
          if(errc.eq.DSVP_SUCCESS) tm=thread_wtime(this%time_start)
         else
          errc=DSVP_ERR_INVALID_REQ
         endif
         if(present(ierr)) ierr=errc
         return
        end function DSVPTimeActive
!----------------------------------------------------
        function DSVPGetStatus(this,ierr) result(sts)
!Returns the current status of the DSVP.
         implicit none
         integer(INTD):: sts                         !out: DSVP current status
         class(dsvp_t), intent(in):: this            !in: DSVP
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc

         errc=DSVP_SUCCESS
         sts=this%stat
         if(present(ierr)) ierr=errc
         return
        end function DSVPGetStatus
!----------------------------------------------
        subroutine DSVPSetStatus(this,sts,ierr)
!Sets the DSVP status.
         implicit none
         class(dsvp_t), intent(inout):: this         !inout: DSVP
         integer(INTD), intent(in):: sts             !in: status
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc

         errc=DSVP_SUCCESS
         this%stat=sts; if(this%stat.eq.DSVP_STAT_ON) call this%start_time(errc)
         if(present(ierr)) ierr=errc
         return
        end subroutine DSVPSetStatus
!------------------------------------------
        subroutine DSVPStartTime(this,ierr)
!Starts the time after initializing the DSVP.
         implicit none
         class(dsvp_t), intent(inout):: this         !inout: DSVP
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc

         errc=DSVP_SUCCESS
         this%time_start=thread_wtime()
         if(present(ierr)) ierr=errc
         return
        end subroutine DSVPStartTime

       end module dsvp_base
