!Domain-specific virtual processor (DSVP): Abstract base module.
!This module provides the infrastructure for the tensor algebra processor.
!AUTHOR: Dmitry I. Lyakh (Liakh): quant4me@gmail.com
!REVISION: 2016/11/09

!Copyright (C) 2014-2016 Dmitry I. Lyakh (Liakh)
!Copyright (C) 2014-2016 Oak Ridge National Laboratory (UT-Battelle)

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

       module dsvp_base
        use dil_basic
        implicit none
        private
!PARAMETERS:
 !Basic:
        integer(INTD), private:: CONS_OUT=6 !default output for this module
        integer(INTD), private:: DEBUG=0    !debugging mode
        logical, private:: VERBOSE=.true.   !verbosity for errors
 !Error codes:
        integer(INTD), parameter, public:: DSVP_SUCCESS=SUCCESS       !success
        integer(INTD), parameter, public:: DSVP_ERROR=GENERIC_ERROR   !generic error code
        integer(INTD), parameter, public:: DSVP_ERROR_INVALID_ARGS=1  !specific error code: invalid procedure arguments
 !DSVP status:
        integer(INTD), parameter, public:: DSVP_STAT_CLEAN=0  !DSVP has not been initialized
        integer(INTD), parameter, public:: DSVP_STAT_ACTIVE=1 !DSVP has been initialized and is active now
        integer(INTD), parameter, public:: DSVP_STAT_ERROR=2  !DSVP encountered an error
!DERIVED TYPES:
#if 0
 !Domain-specific instruction:
        type, abstract, public:: ds_instr_t
         
         contains
        end type ds_instr_t
 !Domain-specific virtual processor:
        type, abstract, public:: dsvp_t
         integer(INTL), private:: id=-1                       !DSVP unique ID
         integer(INTD), private:: stat=DSVP_STAT_CLEAN        !current DSVP status
         integer(INTL), private:: instr_received=0_INTL       !total number of received instructions
         integer(INTL), private:: instr_processed=0_INTL      !total number of processed (retired) instructions
         integer(INTL), private:: instr_failed=0_INTL         !total number of failed instructions
         real(8), private:: time_start                        !start time
         contains
          procedure, public:: init=>DSVPInit                                                         !initializes DSVP to an active state
          procedure, public:: shutdown=>DSVPShutdown                                                 !shutdowns DSVP
          procedure, public:: incr_instruction_counter=>DSVPIncrInstructionCounter                   !increments the receieved instruction counter
          procedure, public:: incr_retired_counter=>DSVPIncrRetiredCounter                           !indrements the processed (retired) instruction counter
          procedure, public:: incr_error_counter=>DSVPIncrErrorCounter                               !increments the failed instruction counter
          procedure, public:: time_active=>DSVPTimeActive                                            !returns the time DSVP is active in seconds
          procedure, deferred, public:: fetch_instructions=>DSVPFetchInstructions                    !fetches a block of domain-specific instructions from another DSVP
          procedure, deferred, public:: return_retired_instructions=>DSVPReturnRetiredInstructions   !returns back a block retired instructions with their statuses
          procedure, deferred, public:: send_instructions=>DSVPSendInstructions                      !sends a block of domain-specific instructions to another DSVP for execution
          procedure, deferred, public:: receive_retired_instructions=>DSVPReceiveRetiredInstructions !receives back a block of retired instructions with their statuses
        end type dsvp_t
!DATA:

!VISIBILITY:
#endif
       contains
!IMPLEMENTATION:
!-------------------------------------------------

       end module dsvp_base
