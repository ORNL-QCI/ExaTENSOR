/** TAL-SH: Byte packet
REVISION: 2019/02/26

Copyright (C) 2018-2019 Dmitry I. Lyakh (Liakh)
Copyright (C) 2018-2019 Oak Ridge National Laboratory (UT-Battelle) **/

#include "byte_packet.h"

#include <assert.h>
#include <cstdlib>

#ifdef __cplusplus
#include <cstddef>
#endif

void initBytePacket(BytePacket * packet)
{
 packet->capacity = 0;
 packet->base_addr = malloc(BYTE_PACKET_CAPACITY);
 if(packet->base_addr != NULL) packet->capacity = BYTE_PACKET_CAPACITY;
 packet->size_bytes = 0;
 packet->position = 0;
 return;
}

void clearBytePacket(BytePacket * packet)
{
 packet->capacity = 0;
 free(packet->base_addr);
 packet->base_addr = NULL;
 packet->size_bytes = 0;
 packet->position = 0;
 return;
}

void resetBytePacket(BytePacket * packet)
{
 packet->position = 0;
 return;
}

#ifdef __cplusplus
template <typename T>
void appendToBytePacket(BytePacket * packet, const T & item)
{
 char * dst_ptr = &(((char*)(packet->base_addr))[packet->position]);
 char * src_ptr = ((char*)(&item));
 unsigned long long type_size = sizeof(T);
 assert(packet->position + type_size <= packet->capacity);
 for(unsigned long long i = 0; i < type_size; ++i) dst_ptr[i] = src_ptr[i];
 packet->position += type_size;
 if(packet->position > packet->size_bytes) packet->size_bytes = packet->position;
 return;
}

template <typename T>
void extractFromBytePacket(BytePacket * packet, T & item)
{
 char * src_ptr = &(((char*)(packet->base_addr))[packet->position]);
 char * dst_ptr = ((char*)(&item));
 unsigned long long type_size = sizeof(T);
 assert(packet->position + type_size <= packet->size_bytes);
 for(unsigned long long i = 0; i < type_size; ++i) dst_ptr[i] = src_ptr[i];
 packet->position += type_size;
 return;
}
#endif
