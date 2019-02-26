/** TAL-SH: Byte packet
REVISION: 2019/02/26

Copyright (C) 2018-2019 Dmitry I. Lyakh (Liakh)
Copyright (C) 2018-2019 Oak Ridge National Laboratory (UT-Battelle) **/

#ifndef BYTE_PACKET_H_
#define BYTE_PACKET_H_

#define BYTE_PACKET_CAPACITY 1048576 //default byte packet capacity in bytes

//Byte packet (interoperable):
typedef struct{
 void * base_addr;              //base address (owning pointer)
 unsigned long long capacity;   //byte packet capacity in bytes
 unsigned long long size_bytes; //actual size of the byte packet in bytes
 unsigned long long position;   //current position inside the byte packet
} BytePacket;


void initBytePacket(BytePacket * packet);
void clearBytePacket(BytePacket * packet);
void resetBytePacket(BytePacket * packet);
#ifdef __cplusplus
template <typename T> void appendToBytePacket(BytePacket * packet, const T & item);
template <typename T> void extractFromBytePacket(BytePacket * packet, T & item);
#endif

#endif //BYTE_PACKET_H_
