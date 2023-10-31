#include <xxhash.h>
#include <string.h>
#include <stdio.h>

void xxhash(const char* data, unsigned long size, unsigned char *result)
{

   XXH128_hash_t tmp;
   
   tmp = XXH3_128bits(data,size);
/* This does not preserve endianness but we only need a consistent one */
   memcpy ( result, &tmp, 16 );
}

