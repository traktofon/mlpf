/*
 * Fortran interfaces to some POSIX I/O routines.
 * In case you want to avoid the record markers.
 * 
 * FO 02/2013
 */


#if (F2C_STYLE == 1)
#   define C_OPEN c_open
#   define C_CLOSE c_close
#   define C_WRITE_AT c_write_at
#elif (F2C_STYLE == 2)
#   define C_OPEN C_OPEN
#   define C_CLOSE C_CLOSE
#   define C_WRITE_AT C_WRITE_AT
#else
#   define C_OPEN c_open_
#   define C_CLOSE c_close_
#   define C_WRITE_AT c_write_at_
#endif


#include <sys/mman.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>


/*---------------------------------------------------------------------- 
 *
 * Opening a file.
 * Access mode is determined as follows:
 *   lrdwr = 0 : read-only
 *   lrdwr = 1 : read-write, create if necessary
 *   lrdwr = 2 : read-write, create if necessary, truncate existing file
 * 
 * FORTRAN:
 *    integer          :: fd,lrdwr
 *    character(len=*) :: fname
 *    call c_open(fd, fname, lrdwr)
 */
void C_OPEN(int* fd, const char* fname, int* lrdwr, int flen) {
    int    flags;
    mode_t mode;
    char*  buf;

    switch (*lrdwr) {
        case 1:
            flags = O_RDWR|O_CREAT;
            break;
        case 2:
            flags = O_RDWR|O_CREAT|O_TRUNC;
            break;
        default:
            flags = O_RDONLY;
    }
    mode = S_IRUSR|S_IWUSR|S_IRGRP|S_IROTH;
    buf = (char*) malloc(flen+1);
    memcpy(buf, fname, flen);
    buf[flen] = '\0';
    *fd = open(buf, flags, mode);
    if (*fd == -1) {
        perror(buf);
        exit(EXIT_FAILURE);
    }
    return;
}


/*----------------------------------------------------------------------
 *
 * Closing a file.
 *
 * FORTRAN:
 *    integer :: fd
 *    call c_close(fd)
 */
void C_CLOSE(int* fd) {
    int ret;
    ret = close(*fd);
    if (ret != 0) {
        perror("close");
        exit(EXIT_FAILURE);
    }
    return;
}


/*----------------------------------------------------------------------
 *
 * Writing to a file at a given position.
 * The file will be expanded as needed.
 *
 * Writes _nelem_ elements, each of size _size_, from the buffer _buf_
 * to the file specified by _fd_ at the position _where_.
 * _where_ is measured in elements, so the actual file offset in bytes
 * is _where_*_size_ .
 *
 * FORTRAN:
 *    integer   :: fd, size
 *    integer*8 :: where, nelem
 *    type(*)   :: buf(*)
 *    call c_write_at(fd, buf, size, where, nelem)
 */
void C_WRITE_AT(int* fd, void* buf, int* size, int64_t* where, int64_t* nelem) {
    off_t   offset, offret;
    size_t  count;
    ssize_t nwritten;

    offset = (*where)*(*size);
    offret = lseek(*fd, offset, SEEK_SET);
    if (offret != offset) {
        perror("lseek");
        exit(EXIT_FAILURE);
    }
    
    count = (*nelem)*(*size);
    nwritten = write(*fd, buf, count);
    if (nwritten != count) {
        perror("write");
        exit(EXIT_FAILURE);
    }
    return;
}


/*----------------------------------------------------------------------
 *
 * Mapping a file into memory.
 *
 * Maps _length_ bytes at the given _offset_ from the file specified
 * by _fd_ into memory and returns a pointer to the mapped area. The
 * access mode is specified by _lrdwr_:
 *   lrdwr = 0 : read-only
 *   lrdwr = 1 : read-write
 *
 * Needs a Fortran2003 interface block to be called from Fortran,
 * see mmap_m.f90.
 */
void* c_mmap(int fd, int lrdwr, size_t length, size_t offset) {
    void  *addr;
    int    prot, flags;

    switch (lrdwr) {
        case 1:
            prot  = PROT_READ | PROT_WRITE;
            flags = MAP_SHARED;
            break;
        default:
            prot  = PROT_READ;
            flags = MAP_SHARED;
    }
    addr = mmap(NULL, length, prot, flags, fd, offset);
    if (addr == MAP_FAILED) {
        perror("mmap");
        exit(EXIT_FAILURE);
    }
    return addr;
}

