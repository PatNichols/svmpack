/*
 * SVMBarrier.h
 *
 *  Created on: Jul 13, 2010
 *      Author: d3p708
 */

#ifndef SVMBARRIERPTH_H_
#define SVMBARRIERPTH_H_
#include <cstdlib>
#include <pthread.h>

using namespace std;

namespace svmpack
{
//#define _USE_TBARRIER_
#ifdef _USE_TBARRIER_
template <class C_t>
struct CyclicBarrierPth {
private:
    mem_fun_t<void, C_t> myfun;
    pthread_mutex_t *mutex;
    pthread_cond_t *cond;
    volatile int *cnt;
    int *block_cnt;
    int nth, tid;
    int nlevels, nnodes;
public:
    CyclicBarrierPth ( CyclicBarrierPth<C_t>& b, int threadId ) :
        myfun( b.myfun ), mutex(b.mutex), cond(b.cond), cnt(b.cnt), block_cnt(b.block_cnt),
        nth(b.nth), tid(threadId), nlevels( b.nlevels), nnodes(b.nnodes)  {};

    CyclicBarrierPth ( int numThreads, mem_fun_t<void, C_t> fun ) :
        myfun(fun), mutex(0x0), cond(0x0), cnt(0x0), block_cnt(0x0),
        nth(numThreads), tid(0), nlevels(0), nnodes(0) {
        int block_size = 2;
        if ( nth > 2 ) {
            while ( block_size <= nth ) {
                ++nlevels;
                nnodes += nth / block_size;
                block_size <<= 1;
            }
            cnt = new int[nnodes];
            block_cnt = new int[nnodes];
            mutex = new pthread_mutex_t[nnodes];
            cond = new pthread_cond_t[nnodes];
            int k;
            int offset = 0;
            int npt=nth;
            block_size = 2;
            for ( int i = 0; i < nlevels; ++i ) {
                int level_size = nth / block_size;
                for ( k = 0; k < level_size; ++k ) {
                    pthread_mutex_init ( mutex + offset, 0x0 );
                    pthread_cond_init ( cond + offset, 0x0 );
                    cnt[offset] = 2;
                    block_cnt[offset] = 2;
                    ++offset;
                }
                if ( npt % 2 ) {
                    k=offset-1;
                    cnt[k] += 1;
                    block_cnt[k] += 1;
                }
                npt= level_size;
                block_size <<= 1;
            }
        } else {
            nnodes = 1;
            cnt = new int[nnodes];
            block_cnt = new int[nnodes];
            mutex = new pthread_mutex_t[nnodes];
            cond = new pthread_cond_t[nnodes];
            pthread_mutex_init ( mutex, 0x0 );
            pthread_cond_init ( cond, 0x0 );
            cnt[0] = nth;
            block_cnt[0] = nth;
        }
    };

    ~CyclicBarrierPth() {
        if ( tid == 0 ) {
            for ( int k = 0; k < nnodes; ++k ) {
                pthread_mutex_destroy ( mutex + k );
                pthread_cond_destroy ( cond + k );
            }
            delete [] block_cnt;
            delete [] cnt;
            delete [] cond;
            delete [] mutex;
        }
    };

    void reduce ( C_t *c ) {
        int winner = 1;
        int level_base = 0;
        int level_size, level_id;
        int block_size = 2;
        while ( block_size <= nth ) {
            level_size = nth / block_size;
            level_id = tid / block_size;
            if ( level_id >= level_size ) --level_id;
            level_id += level_base;
            pthread_mutex_lock ( mutex + level_id );
            int r = --cnt[level_id];
            if ( r ) {
                // we lost
                winner = 0;
                pthread_cond_wait ( cond + level_id, mutex + level_id );
                pthread_mutex_unlock ( mutex + level_id );
                break;
            } else {
                // we won this round reset count to proper value for next round
                cnt[level_id] = block_cnt[level_id];
            }
            level_base += level_size;
            block_size <<= 1;
        }
        if ( winner ) myfun ( c );
        // we need to release all of the losers below us
        block_size >>= 1;
        level_size = nth / block_size;
        level_base -= level_size;
        while ( block_size > 1 ) {
            level_id = tid / block_size;
            if ( level_id >= level_size ) --level_id;
            level_id += level_base;
            pthread_cond_broadcast ( cond + level_id );
            pthread_mutex_unlock ( mutex + level_id );
            block_size >>= 1;
            level_size = nth / block_size;
            level_base -= level_size;
        }
    };

    void await ( ) {
        int winner = 1;
        int level_base = 0;
        int level_size, level_id;
        int block_size = 2;
        while ( block_size <= nth ) {
            level_size = nth / block_size;
            level_id = tid / block_size;
            if ( level_id >= level_size ) --level_id;
            level_id += level_base;
            pthread_mutex_lock ( mutex + level_id );
            int r = --cnt[level_id];
            if ( r ) {
                // we lost
                winner = 0;
                pthread_cond_wait ( cond + level_id, mutex + level_id );
                pthread_mutex_unlock ( mutex + level_id );
                break;
            } else {
                // we won this round reset count to proper value for next round
                cnt[level_id] = block_cnt[level_id];
            }
            block_size <<= 1;
            level_base += level_size;
        }
        // we need to release all of the losers below us
        block_size >>= 1;
        level_size = nth / block_size;
        level_base -= level_size;
        while ( block_size > 1 ) {
            level_id = tid / block_size;
            if ( level_id >= level_size ) --level_id;
            level_id += level_base;
            pthread_cond_broadcast ( cond + level_id );
            pthread_mutex_unlock ( mutex + level_id );
            block_size >>= 1;
            level_size = nth / block_size;
            level_base -= level_size;
        }
    };

};
#else
template <class c_t> struct CyclicBarrierPth {
    mem_fun_t<void, c_t> fun;
    pthread_mutex_t *mutex;
    pthread_cond_t *cond;
    volatile int *cnt;
    int nth;
    int tid;

    CyclicBarrierPth ( int nparty, mem_fun_t<void, c_t>  mfun_in ) :
        fun ( mfun_in ), mutex ( new pthread_mutex_t() ),
        cond ( new pthread_cond_t() ),
        cnt ( new int ( nparty ) ), nth ( nparty ), tid ( 0 ) {
        pthread_mutex_init ( mutex, 0x0 );
        pthread_cond_init ( cond, 0x0 );
    };

    CyclicBarrierPth ( CyclicBarrierPth<c_t>& c, int thID ) :
        fun ( c.fun ), mutex ( c.mutex ), cond ( c.cond ),
        cnt ( c.cnt ), nth ( c.nth ), tid ( thID ) {
    };

    ~CyclicBarrierPth() {
        if ( !tid ) {
            pthread_cond_destroy ( cond );
            pthread_mutex_destroy ( mutex );
            delete cnt;
            delete cond;
            delete mutex;
        }
    };

    void reduce ( c_t * t ) {
        pthread_mutex_lock ( mutex );
        if ( --*cnt ) pthread_cond_wait ( cond, mutex );
        else {
            *cnt = nth;
            fun ( t );
            pthread_cond_broadcast ( cond );
        }
        pthread_mutex_unlock ( mutex );
    };
    void await () {
        pthread_mutex_lock ( mutex );
        if ( --*cnt ) pthread_cond_wait ( cond, mutex );
        else {
            *cnt = nth;
            pthread_cond_broadcast ( cond );
        }
        pthread_mutex_unlock ( mutex );
    };
};
#endif

}
#endif /* SVMBARRIER_H_ */
