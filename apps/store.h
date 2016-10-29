#ifndef _STORE_H_
#define _STORE_H_

inline
symmetricVertex* getStore(uint64_t i) {
    return objStore[i];
}

inline void putStore(uint64_t i, symmetricVertex* t) {
    objStore[i] = t;
}

#endif // _STORE_H_
