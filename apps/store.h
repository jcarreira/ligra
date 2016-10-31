#ifndef _STORE_H_
#define _STORE_H_

#define MAX 10000000
    
struct Store {
    symmetricVertex** objStore;
    char* bitset;
};

extern Store my_store;

inline
void init_store() {
    my_store.objStore = new symmetricVertex*[MAX];
    my_store.bitset = new char[MAX];
}

inline
symmetricVertex* getStore(uint64_t i) {
    return my_store.objStore[i];
}

inline void putStore(uint64_t i, symmetricVertex* t) {
    my_store.objStore[i] = t;
}

inline void deleteStore(Store*store, uint64_t i) {
}

#endif // _STORE_H_
