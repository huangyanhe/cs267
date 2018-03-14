#pragma once

#include <upcxx/upcxx.hpp>
#include "kmer_t.hpp"

struct HashMap {
  std::vector <kmer_pair> data;
  std::vector <int> used;

//  // pointer to an array of pointers that point to local shared used
//  global_ptr<int> ptr2used; 
  std::vector<upcxx::global_ptr<int> > ptr2used;    
  std::vector<upcxx::global_ptr<kmer_pair> > ptr2used;    
  
  size_t my_size;

  size_t size() const noexcept;

  HashMap(size_t size, int &ptr);

  // Most important functions: insert and retrieve
  // k-mers from the hash table.
  bool insert(const kmer_pair &kmer);
  bool find(const pkmer_t &key_kmer, kmer_pair &val_kmer);

  // Helper functions
  
  // Write and read to a logical data slot in the table.
  //void write_slot(uint64_t slot, const kmer_pair &kmer);
  void HashMap::write_slot(kmer_pair *slot, const kmer_pair &kmer);
  //kmer_pair read_slot(uint64_t slot);
  kmer_pair read_slot(kmer_pair *slot);
  // Request a slot or check if it's already used.
  //bool request_slot(uint64_t slot);
  bool HashMap::request_slot(int *slot);
  bool HashMap::slot_used(int *slot);
 //bool slot_used(uint64_t slot);
  //My own added functions
  int64_t HashMap::refdata();
  int64_t HashMap::refused();
  int HashMap::initPtrs(std::vector<upcxx::global_ptr<int> > ptrInitused, std::vector<upcxx::global_ptr<kmer_pair> >ptrInitdata);

};

HashMap::HashMap(size_t size) {
  my_size = size;
  data.resize(size);
  used.resize(size, 0);
}

bool HashMap::insert(const kmer_pair &kmer) {
  uint64_t hash = kmer.hash();
  uint64_t probe = 0;
  bool success = false;
  do {
    uint64_t slot = (hash + probe++) % size();
    uint64_t IndexRank = std::floor(slot/size);
    uint64_t IndexOnLocalVector = slot- IndexRank*size;
    int *ptr2writeused = ptr2used[IndexRank] +  IndexOnLocalVector; 
    success = request_slot(ptr2writeused);
    kmer_pair *ptr2writedata = ptr2data[IndexRank] +  IndexOnLocalVector; 
    if (success) {
      write_slot(ptr2writedata, kmer);
    }
  } while (!success && probe < size());
  return success;
}

bool HashMap::find(const pkmer_t &key_kmer, kmer_pair &val_kmer) {
  uint64_t hash = key_kmer.hash();
  uint64_t probe = 0;
  bool success = false;
  do {
    uint64_t slot = (hash + probe++) % size();
    uint64_t IndexRank = std::floor(slot/size);
    uint64_t IndexOnLocalVector = slot- IndexRank*size;
    int *ptr2writeused = ptr2used[IndexRank] +  IndexOnLocalVector; 
    if (slot_used(ptr2writeused)) {
        kmer_pair *ptr2writedata = ptr2data[IndexRank] +  IndexOnLocalVector; 
        val_kmer = read_slot(ptr2writedata);
      if (val_kmer.kmer == key_kmer) {
        success = true;
      }
    }
  } while (!success && probe < size());
  return success;
}

//bool HashMap::slot_used(uint64_t slot) {
//  return used[slot] != 0;
//}
//
//void HashMap::write_slot(uint64_t slot, const kmer_pair &kmer) {
//  data[slot] = kmer;
//}
//
//kmer_pair HashMap::read_slot(uint64_t slot) {
//  return data[slot];
//}
bool HashMap::slot_used(int *slot) {
  return *slot != 0;
}

void HashMap::write_slot(kmer_pair *slot, const kmer_pair &kmer) {
  *slot = kmer;
}

kmer_pair HashMap::read_slot(kmer_pair *slot) {
  return *slot;
}
//bool HashMap::request_slot(uint64_t slot) {
//  if (used[slot] != 0) {
//    return false;
//  } else {
//    used[slot] = 1;
//    return true;
//  }
//  }
bool HashMap::request_slot(int *slot) {
int result = upcxx::atomic_fetch_add(slot, 1, memory_order_relaxed).wait()
if (result == 0)
{
    return true; 
}
else
{
    return false;
}

}

size_t HashMap::size() const noexcept {
  return my_size;
}
int HashMap::refused()
{
return used[0];
}
kmer_pair HashMap::refdata()
{
return data[0];
}
int HashMap::initPtrs(vector<upcxx::global_ptr<int> > ptrInitused, vector<upcxx::global_ptr<kmer_pair> > ptrInitdata) {
    ptr2used = ptrInitused;    
    ptr2data = ptrInitdata; 
    return 0; 
}
