#pragma once

#include <upcxx/upcxx.hpp>
#include <vector>
#include "kmer_t.hpp"

struct HashMap {
  //std::vector <kmer_pair> data;
  //std::vector <int> used;
  std::vector<upcxx::global_ptr<kmer_pair> > data;
  std::vector<upcxx::global_ptr<int> > used;

  size_t my_size;
  size_t total_size;

  size_t local_size() const noexcept;
  size_t size() const noexcept;

  HashMap(size_t size);
  ~HashMap();

  // Most important functions: insert and retrieve
  // k-mers from the hash table.
  bool insert(const kmer_pair &kmer);
  bool find(const pkmer_t &key_kmer, kmer_pair &val_kmer);

  // Helper functions

  // Write and read to a logical data slot in the table.
  void write_slot(uint64_t slot, const kmer_pair &kmer);
  kmer_pair read_slot(uint64_t slot);

  // Request a slot or check if it's already used.
  bool request_slot(uint64_t slot);
  bool slot_used(uint64_t slot);
};

HashMap::HashMap(size_t size) {
    total_size = size * upcxx::rank_n();
    my_size = (total_size + upcxx::rank_n() - 1)/upcxx::rank_n();
    total_size = my_size * upcxx::rank_n();
    data.assign(upcxx::rank_n(), nullptr);
    used.assign(upcxx::rank_n(), nullptr);
    data[upcxx::rank_me()] = upcxx::new_array<kmer_pair>(my_size);
    used[upcxx::rank_me()] = upcxx::new_array<int>(my_size);
    upcxx::barrier();
    int total_rank = upcxx::rank_n();
    for(int i = 0; i < total_rank; i++){
        data[i] = upcxx::broadcast(data[i], i).wait();
        used[i] = upcxx::broadcast(used[i], i).wait();
    }
}

HashMap::~HashMap(){
    upcxx::delete_array(data[upcxx::rank_me()]);
    upcxx::delete_array(used[upcxx::rank_me()]);
    int total_rank = upcxx::rank_n();
    for(int i = 0; i < total_rank; i++){
        data[i] = nullptr;
        used[i] = nullptr;
    }
}

bool HashMap::insert(const kmer_pair &kmer) {
  uint64_t hash = kmer.hash();
  uint64_t probe = 0;
  bool success = false;
  do {
    uint64_t slot = (hash + probe++) % size();
    success = request_slot(slot);
    if (success) {
      write_slot(slot, kmer);
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
    if (slot_used(slot)) {
      val_kmer = read_slot(slot);
      if (val_kmer.kmer == key_kmer) {
        success = true;
      }
    }
  } while (!success && probe < size());
  return success;
}

bool HashMap::slot_used(uint64_t slot) {
  //return used[slot] != 0;
    int rank = slot/my_size, local_slot = slot%my_size;
    upcxx::future<int> f = upcxx::rget(used[rank]+local_slot);
    f.wait();
    return f.result() != 0;
}

void HashMap::write_slot(uint64_t slot, const kmer_pair &kmer) {
  //data[slot] = kmer;
    int rank = slot/my_size, local_slot = slot%my_size;
    upcxx::rput(kmer, data[rank]+local_slot).wait();
}

kmer_pair HashMap::read_slot(uint64_t slot) {
  //return data[slot];
    int rank = slot/my_size, local_slot = slot%my_size;
    upcxx::future<kmer_pair> f = upcxx::rget(data[rank]+local_slot);
    f.wait();
    return f.result();
}

bool HashMap::request_slot(uint64_t slot) {
  //if (used[slot] != 0) {
    //return false;
  //} else {
    //used[slot] = 1;
    //return true;
  //}
    int rank = slot/my_size, local_slot = slot%my_size;
    upcxx::future<int> f = upcxx::atomic_fetch_add(used[rank]+local_slot, 1, std::memory_order_relaxed);
    f.wait();
    return f.result() == 0;
}

size_t HashMap::local_size() const noexcept{
  return my_size;
}

size_t HashMap::size() const noexcept{
  return total_size;
}
