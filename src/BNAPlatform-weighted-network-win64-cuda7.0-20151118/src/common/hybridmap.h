#ifndef _HYBRIDMAP_H_
#define _HYBRIDMAP_H_
#include <stdlib.h>
#include <iostream>
#include "hashset.h"
#include "bitmap.h"
#include "data_type.h"
using namespace std;

typedef HashSet<_ULonglong> HashGraph;

class HybridMap {
private:
	double loadfactor_inv_;
	u_int num_vertices;
	BitMap *bitmap_;
	HashGraph *hashmap_;
	vector<u_int> index_;
	
public:
	HybridMap():bitmap_(NULL),hashmap_(NULL),index_(vector<u_int>(0)){}
	
	HybridMap(R_type *R, C_type *C, size_t N, size_t Clength, double loadfactor_inv = 2.0): 
		num_vertices(N), index_(vector<u_int>(N,N)), loadfactor_inv_(loadfactor_inv)
	{		
		u_int idx = 0;
		size_t hashmap_size = Clength;
		for(size_t i = 0; i < N; i++)
		{
			int deg = R[i+1]-R[i];			
			if (deg*loadfactor_inv_*sizeof(_ULonglong) > ((N-1)/8+1)){ 
				index_[i] = idx++;	
				hashmap_size -= deg;
			}
		}
		hashmap_ = new HashGraph(hashmap_size, loadfactor_inv_);
		size_t bitmap_size = ((size_t)idx*N - 1)/8 +1;
		bitmap_ = new BitMap(bitmap_size);
		for(size_t i = 0; i < N; i++)
		{
			if(index_[i] != N){
				for(size_t j = R[i]; j < R[i+1]; j++)
					bitmap_->bitmapSet(index_[i]*N+C[j]);
			}
			else{
				for(size_t j = R[i]; j < R[i+1]; j++){
					if(C[j] > i && index_[C[j]] == N){
						_ULonglong key = GetKey(i,C[j]);
						if(!hashmap_->Insert(key)) {
							cout<<"hashtable insert failure!";
						}
					}
				}
			}
		}
	}

	~HybridMap(){
		if (hashmap_ != nullptr)
			delete hashmap_;
		hashmap_ = nullptr;
		if (bitmap_ != nullptr)
			delete bitmap_;
		bitmap_ = nullptr;
	}

	bool Set(u_int x, u_int y){
		if(index_[x] == num_vertices && index_[y] == num_vertices){
			_ULonglong key = GetKey(x, y);
			if(!hashmap_->Insert(key)) {
				cout<<"hashtable insert failure!";
				return false;
			}
			return true;
		}
		bool success = true;
		if(index_[x] != num_vertices){
			if(!bitmap_->bitmapSet(index_[x]*num_vertices + y))
				success = false;
		}
		if(index_[y] != num_vertices){
			if(!bitmap_->bitmapSet(index_[y]*num_vertices + x))
				success = false;
		}
		return success;
	}
	
	bool Get(u_int x, u_int y){
		if(index_[x] == num_vertices && index_[y] == num_vertices){
			_ULonglong key = GetKey(x, y);
			return (hashmap_->exist(key));
		}
		if(index_[x] != num_vertices){
			return (bitmap_->bitmapGet(index_[x]*num_vertices + y));
		}
		if(index_[y] != num_vertices){
			return (bitmap_->bitmapGet(index_[y]*num_vertices + x));
		}
	}

	bool Del(u_int x, u_int y) {
		if(index_[x] == num_vertices && index_[y] == num_vertices){
			_ULonglong key = GetKey(x, y);
			return (hashmap_->Remove(key));
		}
		bool success = true;
		if(index_[x] != num_vertices){
			success = (bitmap_->bitmapDel(index_[x]*num_vertices + y));
		}
		if(index_[y] != num_vertices){
			success = (bitmap_->bitmapDel(index_[y]*num_vertices + x));
		}
		return success;
	}
};



#endif