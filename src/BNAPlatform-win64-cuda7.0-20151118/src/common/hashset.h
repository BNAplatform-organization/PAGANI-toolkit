#ifndef  _HASH_SET_H_
#define  _HASH_SET_H_ 
#include <stdlib.h>
#include <string>
#include <iostream>
#include <vector>
using namespace std;

const _ULonglong a = 11400714819323198485;

//enum Status//设置状态的数组
//{
//    EMPTY,
//	EXIST,        
//};

//template<class K, class V>
//struct KeyValue//字典
//{
//	    K _key;
//		V _value;
//		KeyValue(const K& key = K(), const V& value = V())
//		//设置K()和V()为了无参初始化
//	    :_key(key), _value(value) {}
//};

//static size_t BKDRHash(const char * str)//字符串哈希算法
//{
//    unsigned int seed = 131; // 31 131 1313 13131 131313
//    unsigned int hash = 0;
//    while (*str)
//    {
//        hash = hash * seed + (unsigned int)(*str++);
//    }
//    return (hash & 0x7FFFFFFF);
//}

_ULonglong GetKey(_ULonglong v1, _ULonglong v2)
{
	return v1<v2 ? ((v1<<32) | v2) : ((v2<<32) | v1);
}

static size_t FibonacciHash(const _ULonglong & x)
{
	return a*x;
}

//默认哈希函数实现
template<class K>
struct DefaultHashFuncer//基本类型
{
    size_t operator()(const K& key)
    {
        return key;
    }
};

template<>
struct DefaultHashFuncer<_ULonglong>      //string类型--模板的特化
{
    size_t operator()(const _ULonglong & x)
    {
        return FibonacciHash(x);
    }
};



template<class K, class HashFuncer = DefaultHashFuncer<K>>
class HashSet
{
    protected:
        size_t _size;//哈希表中哈希数的个数
       	size_t _capacity;//哈希表的大小
		K* _table;//存放哈希数
        //Status* _status;//存放状态的数组
    
	public://进行各函数的实现--进行增删查改
		HashSet():_table(NULL), _size(0), _capacity(0), total_count(0),find_cnt(0) {}
		HashSet(size_t size, double loadfactor_inv = 2.0): 
			_size(0), 
			_capacity((size_t)(loadfactor_inv*size)),
			_table(new K[_capacity]),
			total_count(0),
			find_cnt(0)
		{
			memset(_table, 0, sizeof(K)*_capacity);
		}
		
		~HashSet()
		{
		    if (_table)
		    {
		        delete[] _table;
		        //delete[] _status;
		        _size = 0;
		        _capacity = 0;
		    }
		}
		
		size_t HashFunc(const K& key)             //求出key在哈希表中的位置
		{
    		HashFuncer hp;
    		return hp(key)%_capacity;//hp(key)调用仿函数
		}
		
		/*size_t HashFunc0(const K& key)
		{
    			return HashFunc(key);//调用HashFunc函数，找到二次探测最初的位置
		}		
		
		size_t HashFunci(size_t index, size_t i)
		{
    			return index + (2 * i - 1);//优化后的算法
		}*/
		
		size_t Find(const K& key)
		{
			//size_t cnt = 0;
			find_cnt++;
			size_t index = HashFunc(key);
			while (_table[index] != 0)
			{
        		if ( _table[index] == key )
        		{
            		return index;
        		}
				++index;
				total_count++;
				if (index == _capacity)//如果哈希到了数组最后一位，就返回到第一位进行哈希
    			{
    			    index = 0;
    			}
			}
			return _capacity;
		}		
		
		bool exist(const K& key){
			return (Find(key) != _capacity);
		}

		bool Insert(const K& key) //防止冗余用bool
		{
    		//CheckCapacity(_size + 1);//检查容量，不足增容
			//线性探测
			find_cnt++;
			size_t index = HashFunc(key);
			while (_table[index] != 0 && _table[index] != _ULLONG_MAX)//如果不为EMPTY或DELETE就不能在此位置处存入key，存在哈希冲突
			{
				if (_table[index] == key)//如果key已存在，就插入失败
    			{
    			    return false;
    			}
    			++index;
				++total_count;
    			if (index == _capacity)//如果哈希到了数组最后一位，就返回到第一位进行哈希
    			{
    			    index = 0;
    			}
    		}
    		++_size;
    		_table[index] = key;
    		//_table[index]._value = value;
    		//_status[index] = EXIST;
    		return true;
    		//二次探测
    		//size_t i = 0;
			//size_t index = HashFunc0(key);
			//while (_status[index] == EXIST)//如果不为EMPTY或DELETE就不能在此位置处存入key，存在哈希冲突
    		//{
        	//	if (_table[index]._key == key && _table[index]._value == value)//如果key已存在，就插入失败
        	//	{
            //			return false;
        	//	}
        	//	index = HashFunci(index, ++i);
        	//	if (index >= _capacity)//如果哈希到的位置超过了数组最后一位，就从首位开始求出对应位置
        	//	{
            //			index = index - _capacity;
        	//	}
    		//}
    		//_table[index]._key = key;
    		//_table[index]._value = value;
    		//_status[index] = EXIST;
    		//_size++;
    		//return true;;
		}	

		bool Remove(const K& key)
		{
    			size_t pos = Find(key);
				if (pos == _capacity)//返回_capacity表示查找失败
    			{
        			return false;
    			}
				//_status[pos] = EMPTY;
				_table[pos] = _ULLONG_MAX;
				--_size;
    			return true;
		}
		
		size_t size(){
			return _size;
		}

		void GetValVector(vector<K> &vec)
		{
			//vector<K> vec;
			//vec.reserve(_size+1);
			for(size_t i = 0; i < _capacity; i++)
			{
				if(_table[i] != 0 && _table[i] != _ULLONG_MAX){
					vec.push_back(_table[i]);	
					//cout<<_table[i]._value<<endl;
				}
			}
			//return vec;
			return;
		}		

		size_t total_count;
		size_t find_cnt;

};


#endif
