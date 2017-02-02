#ifndef HISTOGRAM_H
#define HISTOGRAM_H

#include <thrust/device_vector.h>
#include <thrust/host_vector.h>
#include <thrust/sort.h>
#include <thrust/copy.h>
#include <thrust/random.h>
#include <thrust/inner_product.h>
#include <thrust/binary_search.h>
#include <thrust/adjacent_difference.h>
#include <thrust/iterator/constant_iterator.h>
#include <thrust/iterator/counting_iterator.h>
#include <iostream>
#include <iomanip>
#include <iterator>

typedef unsigned int uint__t;
// simple routine to print contents of a vector
template<class T>
struct custom_minus{
	T a;
    
	custom_minus(T _a){
        a=_a;
    }

    __host__ __device__ T operator()(T &x) const{
        return x - a;
    }
};

template<class T, class VecType>
struct ax_functor{
	T a;
    
	ax_functor(T _a){
        a=_a;
    }

    __host__ __device__ T operator()(const VecType &x) const{
        return  a*x;
    }
};

template <typename Vector>
void print_vector(const std::string& name, const Vector& v)
{
  typedef typename Vector::value_type T;
  std::cout << "  " << std::setw(20) << name << "  ";
  thrust::copy(v.begin(), v.end(), std::ostream_iterator<T>(std::cout, " "));
  std::cout << std::endl;
}

//slow version x = scalar * x
template <typename Vector>
void sax_slow(float A, Vector& X) 
{   
	typedef typename Vector::value_type IndexType;
	//thrust::device_vector<IndexType> temp(X.size()); 
	// temp <- A 
	//thrust::fill(temp.begin(), temp.end(), A); 
//	print_vector("temp",temp);
	// temp <- A * X 
	thrust::transform(X.begin(), X.end(), X.begin(), ax_functor<float, IndexType>(A)); 
//	print_vector("x",	X);
}

template <typename Vector>
void substract(Vector& X, long long a)
{
	typedef typename Vector::value_type IndexType;
	custom_minus<long long> f(a);
	//thrust::device_vector<IndexType> temp(X.size()); 
	//thrust::fill(temp.begin(), temp.end(), a); //flag~
	//print_vector("temp",temp);
	//print_vector("prehistogram",X);
	thrust::transform(X.begin(), X.end(), X.begin(), f);
	//thrust::transform(X.begin(), X.end(),), temp.begin(),  thrust::minus<IndexType>());
//	print_vector("tempafter",temp);
	//thrust::copy(temp.begin(), temp.end(), X.begin()); //
  // print_vector("histogramafter",X);
	//thrust::adjacent_difference(X.begin(), X.end(), X.begin());
	//print_vector("difference_histogram",X);
}

// dense histogram using binary search
template <typename Ptr1, 
          typename Vector2>
void dense_histogram(const Ptr1& input,unsigned long long Number, float width, uint__t Num_Bins,
                           Vector2& histogram)
{
  typedef typename Ptr1::value_type ValueType; // input value type
  typedef typename Vector2::value_type IndexType; // histogram index type
  thrust::device_vector<IndexType> temphisto(Num_Bins,0);
  
  thrust::sort(input, input + Number);
  
  thrust::counting_iterator<IndexType> search_begin(1);
  thrust::device_vector<ValueType> absciss(search_begin,search_begin + Num_Bins);
  sax_slow(width, absciss);
  thrust::upper_bound(input, input + Number,
	                  absciss.begin(), absciss.end(),
                      temphisto.begin());
  
  thrust::transform(temphisto.begin(), temphisto.end(), histogram.begin(), histogram.begin(), thrust::plus<IndexType>());
   //thrust::adjacent_difference(histogram.begin(), histogram.end(), histogram.begin());

  thrust::device_vector<IndexType>().swap(temphisto);
  thrust::device_vector<ValueType>().swap(absciss);
}

// dense histogram using binary search and result * 2
template <typename Ptr1, 
          typename Vector2>
void dense_histogram2(const Ptr1& input,unsigned long long Number, float width, uint__t Num_Bins,
                           Vector2& histogram)
{
  typedef typename Ptr1::value_type ValueType; // input value type
  typedef typename Vector2::value_type IndexType; // histogram index type 
  thrust::device_vector<IndexType> temphisto(Num_Bins,0);
  //1.sort
  thrust::sort(input, input + Number);
  //2.generate intervals and upperBound it
  thrust::counting_iterator<ValueType> search_begin(1);
  thrust::device_vector<ValueType> absciss(search_begin,search_begin + Num_Bins);
  sax_slow(width, absciss);
  thrust::upper_bound(input, input + Number, absciss.begin(), absciss.end(),temphisto.begin()); 
  //3. multiply temphisto by 2
  sax_slow(2, temphisto);
  //4.plus previous histogram
  thrust::transform(temphisto.begin(), temphisto.end(), histogram.begin(), histogram.begin(), thrust::plus<IndexType>());
  //thrust::adjacent_difference(temphisto.begin(), temphisto.end(), temphisto.begin()); 
  thrust::device_vector<IndexType>().swap(temphisto);
  thrust::device_vector<ValueType>().swap(absciss);
}

#endif