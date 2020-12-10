#ifndef SJTU_MATRIX_HPP
#define SJTU_MATRIX_HPP

#include <cstddef>
#include <initializer_list>
#include <utility>
#include <exception>
#include <iterator>

using std::size_t;

namespace sjtu
{
	
	template <class T>
	class Matrix
	{
	private:
		size_t N = 0;
		size_t M = 0;
		size_t siz = 0;
		T *matrix = nullptr;
		
	public:
		template<class T_,class U>
		friend auto operator*(const Matrix<T_> &mat, const U &x);
		template<class T_,class U>
		friend auto operator*(const U &x, const Matrix<T_> &mat);
		template<class U, class V>
		friend auto operator*(const Matrix<U> &a, const Matrix<V> &b);
		template<class U, class V>
		friend auto operator+(const Matrix<U> &a, const Matrix<V> &b);
		template<class U, class V>
		friend auto operator-(const Matrix<U> &a, const Matrix<V> &b);
		Matrix() = default;
		
		Matrix(size_t n, size_t m, T _init = T())
		{
			N = n;
			M = m;
			siz = n * m;
			matrix = new T[siz];
			for(int i=0; i<siz; i++) {
				matrix[i] = _init;
			}
		}
		
		explicit Matrix(std::pair<size_t, size_t> sz, T _init = T())
		{
			N = sz.first;
			M = sz.second;
			matrix = new T[siz = N * M];
			for(int i=0; i<siz; i++) {
				matrix[i] = _init;
			}
		}
		
		Matrix(const Matrix &o)
		{
			N = o.rowLength();
			M = o.columnLength();
			matrix = new T[siz = N * M];
			for (int i=0;i<siz;i++) {
				matrix[i] = o.matrix[i];
			} 
		}
		
		template <class U>
		Matrix(const Matrix<U> &o)
		{
			N = o.rowLength();
			M = o.columnLength();
			siz = N * M;
			matrix = new T[siz];
			for(int i=0;i<N;i++){
				for(int j=0;j<M;j++){
					matrix[i * M + j] = static_cast<T>(o(i, j));
				}
			}
		}
		
		Matrix &operator=(const Matrix &o)
		{
			if(this == &o) return *this;
			delete [] matrix;
			siz = o.totalLength();
			N = o.rowLength();
			M = o.columnLength();
			matrix = new T[siz];
			for(int i=0;i<N;i++){
				for(int j=0;j<M;j++){
					matrix[i*M+j] = o(i, j);
				}
			}
			return *this;
		}
		
		template <class U>
		Matrix &operator=(const Matrix<U> &o)
		{
			delete [] matrix;
			siz = o.totalLength();
			N = o.rowLength();
			M = o.columnLength();
			matrix = new T[siz];
			for(int i=0;i<N;i++){
				for(int j=0;j<M;j++){
					matrix[i*M+j] = static_cast<T>(o(i, j));
				}
			}
			return *this;
		}
		
		Matrix(Matrix &&o) noexcept
		{
			matrix = o.matrix;
			N = o.rowLength();
			M = o.columnLength();
			siz = o.totalLength();
			o.matrix = nullptr;
		}
		
		Matrix &operator=(Matrix &&o) noexcept
		{
			if(this == &o) return *this;
			delete [] matrix;
			matrix = o.matrix;
			N = o.rowLength();
			M = o.columnLength();
			siz = o.totalLength();
			o.matrix = nullptr;
			return *this;
		}
		
		~Matrix() 
		{
			if(matrix != nullptr){
				delete [] matrix;
			}
		}
		
		Matrix(std::initializer_list<std::initializer_list<T>> il)
		{
			size_t num = 0;
			N = il.size();
			M = (*std::begin(il)).size();
			siz = N * M;
			matrix = new T[siz];
			for(auto i0=std::begin(il);i0<std::end(il);i0++){
				if(i0->size() != M){
					delete [] matrix;
					throw std::invalid_argument("init_list");
				}
				for(auto j0=std::begin(*i0);j0<std::end(*i0);j0++){
					matrix[num] = *j0;
					num += 1;
				}
			}
		}
		
	public:
		size_t rowLength() const {return N;}
		
		size_t columnLength() const {return M;}
		
		size_t totalLength() const {return siz;}
		
		void resize(size_t _n, size_t _m, T _init = T())
		{
			size_t newsiz = _n * _m;
			if(newsiz == siz){
				N = _n;
				M = _m;
			}
			else{
				T *temp = new T[newsiz];
				for(int i=0;i<newsiz;i++){
					if(i < siz) temp[i] = matrix[i];
					else temp[i] = _init;
				}
				T *tmpmat = matrix;
				matrix = temp;
				delete [] tmpmat;
				N = _n;
				M = _m;
				siz = newsiz;
			}
		}
		
		void resize(std::pair<size_t, size_t> sz, T _init = T())
		{
			size_t n = sz.first;
			size_t m = sz.second;
			size_t newsiz = n * m;
			if(newsiz == siz){
				N = n;
				M = m;
			}
			else{
				T *temp = new T[newsiz];
				for(int i=0;i<newsiz;i++){
					if(i < siz) temp[i] = matrix[i];
					else temp[i] = _init;
				}
				T *tmpmat = matrix;
				matrix = temp;
				delete [] tmpmat;
				N = n;
				M = m;
				siz = newsiz;
			}
		}
		
		std::pair<size_t, size_t> size() const
		{
			std::pair<size_t, size_t> mat;
			mat.first = N;
			mat.second = M;
			return mat;
		};
		
		void clear()
		{
			if(matrix != nullptr)
			delete [] matrix;
			matrix = nullptr;
			N = 0;
			M = 0;
			siz = 0;
		}
		
	public:
		const T &operator()(size_t i, size_t j) const
		{
			if(i<0 || j<0 || i>=N || j>=M) throw std::invalid_argument("()");
			return matrix[i * M + j];
		}
		
		T &operator()(size_t i, size_t j)
		{
			if(i<0 || j<0 || i>=N || j>=M) throw std::invalid_argument("()");
			return matrix[i * M + j];
		}
		
		Matrix<T> row(size_t i) const
		{
			if(i<0 || i>=N) throw std::invalid_argument("row");
			Matrix<T> newMat;
			newMat.N = 1;
			newMat.M = M;
			newMat.siz = M;
			newMat.matrix = new T[M];
			for(int k=0;k<newMat.siz;k++){
				newMat.matrix[k] = matrix[i*M+k];
			} 
			return newMat;
		}
		
		Matrix<T> column(size_t i) const
		{
			if(i<0 || i>=M) throw std::invalid_argument("column");
			Matrix<T> newMat;
			newMat.N = N;
			newMat.M = 1;
			newMat.siz = N;
			newMat.matrix = new T[N];
			for(int k=0;k<newMat.siz;k++){
				newMat.matrix[k] = matrix[k*M+i];
			} 
			return newMat;
		}
		
		
	public:
		template <class U>
		bool operator==(const Matrix<U> &o) const
		{
			if(M != o.columnLength() || N != o.rowLength()) return false;
			for(int i=0;i<N;i++){
				for(int j=0;j<M;j++){
					if(matrix[i*M+j] != o(i, j))return false;
				}
			}
			return true;
		}
		
		template <class U>
		bool operator!=(const Matrix<U> &o) const
		{
			if(M != o.columnLength() || N != o.rowLength()) return true;
			for(int i=0;i<N;i++){
				for(int j=0;j<M;j++){
					if(matrix[i*M+j] != o(i, j))return true;
				}
			}
			return false;
		}
		
		Matrix operator-() const
		{
			Matrix<T> newMat;
			newMat.N = N;
			newMat.M = M;
			newMat.siz = siz;
			newMat.matrix = new T[siz];
			for(int i=0;i<siz;i++){
				newMat.matrix[i] = -matrix[i];
			}
			return newMat;
		}
		
		template <class U>
		Matrix &operator+=(const Matrix<U> &o)
		{
			if(N != o.rowLength() || M != o.columnLength()) throw std::invalid_argument("+="); 
			for(int i=0;i<N;i++){
				for(int j=0;j<M;j++){
					matrix[i*M+j] += o(i, j);
				}
			}
			return *this;
		}
		
		template <class U>
		Matrix &operator-=(const Matrix<U> &o)
		{
			if(N != o.rowLength() || M != o.columnLength()) throw std::invalid_argument("-="); 
			for(int i=0;i<N;i++){
				for(int j=0;j<M;j++){
					matrix[i*M+j] -= o(i, j);
				}
			}
			return *this;
		}
		
		template <class U>
		Matrix &operator*=(const U &x)
		{
			for(int i=0;i<siz;i++){
				//matrix[i] = matrix[i] * static_cast<T>(x);
				matrix[i] *= static_cast<T>(x); 
			}
			//*this = *this * x;
			return *this;
		}
		
		Matrix tran() const
		{
			Matrix<T> trMat;
			trMat.N = M;
			trMat.M = N;
			trMat.siz = siz;
			trMat.matrix = new T[siz];
			for(int i=0;i<M;i++){
				for(int j=0;j<N;j++){
					trMat.matrix[i*N+j] = matrix[j*M+i];
				}
			}
			return trMat;
		}
		
	public: // iterator
		class iterator
		{
		public:
			using iterator_category = std::random_access_iterator_tag;
			using value_type        = T;
			using pointer           = T *;
			using reference         = T &;
			using size_type         = size_t;
			using difference_type   = std::ptrdiff_t;
			
			iterator() = default;
			
			iterator(const iterator &) = default;
			
			iterator &operator=(const iterator &) = default;
			
		private:
			size_t now_N,now_M,lst_N,lst_M,nxt_N,nxt_M,mat_N,mat_M,sub_N,sub_M;
			T* p;
			friend class Matrix;
			
		public:
			difference_type operator-(const iterator &o)
			{
				if(o.lst_N != lst_N || o.lst_M != lst_M || o.nxt_M != nxt_M || o.nxt_N != nxt_N) throw std::invalid_argument("operator-");
				return (((now_N - lst_N) * sub_M + (now_M - lst_M))-((o.now_N - o.lst_N) * o.sub_M + (o.now_M - o.lst_M)));
			}
			
			iterator &operator+=(difference_type offset)
			{
				p -= (now_N * mat_M + now_M);
				int s0 = ((now_N - lst_N) * sub_M +(now_M-lst_M)) + offset;
				now_N = s0 / sub_M + lst_N;
				now_M = s0 % sub_M + lst_M;
				p += (now_N * mat_M + now_M);
				return *this;
			}
			
			iterator operator+(difference_type offset) const
			{
				iterator It = *this;
				It += offset;
				return It;
			}
			
			iterator &operator-=(difference_type offset)
			{
				p -= (now_N * mat_M + now_M);
				int s0 = ((now_N - lst_N) * sub_M +(now_M-lst_M)) - offset;
				now_N = s0 / mat_M + lst_N;
				now_M = s0 % mat_M + lst_M;
				p += (now_N * mat_M + now_M);
				return *this;
			}
			
			iterator operator-(difference_type offset) const
			{
				iterator It = *this;
				It -= offset;
				return It;
			}
			
			iterator &operator++()
			{
				*this += 1;
				return *this;
			}
			
			iterator operator++(int)
			{
				iterator It = *this;
				*this += 1;
				return It;
			}
			
			iterator &operator--()
			{
				*this -= 1;
				return *this;
			}
			
			iterator operator--(int)
			{
				iterator It = *this;
				*this -= 1;
				return It;
			}
			
			reference operator*() const
			{
				return *p;
			}
			
			pointer operator->() const
			{
				return p;
			}
			
			bool operator==(const iterator &o) const
			{
				return (p == o.p);
			}
			
			bool operator!=(const iterator &o) const
			{
				return (p != o.p);
			}
		};
		
		iterator begin()
		{
			iterator It;
			It.now_N = It.now_M = It.lst_N = It.now_M = 0;
			It.nxt_N = N-1;
			It.nxt_M = M-1;
			It.sub_N = It.mat_N = N;
			It.sub_M = It.mat_M = M;
			It.p = matrix;
			return It;
		}
		
		iterator end()
		{
			iterator It;
			It.lst_N = It.now_M = 0;
			It.nxt_N = N-1;
			It.nxt_M = M-1;
			It.now_N = N;
			It.now_M = 0;
			It.sub_N = It.mat_N = N;
			It.sub_M = It.mat_M = M;
			It.p = matrix;
			return It;
		}
		
		std::pair<iterator, iterator> subMatrix(std::pair<size_t, size_t> l, std::pair<size_t, size_t> r)
		{
			iterator l0, r0;
			l0.mat_N = r0.mat_N = N;
			l0.mat_M = r0.mat_M = M;
			l0.sub_N = r0.sub_N = r.first - l.first + 1;
			l0.sub_M = r0.sub_M = r.second - l.second + 1;
			l0.now_N = l0.lst_N = r0.lst_N = l.first;
			l0.now_M = l0.lst_M = r0.lst_M = l.second;
			r0.now_N = l0.nxt_N = r0.nxt_N = r.first;
			r0.now_M = l0.nxt_M = r0.nxt_M = r.second;
			l0.p = matrix + l.first * M + l.second;
			r0.p = matrix + r.first * M + r.second;
			return {l0, ++r0};
        }
	};
		
}

//
namespace sjtu
{ 
	
	template <class T, class U>
	auto operator*(const Matrix<T> &mat, const U &x)
	{
		size_t n = mat.N;
		size_t m = mat.M;
		Matrix<decltype(mat.matrix[0]*x)> newMat(n, m);
		for(int i=0;i<n;i++){
			for(int j=0;j<m;j++){
				newMat.matrix[i*m+j]=mat.matrix[i*m+j] * x;
			}
		}
		return newMat;
	}
	
	template <class T, class U>
	auto operator*(const U &x, const Matrix<T> &mat)
	{
		size_t n = mat.N;
		size_t m = mat.M;
		Matrix<decltype(mat.matrix[0]*x)> newMat(n, m);
		for(int i=0;i<n;i++){
			for(int j=0;j<m;j++){
				newMat.matrix[i*m+j]=mat.matrix[i*m+j] * x;
			}
		}
		return newMat;
	}
	
	template <class U, class V>
	auto operator*(const Matrix<U> &a, const Matrix<V> &b)
	{
		if(a.M != b.N) throw std::invalid_argument("Matrix*");
		size_t n = a.N;
		size_t m = b.M;
		Matrix<decltype(a.matrix[0]*b.matrix[0])> newMat(n, m, 0);
		for(int i=0;i<n;i++){
			for(int j=0;j<m;j++){
				for(int k=0;k<a.M;k++){
					newMat.matrix[i*m+j] += a.matrix[i*a.M+k] * b.matrix[k*b.M+j];
				}
			}
		}
		return newMat;
	}
	
	template <class U, class V>
	auto operator+(const Matrix<U> &a, const Matrix<V> &b)
	{
		if(a.N != b.N || a.M != b.M) throw std::invalid_argument("Matrix+");
		size_t n = a.N;
		size_t m = a.M;
		Matrix<decltype(a.matrix[0] + b.matrix[0])> newMat(n, m, 0);
		for(int i=0;i<n;i++){
			for(int j=0;j<m;j++){
				newMat.matrix[i*m+j] = a.matrix[i*m+j] + b.matrix[i*m+j];
			}
		}
		return newMat;
	}
	
	template <class U, class V>
	auto operator-(const Matrix<U> &a, const Matrix<V> &b)
	{
		if(a.N != b.N || a.M != b.M) throw std::invalid_argument("Matrix-");
		size_t n = a.N;
		size_t m = a.M;
		Matrix<decltype(a.matrix[0] - b.matrix[0])> newMat(n, m, 0);
		for(int i=0;i<n;i++){
			for(int j=0;j<m;j++){
				newMat.matrix[i*m+j] = a.matrix[i*m+j] - b.matrix[i*m+j];
			}
		}
		return newMat;
	}
}

#endif //SJTU_MATRIX_HPP
