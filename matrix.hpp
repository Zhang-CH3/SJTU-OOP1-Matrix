#ifndef SJTU_MATRIX_HPP
#define SJTU_MATRIX_HPP

#include <cstddef>
#include <initializer_list>
#include <utility>
#include <iterator>
#include <exception>

using std::size_t;

namespace sjtu {

    template<class T>
    class Matrix {
    private:
        // your private member variables here.
        const std::pair<size_t, std::size_t> init_size={0,0};
        bool temp=false;
        std::pair<size_t, std::size_t> size_Matrix=init_size;
        T *Core = nullptr;
        size_t len = 0,true_len=0;
    public:
        template<class T1,class U>
        friend auto operator*(const Matrix<T1> &mat, const U &x);

        template<class T1,class U>
        friend auto operator*(const U &x, const Matrix<T1> &mat);

        template<class U, class V>
        friend auto operator*(const Matrix<U> &a, const Matrix<V> &b);

        template<class U, class V>
        friend auto operator+(const Matrix<U> &a, const Matrix<V> &b);

        template<class U, class V>
        friend auto operator-(const Matrix<U> &a, const Matrix<V> &b);

        Matrix() = default;

        Matrix(size_t n, size_t m, T _init = T()) {
            size_Matrix.first = n;
            size_Matrix.second = m;
            true_len=len = n * m;
            Core = new T[len];
            for (int i = 0; i < len; ++i) Core[i] = _init;
        }

        explicit Matrix(std::pair<size_t, size_t> sz, T _init = T()) {
            size_Matrix = sz;
            true_len=len = sz.first * sz.second;
            Core = new T[len];
            for (int i = 0; i < len; ++i) Core[i] = _init;
        }

        Matrix(const Matrix &o) {
            size_Matrix = o.size_Matrix;
            true_len=len = o.len;
            Core = new T[len];
            for (int i = 0; i < len; ++i) Core[i] = o.Core[i];
        }

        template<class U>
        Matrix(const Matrix<U> &o) {
            std::pair<size_t, std::size_t> sz=o.size();
            size_Matrix = sz;
            true_len=len = sz.second*sz.first;
            Core = new T[len];
            for (int i = 0; i < sz.first; ++i)
                for (int j=0;j<sz.second;++j)
                    Core[i*sz.second+j] =static_cast<T>(o(i,j));
        }

        Matrix &operator=(const Matrix &o) {
            if (this==&o) return *this;
            delete[] Core;
            true_len=len = o.len;
            size_Matrix = o.size_Matrix;
            Core = new T[len];
            for (int i = 0; i < len; ++i) Core[i] = o.Core[i];
            return *this;
        }

        template<class U>
        Matrix &operator=(const Matrix<U> &o) {
            delete [] Core;
            std::pair<size_t, std::size_t> sz=o.size();
            true_len=len=sz.first*sz.second;
            size_Matrix = sz;
            Core=new T[len];
            for (int i=0;i<sz.first;++i)
                for (int j=0;j<sz.second;++j)
                    Core[i*sz.second+j]=static_cast<T>(o(i,j));
            return *this;
        }

        Matrix(Matrix &&o) noexcept {
            Core = o.Core;
            size_Matrix = o.size_Matrix;
            true_len=len = o.len;
            o.Core = nullptr;
        }

        Matrix &operator=(Matrix &&o) noexcept {
            if (this==&o) return *this;
            delete [] Core;
            Core = o.Core;
            size_Matrix = o.size_Matrix;
            true_len=len = o.len;
            o.Core = nullptr;
            return *this;
        }

        ~Matrix() {
            if (Core != nullptr)
                delete[] Core;
        }

        Matrix(std::initializer_list<std::initializer_list<T>> il) {
            size_Matrix.first=il.size();
            size_Matrix.second=(*std::begin(il)).size();
            true_len=len=size_Matrix.first*size_Matrix.second;
            Core=new T [len];
            int i=-1;
            for (auto it_i=std::begin(il);it_i!=std::end(il);++it_i){
                if (it_i->size()!=size_Matrix.second) {
                    delete [] Core;
                    throw std::invalid_argument("init_list");
                }
                for (auto it_j=std::begin(*it_i);it_j!=std::end(*it_i);++it_j){
                    Core[++i]=*it_j;
                }
            }
        }

    public:
        size_t rowLength() const {
            return size_Matrix.first;
        }

        size_t columnLength() const {
            return size_Matrix.second;
        }

        void resize(size_t _n, size_t _m, T _init = T()) {
            if (_n*_m==len){
                size_Matrix={_n,_m};
            }else {
                T* tmp=new T [_n*_m];
                for (int i=0;i<_n*_m;++i)
                    if (i<len) tmp[i]=Core[i];
                    else tmp[i]=_init;
                T* _tmp=Core;
                Core=tmp;
                delete [] _tmp;
                size_Matrix={_n,_m};
                true_len=len=_n*_m;
            }
        }

        void resize(std::pair<size_t, size_t> sz, T _init = T()) {
            std::size_t _n=sz.first,_m=sz.second;
            if (_n*_m==len){
                size_Matrix={_n,_m};
            }else {
                T* tmp=new T [_n*_m];
                for (int i=0;i<_n*_m;++i)
                    if (i<len) tmp[i]=Core[i];
                    else tmp[i]=_init;
                T* _tmp=Core;
                Core=tmp;
                delete [] _tmp;
                size_Matrix={_n,_m};
                true_len=len=_n*_m;
            }
        }

        std::pair<size_t, size_t> size() const {
            return size_Matrix;
        };

        void clear() {
            if (Core != nullptr)
            delete[] Core;
            Core = nullptr;
            true_len=len = 0;
            size_Matrix=init_size;
        }

    public:
        const T &operator()(size_t i, size_t j) const {
            if (i<0 || j<0 || i>=size_Matrix.first || j>=size_Matrix.second) throw std::invalid_argument("()");
            return Core[i*size_Matrix.second+j];
        }

        T &operator()(size_t i, size_t j) {
            if (i<0 || j<0 || i>=size_Matrix.first || j>=size_Matrix.second) throw std::invalid_argument("()");
            return Core[i*size_Matrix.second+j];
        }

        Matrix<T> row(size_t i) const {
            if (i<0||i>=size_Matrix.first) throw std::invalid_argument("row");
            Matrix<T> res;
            res.size_Matrix={1,size_Matrix.second};
            res.true_len=res.len=size_Matrix.second;
            res.Core=new T [res.len];
            for (int j=0;j<res.len;++j) res.Core[j]=Core[i*size_Matrix.second+j];
            return res;
        }

        Matrix<T> column(size_t i) const {
            if (i<0||i>=size_Matrix.second) throw std::invalid_argument("column");
            Matrix<T> res;
            res.size_Matrix={size_Matrix.first,1};
            res.true_len=res.len=size_Matrix.first;
            res.Core=new T [res.len];
            for (int j=0;j<res.len;++j) res.Core[j]=Core[j*size_Matrix.second+i];
            return res;
        }


    public:
        template<class U>
        bool operator==(const Matrix<U> &o) const {
            if (size_Matrix!=o.size()) return 0;
            for (int i=0;i<size_Matrix.first;++i)
                for (int j=0;j<size_Matrix.second;++j)
                    if (Core[i*size_Matrix.second+j]!=o(i,j)) return 0;
            return 1;
        }

        template<class U>
        bool operator!=(const Matrix<U> &o) const {
            if (size_Matrix!=o.size()) return 1;
            for (int i=0;i<size_Matrix.first;++i)
                for (int j=0;j<size_Matrix.second;++j)
                if (Core[i*size_Matrix.second+j]!=o(i,j)) return 1;
            return 0;
        }

        Matrix operator-() const {
            Matrix<T> res;
            res.size_Matrix=size_Matrix;
            res.true_len=res.len=len;
            res.Core=new T [len];
            for (int i=0;i<len;++i) res.Core[i]=-Core[i];
            return res;
        }

        template<class U>
        Matrix &operator+=(const Matrix<U> &o) {
            if (size_Matrix!=o.size()) throw std::invalid_argument("+=");
            for (int i=0;i<size_Matrix.first;++i)
                for (int j=0;j<size_Matrix.second;++j)
                    Core[i*size_Matrix.second+j]+=o(i,j);
//            for (int i=0;i<len;++i) Core[i]+=static_cast<T>(o.Core[i]);
            return *this;
        }

        template<class U>
        Matrix &operator-=(const Matrix<U> &o) {
            if (size_Matrix!=o.size()) throw std::invalid_argument("-=");
            for (int i=0;i<size_Matrix.first;++i)
                for (int j=0;j<size_Matrix.second;++j)
                    Core[i*size_Matrix.second+j]-=o(i,j);
//            for (int i=0;i<len;++i) Core[i]-=static_cast<T>(o.Core[i]);
            return *this;
        }

        template<class U>
        Matrix &operator*=(const U &x) {
            for (int i=0;i<len;++i) Core[i]*=static_cast<T>(x);
            return *this;
        }

        Matrix tran() const {
            Matrix<T> res;
            res.size_Matrix={size_Matrix.second,size_Matrix.first};
            res.true_len=res.len=len;
            res.Core=new T [len];
            for (int i=0;i<size_Matrix.second;++i)
                for (int j=0;j<size_Matrix.first;++j)
                    res.Core[i*size_Matrix.first+j]=Core[j*size_Matrix.second+i];
            return res;
        }

    public: // iterator
        class iterator {
        public:
            using iterator_category = std::random_access_iterator_tag;
            using value_type = T;
            using pointer = T *;
            using reference = T &;
            using size_type = size_t;
            using difference_type = std::ptrdiff_t;

            iterator() = default;

            iterator(const iterator &) = default;

            iterator &operator=(const iterator &) = default;

        private:
            T* ptr;
            std::pair<size_type,size_type> now,l_up,r_down,size_mat,size_sub;
            friend class Matrix;
        public:
            difference_type operator-(const iterator &o) {
                if (o.l_up!=l_up || o.r_down!=r_down) throw std::invalid_argument("operator-");
                return (((now.first-l_up.first)*size_sub.second+(now.second-l_up.second))-((o.now.first-o.l_up.first)*o.size_sub.second+(o.now.second-o.l_up.second)));
            }
            iterator &operator+=(difference_type offset) {
                ptr-=(now.first*size_mat.second+now.second);
                int stat=((now.first-l_up.first)*size_sub.second+(now.second-l_up.second))+offset;
                now.first=l_up.first+stat/size_sub.second;
                now.second=l_up.second+stat%size_sub.second;
                ptr+=(now.first*size_mat.second+now.second);
                return *this;
            }

            iterator operator+(difference_type offset) const {
                iterator res=*this;
                res+=offset;
                return res;
            }

            iterator &operator-=(difference_type offset) {
                ptr-=(now.first*size_mat.second+now.second);
                int stat=((now.first-l_up.first)*size_sub.second+(now.second-l_up.second))-offset;
                now.first=l_up.first+stat/size_mat.second;
                now.second=l_up.second+stat%size_mat.second;
                ptr+=(now.first*size_mat.second+now.second);
                return *this;
            }

            iterator operator-(difference_type offset) const {
                iterator res=*this;
                res-=offset;
                return res;
            }

            iterator &operator++() {
                *this+=1;
                return *this;
            }

            iterator operator++(int) {
                iterator res=*this;
                *this+=1;
                return res;
            }

            iterator &operator--() {
                *this-=1;
                return *this;
            }

            iterator operator--(int) {
                iterator res=*this;
                *this-=1;
                return res;
            }


            reference operator*() const {
                return *ptr;
            }

            pointer operator->() const {
                return ptr;
            }

            bool operator==(const iterator &o) const {
                return (ptr==o.ptr);
            }

            bool operator!=(const iterator &o) const {
                return (ptr!=o.ptr);
            }
        };

        iterator begin() {
            iterator res;
            res.now=res.l_up={0,0};
            res.r_down={size_Matrix.first-1,size_Matrix.second-1};
            res.size_sub=res.size_mat=size_Matrix;
            res.ptr=Core;
            return res;
        }

        iterator end() {
            iterator res;
            res.l_up={0,0};
            res.r_down={size_Matrix.first-1,size_Matrix.second-1};
            res.now={size_Matrix.first,0};
            res.size_sub=res.size_mat=size_Matrix;
            res.ptr=Core+len;
            return res;
        }

        std::pair<iterator, iterator> subMatrix(std::pair<size_t, size_t> l, std::pair<size_t, size_t> r) {
            iterator ll,rr;
            ll.size_mat=rr.size_mat=size_Matrix;
            ll.size_sub=rr.size_sub={r.first-l.first+1,r.second-l.second+1};
            ll.l_up=rr.l_up=l;
            rr.r_down=rr.r_down=r;
            ll.now=l;
            rr.now=r;
            ll.ptr=Core+l.first*size_Matrix.second+l.second;
            rr.ptr=Core+r.first*size_Matrix.second+r.second;
            return {ll,++rr};
        }
    };

}

//
namespace sjtu {
    template<class T, class U>
    auto operator*(const Matrix<T> &mat, const U &x) {
        std::pair<size_t,size_t> size_=mat.size_Matrix;
        Matrix<decltype(mat.Core[0]*x)> res(size_);
        for (int i=0;i<size_.first;++i)
            for (int j=0;j<size_.second;++j)
                res.Core[i*size_.second+j]=mat.Core[i*size_.second+j]*x;
        return res;
    }

    template<class T, class U>
    auto operator*(const U &x, const Matrix<T> &mat) {
        std::pair<size_t,size_t> size_=mat.size_Matrix;
        Matrix<decltype(mat.Core[0]*x)> res(size_);
        for (int i=0;i<size_.first;++i)
            for (int j=0;j<size_.second;++j)
                res.Core[i*size_.second+j]=mat.Core[i*size_.second+j]*x;
        return res;
    }

    template<class U, class V>
    auto operator*(const Matrix<U> &a, const Matrix<V> &b) {
        if (a.size_Matrix.second!=b.size_Matrix.first) throw std::invalid_argument("Matrix*");
        std::pair<size_t,size_t> size_={a.size_Matrix.first,b.size_Matrix.second};
        Matrix<decltype(a.Core[0]*b.Core[0])> res(size_,0);
        for (int i=0;i<a.size_Matrix.first;++i)
            for (int j=0;j<b.size_Matrix.second;++j)
                for (int k=0;k<a.size_Matrix.second;++k)
                    res.Core[i*size_.second+j]+=a.Core[i*a.size_Matrix.second+k]*b.Core[k*b.size_Matrix.second+j];
        return res;
    }

    template<class U, class V>
    auto operator+(const Matrix<U> &a, const Matrix<V> &b) {
        if (a.size()!=b.size()) throw std::invalid_argument("Matrix+");
        std::pair<size_t,size_t> size_={a.size_Matrix.first,a.size_Matrix.second};
        Matrix<decltype(a.Core[0]+b.Core[0])> res(size_,0);
        for (int i=0;i<a.size_Matrix.first;++i)
            for (int j=0;j<a.size_Matrix.second;++j)
                    res.Core[i*size_.second+j]=a.Core[i*size_.second+j]+b.Core[i*size_.second+j];
        return res;
    }

    template<class U, class V>
    auto operator-(const Matrix<U> &a, const Matrix<V> &b) {
        if (a.size()!=b.size()) throw std::invalid_argument("Matrix-");
        std::pair<size_t,size_t> size_={a.size_Matrix.first,a.size_Matrix.second};
        Matrix<decltype(a.Core[0]-b.Core[0])> res(size_,0);
        for (int i=0;i<a.size_Matrix.first;++i)
            for (int j=0;j<a.size_Matrix.second;++j)
                res.Core[i*size_.second+j]=a.Core[i*size_.second+j]-b.Core[i*size_.second+j];
        return res;
    }
}

#endif //SJTU_MATRIX_HPP
