// absorbed into simd.hpp

template <class X>
    class simd_vector
    {
    public:

        typedef simd_traits<X>::value_type value_type;

        // ...

        value_type operator[](size_t index) const
        {
            size_t size = simd_traits<X>::size;
            value_type v[size];
            (*this)().store_u(v);
            return v[index];
        }
    };


