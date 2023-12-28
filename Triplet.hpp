#ifndef TRIPLET_HPP
#define TRIPLET_HPP

#include <complex>
#include <vector>
#include <memory>


class Triplet
{
    public:

        Triplet() = delete;

        Triplet(int i_, int j_, double val_) : i(i_), j(j_), val(val_) {}

        virtual ~Triplet() {}

        bool operator< (const Triplet& t) const
        {
            if (i < t.i)
            {
                return true;
            }
            else if(i == t.i)
            {
                if (j < t.j)
                {
                    return true;
                }
            }
            return false;
        }

        bool operator== (const Triplet& t) const
        {
            if(i == t.i && j == t.j)
            {
                return true;
            }
            return false;
        }

        int i;
        int j;
        double val;
        
};


#endif // TRIPLET_HPP