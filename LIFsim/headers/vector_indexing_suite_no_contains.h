/* 
 * File:   vector_indexing_suite_no_contains.h
 * Author: ranr
 *
 * Created on July 5, 2011, 11:41 AM
 */

#ifndef VECTOR_INDEXING_SUITE_NO_CONTAINS_H
#define	VECTOR_INDEXING_SUITE_NO_CONTAINS_H
#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>

using namespace boost::python;

template < class Container >
class vector_indexing_suite_no_contains : public vector_indexing_suite<Container,false,vector_indexing_suite_no_contains<Container> >{
public:
    typedef typename Container::value_type key_type;
    
    static bool
        contains(Container& container, key_type const& key)
        {
                std::cout<<"ERROR: No contains implementation for this class returning False\n\n"<<std::endl;
                return false;
        }
};



#endif	/* VECTOR_INDEXING_SUITE_NO_CONTAINS_H */

