 #ifndef BOOL_ITERATOR
 #define BOOL_ITERATOR

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <climits>
#include <iostream>
#include <assert.h>
#include <thread>
#include <vector>

#include <iterator>

using namespace std;

class bool_iterator : public std::iterator<std::forward_iterator_tag, bool>{

public:

	uint64_t curr;
	uint64_t size;
    vector<bool>* vect;
    bool initialized = false;
    bool ended = false;
	
	bool_iterator() {
	}



	bool_iterator(bool b){
		ended = b;
		initialized = true;
	}
	


	bool_iterator(const bool_iterator& cr){
		curr = cr.curr;
		vect = cr.vect;
	}
	


	bool_iterator(vector<bool>* vect){
		curr = vect->size();
		vect = vect;
		size = vect->size();
	}
	


	~bool_iterator(){
	}
	


	bool const& operator*(){
		return curr;
	}
	


	bool_iterator& operator++(){
		if(ended){
            return *this;
        }
		if(!initialized){
			initialized = true;
			curr = 0;
		}
		while(!(*vect)[curr]){
			curr ++;
			if(curr > size){
				ended = true;
				return *this;
			}
		}
		cout << curr << endl;
		return *this;
	}



    friend bool operator==(bool_iterator const& lhs, bool_iterator const& rhs){
		if ((lhs.ended && rhs.ended) || (!lhs.initialized && !rhs.initialized)){
			return true;
		}
		return (rhs.curr == lhs.curr) && (lhs.size == rhs.size);
	}
	

	
	friend bool operator!=(bool_iterator const& lhs, bool_iterator const& rhs){
		return !(lhs == rhs);  
	}
	
};

#endif